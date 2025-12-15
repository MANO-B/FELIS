data {
  int<lower=0> Median;         // 前半期間の中央値（群分けの基準）
  int<lower=0> Nobs_early;     // 後半期間の観測死亡数（前半早期群）
  int<lower=0> Nobs_late;      // 後半期間の観測死亡数（前半後期群）
  int<lower=0> Ncen_early;     // 後半期間の打ち切り数（前半早期群）
  int<lower=0> Ncen_late;      // 後半期間の打ち切り数（前半後期群）
  int<lower=0> Nexp;           // 前半期間の症例数
  vector[Nexp] ybef_exp;       // 前半の生存期間（全例観測済み）
  vector[Nobs_early] yobs_early; // 後半期間死亡時間（前半早期群）
  vector[Nobs_late] yobs_late;   // 後半期間死亡時間（前半後期群）
  vector[Ncen_early] ycen_early; // 後半期間打ち切り時間（前半早期群）
  vector[Ncen_late] ycen_late;   // 後半期間打ち切り時間（前半後期群）
}

parameters {
  // Weibull分布パラメータ（前半期間）
  real<lower=0.3, upper=3> shape_exp;      // 形状パラメータ
  real<lower=30, upper=3000> scale_exp;    // 尺度パラメータ
  
  // Log-logistic分布パラメータ（後半期間：前半早期群）
  real<lower=30, upper=3000> mu;           // 基準群の尺度パラメータ
  real<lower=0.3, upper=3> beta;          // 基準群の形状パラメータ
  
  // 群間差のfactor
  real<lower=-2, upper=2> factor;         // 前半後期群の尺度パラメータ差（対数スケール）
}

transformed parameters {
  // 前半後期群の尺度パラメータ
  real<lower=0> mu_late = mu * exp(factor);
}

model {
  // 事前分布
  shape_exp ~ normal(0, 5) T[0.3, 3];
  scale_exp ~ normal(300, 1000) T[30, 3000];
  mu ~ normal(300, 1000) T[30, 3000];
  beta ~ normal(0, 5) T[0.3, 3];
  factor ~ normal(0, 5) T[-2, 2];
  
  // 尤度
  // 前半期間: Weibull分布（打ち切りなし）
  ybef_exp ~ weibull(shape_exp, scale_exp);
  
  // 後半期間: Log-logistic分布（最新Stan）
  // 前半早期群（Median以下）
  for (i in 1:Nobs_early) {
    target += weibull_lpdf(yobs_early[i] | beta, mu);
  }
  for (i in 1:Ncen_early) {
    target += weibull_lccdf(ycen_early[i] | beta, mu);
  }
  
  // 前半後期群（Median超過）
  for (i in 1:Nobs_late) {
    target += weibull_lpdf(yobs_late[i] | beta, mu_late);
  }
  for (i in 1:Ncen_late) {
    target += weibull_lccdf(ycen_late[i] | beta, mu_late);
  }
}

generated quantities {
  // 予測値生成
  array[Nobs_early + Nobs_late + Ncen_early + Ncen_late] real yhat_exp_uncens;
  array[Nobs_early + Nobs_late + Ncen_early + Ncen_late] real yhat_exp_total;
  array[Nobs_early + Nobs_late + Ncen_early + Ncen_late] real y_exptmp;
  array[Nobs_early + Nobs_late + Ncen_early + Ncen_late] real yhat_early;
  array[Nobs_early + Nobs_late + Ncen_early + Ncen_late] real yhat_late;
  
  // 群間比較の要約統計
  real hazard_ratio = exp(factor);      // 前半後期群 vs 前半早期群のハザード比
  real median_survival_early = mu;      // 前半早期群の後半期間中央値
  real median_survival_late = mu_late;  // 前半後期群の後半期間中央値
  
  // 予測値の生成
  for (i in 1:(Nobs_early + Nobs_late + Ncen_early + Ncen_late)) {
    // 前半期間の予測
    y_exptmp[i] = weibull_rng(shape_exp, scale_exp);
    
    // 群分類と後半期間の予測
    yhat_early[i] = weibull_rng(beta, mu);
    yhat_late[i] = weibull_rng(beta, mu_late);
    
    // 前半期間に基づく群分けで後半期間を決定
    if (y_exptmp[i] <= Median) {
      yhat_exp_uncens[i] = yhat_early[i];
    } else {
      yhat_exp_uncens[i] = yhat_late[i];
    }
    
    // 総生存時間
    yhat_exp_total[i] = yhat_exp_uncens[i] + y_exptmp[i];
  }
}
