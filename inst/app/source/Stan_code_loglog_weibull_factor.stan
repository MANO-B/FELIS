data {
  int<lower=0> Median_pos;
  int<lower=0> Median_neg;
  int<lower=0> Nobs_early_pos;
  int<lower=0> Nobs_late_pos;
  int<lower=0> Ncen_early_pos;
  int<lower=0> Ncen_late_pos;
  int<lower=0> Nobs_early_neg;
  int<lower=0> Nobs_late_neg;
  int<lower=0> Ncen_early_neg;
  int<lower=0> Ncen_late_neg;
  int<lower=0> Nexp_pos;
  int<lower=0> Nexp_neg;
  vector[Nexp_pos] ybef_exp_pos;
  vector[Nexp_neg] ybef_exp_neg;
  vector[Nobs_early_pos] yobs_early_pos;
  vector[Nobs_late_pos] yobs_late_pos;
  vector[Ncen_early_pos] ycen_early_pos;
  vector[Ncen_late_pos] ycen_late_pos;
  vector[Nobs_early_neg] yobs_early_neg;
  vector[Nobs_late_neg] yobs_late_neg;
  vector[Ncen_early_neg] ycen_early_neg;
  vector[Ncen_late_neg] ycen_late_neg;
}

parameters {
  // より直感的なパラメータ化
  real<lower=30, upper=3000> mu_base;           // 基準群の尺度パラメータ
  real<lower=0.3, upper=3> beta;               // 後半期間の形状パラメータ
  real<lower=0.3, upper=3> shape_exp;          // 前半期間の形状パラメータ
  real<lower=30, upper=3000> scale_exp_base;   // 前半期間基準尺度パラメータ
  
  // 効果パラメータ（対数スケール）
  real<lower=-2, upper=2> log_effect_pos_early;  // pos群早期効果
  real<lower=-2, upper=2> log_effect_late;       // 後期効果（共通）
  real<lower=-2, upper=2> log_effect_pos_exp;    // pos群前半期間効果
}

transformed parameters {
  // 後半期間の尺度パラメータ
  real<lower=0> mu_early_pos = mu_base * exp(log_effect_pos_early);
  real<lower=0> mu_late_pos = mu_base * exp(log_effect_pos_early + log_effect_late);
  real<lower=0> mu_early_neg = mu_base;
  real<lower=0> mu_late_neg = mu_base * exp(log_effect_late);
  
  // 前半期間の尺度パラメータ
  real<lower=0> scale_exp_pos = scale_exp_base * exp(log_effect_pos_exp);
  real<lower=0> scale_exp_neg = scale_exp_base;
}

model {
  // 事前分布（先に定義）
  mu_base ~ normal(300, 1000) T[30, 3000];
  scale_exp_base ~ normal(300, 1000) T[30, 3000];
  beta ~ normal(0, 5) T[0.3, 3];
  shape_exp ~ normal(0, 5) T[0.3, 3];
  
  log_effect_pos_early ~ normal(0, 5) T[-2, 2];
  log_effect_late ~ normal(0, 5) T[-2, 2];
  log_effect_pos_exp ~ normal(0, 5) T[-2, 2];
  
  // 前半期間の尤度
  ybef_exp_pos ~ weibull(shape_exp, scale_exp_pos);
  ybef_exp_neg ~ weibull(shape_exp, scale_exp_neg);
  
  // 後半期間の尤度
  // 観測されたイベント（死亡）
  target += weibull_lpdf(yobs_early_pos | beta, mu_early_pos);
  target += weibull_lpdf(yobs_late_pos | beta, mu_late_pos);
  target += weibull_lpdf(yobs_early_neg | beta, mu_early_neg);
  target += weibull_lpdf(yobs_late_neg | beta, mu_late_neg);
  
  // 打ち切り（生存）
  target += weibull_lccdf(ycen_early_pos | beta, mu_early_pos);
  target += weibull_lccdf(ycen_late_pos | beta, mu_late_pos);
  target += weibull_lccdf(ycen_early_neg | beta, mu_early_neg);
  target += weibull_lccdf(ycen_late_neg | beta, mu_late_neg);
}

generated quantities {
  array[Nobs_early_pos + Nobs_late_pos + Ncen_early_pos + Ncen_late_pos + Nobs_early_neg + Nobs_late_neg + Ncen_early_neg + Ncen_late_neg] real yhat_total;
  array[Nobs_early_pos + Nobs_late_pos + Ncen_early_pos + Ncen_late_pos + Nobs_early_neg + Nobs_late_neg + Ncen_early_neg + Ncen_late_neg] real yhat_factor_total;
  
  for (i in 1:(Nobs_early_pos + Nobs_late_pos + Ncen_early_pos + Ncen_late_pos + Nobs_early_neg + Nobs_late_neg + Ncen_early_neg + Ncen_late_neg)) {
    real y_bef_neg, y_bef_pos, y_aft_neg, y_aft_pos;
    
    // neg群のシミュレーション
    y_bef_neg = weibull_rng(shape_exp, scale_exp_neg);
    if (y_bef_neg <= Median_neg) {
      y_aft_neg = weibull_rng(beta, mu_early_neg);  // Weibullに統一
    } else {
      y_aft_neg = weibull_rng(beta, mu_late_neg);
    }
    yhat_total[i] = y_bef_neg + y_aft_neg;
    
    // pos群のシミュレーション
    y_bef_pos = weibull_rng(shape_exp, scale_exp_pos);
    if (y_bef_pos <= Median_pos) {
      y_aft_pos = weibull_rng(beta, mu_early_pos);  // Weibullに統一
    } else {
      y_aft_pos = weibull_rng(beta, mu_late_pos);
    }
    yhat_factor_total[i] = y_bef_pos + y_aft_pos;
  }
}
