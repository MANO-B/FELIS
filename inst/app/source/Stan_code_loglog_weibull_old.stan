
    // log-Logistic survival model
    
    functions {
      // Defines the log hazard
      vector log_h (vector t, real shape_tmp, vector scale_tmp) {
        vector[num_elements(t)] log_H;
        for (i in 1:num_elements(t)) {
          log_H[i] = log(shape_tmp)-log(scale_tmp[i])+(shape_tmp-1)*(log(t[i])-log(scale_tmp[i]))-log(1+pow((t[i]/scale_tmp[i]),shape_tmp));
        }
        return log_H;
      }
      
      // Defines the log survival
      vector log_S (vector t, real shape_tmp, vector scale_tmp) {
        vector[num_elements(t)] log_s;
        for (i in 1:num_elements(t)) {
          log_s[i] = -log(1+pow((t[i]/scale_tmp[i]),shape_tmp));
        }
        return log_s;
      }
      
      // Defines the sampling distribution
      real surv_loglogistic_lpdf (vector t, vector d, real shape_tmp, vector scale_tmp) {
        vector[num_elements(t)] log_lik;
        real prob;
        log_lik = d .* log_h(t,shape_tmp,scale_tmp) + log_S(t,shape_tmp,scale_tmp);
        prob = sum(log_lik);
        return prob;
      }
      real log_logistic_rng(real alpha, real beta) {
        real logistic_random;
        logistic_random = logistic_rng(0, 1);
        return alpha * exp(logistic_random / beta);
      }
    }
    
    data {
      int<lower=0> Median;
      int<lower=0> Nobs_early;
      int<lower=0> Nobs_late;
      int<lower=0> Ncen_early;
      int<lower=0> Ncen_late;
      int<lower=0> Nexp;
      int<lower=0> M_bg;
      vector[Nexp] ybef_exp;
      vector[Nobs_early] yobs_early;
      vector[Nobs_late] yobs_late;
      vector[Ncen_early] ycen_early;
      vector[Ncen_late] ycen_late;
    }
    
    transformed data {
      vector[Nobs_early] one_obs_early;
      vector[Nobs_late] one_obs_late;
      vector[Ncen_early] zero_cen_early;
      vector[Ncen_late] zero_cen_late;
    
      for (i in 1:Nobs_early) {
        one_obs_early[i] = 1;
      }
      for (i in 1:Nobs_late) {
        one_obs_late[i] = 1;
      }
      for (i in 1:Ncen_early) {
        zero_cen_early[i] = 0;
      }
      for (i in 1:Ncen_late) {
        zero_cen_late[i] = 0;
      }
    }
    
    parameters {
      //real<lower=3, upper=8> mu;
      //real<lower=3, upper=8> scale_exp;
      real<lower=30, upper=3000> mu;
      real<lower=30, upper=3000> scale_exp;

      real<lower=-2, upper=2> factor;
      real<lower=0.3, upper=3> beta;
      real<lower=0.3, upper=3> shape_exp;
    }
    
    transformed parameters {
      vector[Nobs_early] muobs_early;
      vector[Nobs_late] muobs_late;
      vector[Ncen_early] mucen_early;
      vector[Ncen_late] mucen_late;
      for (i in 1:Nobs_early) {
        muobs_early[i] = mu;
      }
      for (i in 1:Nobs_late) {
        muobs_late[i] = mu * exp(factor);
      }
      for (i in 1:Ncen_early) {
        mucen_early[i] = mu;
      }
      for (i in 1:Ncen_late) {
        mucen_late[i] = mu * exp(factor);
      }
    }
    
    model {
      //ybef_exp ~ weibull(shape_exp, exp(scale_exp));
      ybef_exp ~ weibull(shape_exp, scale_exp);
      
      yobs_early ~ surv_loglogistic(one_obs_early, beta, muobs_early);
      yobs_late ~ surv_loglogistic(one_obs_late, beta, muobs_late);
      ycen_early ~ surv_loglogistic(zero_cen_early, beta, mucen_early);
      ycen_late ~ surv_loglogistic(zero_cen_late, beta, mucen_late);
      
      //mu ~ normal(6,5) T[3, 8];
      //scale_exp ~ normal(6,5) T[3, 8];
      mu ~ normal(300,1000) T[30, 3000];
      scale_exp ~ normal(300,1000) T[30, 3000];
      
      factor ~ normal(0,5) T[-2, 2];
      beta ~ normal(0,5) T[0.3, 3];
      shape_exp ~ normal(0,5) T[0.3, 3];
    }
    
    generated quantities {
      array[Nobs_early + Nobs_late + Ncen_early + Ncen_late] real yhat_exp_uncens;
      array[Nobs_early + Nobs_late + Ncen_early + Ncen_late] real yhat_exp_total;
      array[Nobs_early + Nobs_late + Ncen_early + Ncen_late] real y_exptmp;
      array[Nobs_early + Nobs_late + Ncen_early + Ncen_late] real yhat_early;
      array[Nobs_early + Nobs_late + Ncen_early + Ncen_late] real yhat_late;
      for (i in 1:(Nobs_early + Nobs_late + Ncen_early + Ncen_late)) {
        //y_exptmp[i] = weibull_rng(shape_exp, exp(scale_exp));
        y_exptmp[i] = weibull_rng(shape_exp, scale_exp);
        yhat_early[i] = exp(logistic_rng(log(mu), 1.0/beta));
        //yhat_late[i] = exp(logistic_rng(mu + factor, 1.0/beta));
        yhat_late[i] = exp(logistic_rng(log(mu) + factor, 1.0/beta));
        if (y_exptmp[i] <= Median)
          yhat_exp_uncens[i] = yhat_early[i];
        else
          yhat_exp_uncens[i] = yhat_late[i];
        yhat_exp_total[i] = yhat_exp_uncens[i] + y_exptmp[i];
      }
    }
    
    
