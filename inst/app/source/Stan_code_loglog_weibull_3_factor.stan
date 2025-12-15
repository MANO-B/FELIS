
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
    }
    
    data {
      int<lower=0> Median_2;
      int<lower=0> Median_1;
      int<lower=0> Median_0;
      int<lower=0> Nobs;
      int<lower=0> Ncen;
      int<lower=0> Ntot;
      int<lower=0> Nexp;
      int<lower=0> M_bg;
      vector[Nexp] ybef_exp;
      vector[Nexp] xfactor;
      vector[Nexp] xfactor2;
      vector[Nobs] yobs;
      vector[Ncen] ycen;
      matrix[Nobs, M_bg] Xobs_bg;
      matrix[Ncen, M_bg] Xcen_bg;
    }
    
    transformed data {
      vector[Nobs] one_obs;
      vector[Ncen] zero_cen;
    
      matrix[Nobs, M_bg] Xobs_bg_early_2 = Xobs_bg;
      matrix[Nobs, M_bg] Xobs_bg_late_2 = Xobs_bg;
      matrix[Ncen, M_bg] Xcen_bg_early_2 = Xcen_bg;
      matrix[Ncen, M_bg] Xcen_bg_late_2 = Xcen_bg;
      matrix[Nobs, M_bg] Xobs_bg_early_1 = Xobs_bg;
      matrix[Nobs, M_bg] Xobs_bg_late_1 = Xobs_bg;
      matrix[Ncen, M_bg] Xcen_bg_early_1 = Xcen_bg;
      matrix[Ncen, M_bg] Xcen_bg_late_1 = Xcen_bg;
      matrix[Nobs, M_bg] Xobs_bg_early_0 = Xobs_bg;
      matrix[Nobs, M_bg] Xobs_bg_late_0 = Xobs_bg;
      matrix[Ncen, M_bg] Xcen_bg_early_0 = Xcen_bg;
      matrix[Ncen, M_bg] Xcen_bg_late_0 = Xcen_bg;
    
      for (i in 1:Nobs) {
        Xobs_bg_early_2[i,1] = 0;
        Xobs_bg_late_2[i,1] = 1;
        Xobs_bg_early_1[i,1] = 0;
        Xobs_bg_late_1[i,1] = 1;
        Xobs_bg_early_0[i,1] = 0;
        Xobs_bg_late_0[i,1] = 1;
        Xobs_bg_early_2[i,2] = 0;
        Xobs_bg_late_2[i,2] = 0;
        Xobs_bg_early_1[i,2] = 1;
        Xobs_bg_late_1[i,2] = 1;
        Xobs_bg_early_0[i,2] = 0;
        Xobs_bg_late_0[i,2] = 0;
        Xobs_bg_early_2[i,3] = 1;
        Xobs_bg_late_2[i,3] = 1;
        Xobs_bg_early_1[i,3] = 0;
        Xobs_bg_late_1[i,3] = 0;
        Xobs_bg_early_0[i,3] = 0;
        Xobs_bg_late_0[i,3] = 0;
        one_obs[i] = 1;
      }
      for (i in 1:Ncen) {
        Xcen_bg_early_2[i,1] = 0;
        Xcen_bg_late_2[i,1] = 1;
        Xcen_bg_early_1[i,1] = 0;
        Xcen_bg_late_1[i,1] = 1;
        Xcen_bg_early_0[i,1] = 0;
        Xcen_bg_late_0[i,1] = 1;
        Xcen_bg_early_2[i,2] = 0;
        Xcen_bg_late_2[i,2] = 0;
        Xcen_bg_early_1[i,2] = 1;
        Xcen_bg_late_1[i,2] = 1;
        Xcen_bg_early_0[i,2] = 0;
        Xcen_bg_late_0[i,2] = 0;
        Xcen_bg_early_2[i,3] = 1;
        Xcen_bg_late_2[i,3] = 1;
        Xcen_bg_early_1[i,3] = 0;
        Xcen_bg_late_1[i,3] = 0;
        Xcen_bg_early_0[i,3] = 0;
        Xcen_bg_late_0[i,3] = 0;
        zero_cen[i] = 0;
      }
    }
    
    parameters {
      vector<lower=-3, upper=3>[M_bg] alpha;
      real<lower=0.1, upper=10> beta;
      real<lower=-3, upper=3> factor;
      real<lower=-3, upper=3> factor2;
      real<lower=3, upper=10> mu;
      real<lower=0.1, upper=10> shape_exp;
      real<lower=3, upper=10> scale_exp;
    }
    
    transformed parameters {
      vector[Nobs] linpredobs;
      vector[Ncen] linpredcen;
      vector[Nobs] muobs;
      vector[Ncen] mucen;
      linpredobs = Xobs_bg*alpha;
      for (i in 1:Nobs) {
        muobs[i] = exp(mu + linpredobs[i]);
      }
      linpredcen = Xcen_bg*alpha;
      for (i in 1:Ncen) {
        mucen[i] = exp(mu + linpredcen[i]);
      }
    }
    
    model {
      ybef_exp ~ weibull(shape_exp, exp(scale_exp + xfactor * factor + xfactor2 * factor2));
      yobs ~ surv_loglogistic(one_obs, beta, muobs);
      ycen ~ surv_loglogistic(zero_cen, beta, mucen);
    
      mu ~ normal(6,2);
      alpha ~ normal(0,2);
      beta ~ normal(0,2);
      shape_exp ~ normal(1.0,1.0);
      scale_exp ~ normal(5.0,2.0);
      factor ~ normal(0.0,1.0);
      factor2 ~ normal(0.0,1.0);
    }
    
    generated quantities {
      real yhat_total[Nobs + Ncen];
      real yhat_factor_total[Nobs + Ncen];
      real yhat_factor2_total[Nobs + Ncen];
      real y_bef;
      real y_aft;
      for (i in 1:Nobs) {
        y_bef = weibull_rng(shape_exp, exp(scale_exp));
        if (y_bef <= Median_0) {
          y_aft = exp(logistic_rng(((mu + Xobs_bg_early_0[i,] * alpha)), 1.0/beta));
          yhat_total[i] = y_aft + y_bef;
        } else {
          y_aft = exp(logistic_rng(((mu + Xobs_bg_late_0[i,] * alpha)), 1.0/beta));
          yhat_total[i] = y_aft + y_bef;
        }
        y_bef =  weibull_rng(shape_exp, exp(scale_exp + factor));
        if (y_bef <= Median_1) {
          y_aft = exp(logistic_rng(((mu + Xobs_bg_early_1[i,] * alpha)), 1.0/beta));
          yhat_factor_total[i] = y_aft + y_bef;
        } else {
          y_aft = exp(logistic_rng(((mu + Xobs_bg_late_1[i,] * alpha)), 1.0/beta));
          yhat_factor_total[i] = y_aft + y_bef;
        }
        y_bef =  weibull_rng(shape_exp, exp(scale_exp + factor2));
        if (y_bef <= Median_2) {
          y_aft = exp(logistic_rng(((mu + Xobs_bg_early_2[i,] * alpha)), 1.0/beta));
          yhat_factor2_total[i] = y_aft + y_bef;
        } else {
          y_aft = exp(logistic_rng(((mu + Xobs_bg_late_2[i,] * alpha)), 1.0/beta));
          yhat_factor2_total[i] = y_aft + y_bef;
        }
      }
      for (i in 1:Ncen) {
        y_bef = weibull_rng(shape_exp, exp(scale_exp));
        if (y_bef <= Median_0) {
          y_aft = exp(logistic_rng(((mu + Xcen_bg_early_0[i,] * alpha)), 1.0/beta));
          yhat_total[Nobs + i] = y_aft + y_bef;
        } else {
          y_aft = exp(logistic_rng(((mu + Xcen_bg_late_0[i,] * alpha)), 1.0/beta));
          yhat_total[Nobs + i] = y_aft + y_bef;
        }
        y_bef =  weibull_rng(shape_exp, exp(scale_exp + factor));
        if (y_bef <= Median_1) {
          y_aft = exp(logistic_rng(((mu + Xcen_bg_early_1[i,] * alpha)), 1.0/beta));
          yhat_factor_total[Nobs + i] = y_aft + y_bef;
        } else {
          y_aft = exp(logistic_rng(((mu + Xcen_bg_late_1[i,] * alpha)), 1.0/beta));
          yhat_factor_total[Nobs + i] = y_aft + y_bef;
        }
        y_bef =  weibull_rng(shape_exp, exp(scale_exp + factor2));
        if (y_bef <= Median_2) {
          y_aft = exp(logistic_rng(((mu + Xcen_bg_early_2[i,] * alpha)), 1.0/beta));
          yhat_factor2_total[Nobs + i] = y_aft + y_bef;
        } else {
          y_aft = exp(logistic_rng(((mu + Xcen_bg_late_2[i,] * alpha)), 1.0/beta));
          yhat_factor2_total[Nobs + i] = y_aft + y_bef;
        }
      }
    }
    
    
