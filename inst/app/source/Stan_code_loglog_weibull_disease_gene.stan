
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
      int<lower=0> Median_pos;
      int<lower=0> Median_neg;
      int<lower=0> Nobs;
      int<lower=0> Ncen;
      int<lower=0> Ntot;
      int<lower=0> Nexp;
      int<lower=0> M_bg;
      int n_disease;
      vector[Nexp] ybef_exp;
      vector[Nexp] xfactor;
      matrix[Nexp, n_disease] xdisease;
      vector[Nobs] yobs;
      vector[Ncen] ycen;
      matrix[Nobs, M_bg] Xobs_bg;
      matrix[Ncen, M_bg] Xcen_bg;
    }
    
    transformed data {
      vector[Nobs] one_obs;
      vector[Ncen] zero_cen;
    
      matrix[Nobs, M_bg] Xobs_bg_early_pos = Xobs_bg;
      matrix[Nobs, M_bg] Xobs_bg_late_pos = Xobs_bg;
      matrix[Ncen, M_bg] Xcen_bg_early_pos = Xcen_bg;
      matrix[Ncen, M_bg] Xcen_bg_late_pos = Xcen_bg;
      matrix[Nobs, M_bg] Xobs_bg_early_neg = Xobs_bg;
      matrix[Nobs, M_bg] Xobs_bg_late_neg = Xobs_bg;
      matrix[Ncen, M_bg] Xcen_bg_early_neg = Xcen_bg;
      matrix[Ncen, M_bg] Xcen_bg_late_neg = Xcen_bg;
    
      for (i in 1:Nobs) {
        Xobs_bg_early_pos[i,1] = 0;
        Xobs_bg_late_pos[i,1] = 1;
        Xobs_bg_early_neg[i,1] = 0;
        Xobs_bg_late_neg[i,1] = 1;
        Xobs_bg_early_pos[i,2] = 1;
        Xobs_bg_late_pos[i,2] = 1;
        Xobs_bg_early_neg[i,2] = 0;
        Xobs_bg_late_neg[i,2] = 0;
        one_obs[i] = 1;
      }
      for (i in 1:Ncen) {
        Xcen_bg_early_pos[i,1] = 0;
        Xcen_bg_late_pos[i,1] = 1;
        Xcen_bg_early_neg[i,1] = 0;
        Xcen_bg_late_neg[i,1] = 1;
        Xcen_bg_early_pos[i,2] = 1;
        Xcen_bg_late_pos[i,2] = 1;
        Xcen_bg_early_neg[i,2] = 0;
        Xcen_bg_late_neg[i,2] = 0;
        zero_cen[i] = 0;
      }
    }
    
    parameters {
      vector<lower=-3, upper=3>[M_bg] alpha;
      real<lower=0.1, upper=10> beta;
      real<lower=-3, upper=3> factor;
      vector<lower=-3, upper=3>[n_disease] disease;
      real<lower=3, upper=10> mu;
      real<lower=0.1, upper=10> shape_exp;
      real<lower=3, upper=10> scale_exp;
    }
    
    transformed parameters {
      vector[Nobs] linpredobs;
      vector[Ncen] linpredcen;
      vector[Nobs] muobs;
      vector[Ncen] mucen;
      vector[Nexp] disease_effect;
    
      linpredobs = Xobs_bg*alpha;
      for (i in 1:Nobs) {
        muobs[i] = exp(mu + linpredobs[i]);
      }
      linpredcen = Xcen_bg*alpha;
      for (i in 1:Ncen) {
        mucen[i] = exp(mu + linpredcen[i]);
      }
      disease_effect = xdisease*disease;
    }
    
    model {
      ybef_exp ~ weibull(shape_exp, exp(scale_exp + xfactor * factor + disease_effect));
      yobs ~ surv_loglogistic(one_obs, beta, muobs);
      ycen ~ surv_loglogistic(zero_cen, beta, mucen);
    
      mu ~ normal(6,2);
      alpha ~ normal(0,2);
      beta ~ normal(0,2);
      shape_exp ~ normal(1.0,1.0);
      scale_exp ~ normal(5.0,2.0);
      factor ~ normal(0.0,1.0);
      disease_effect ~ normal(0.0,1.0);
    }
    
    generated quantities {
      real yhat_total[Nobs + Ncen];
      real yhat_factor_total[Nobs + Ncen];
      real y_bef;
      real y_aft;
      for (i in 1:Nobs) {
        y_bef = weibull_rng(shape_exp, exp(scale_exp + disease_effect[i]));
        if (y_bef <= Median_neg) {
          y_aft = exp(logistic_rng(((mu + Xobs_bg_early_neg[i,] * alpha)), 1.0/beta));
          yhat_total[i] = y_aft + y_bef;
        } else {
          y_aft = exp(logistic_rng(((mu + Xobs_bg_late_neg[i,] * alpha)), 1.0/beta));
          yhat_total[i] = y_aft + y_bef;
        }
        y_bef =  weibull_rng(shape_exp, exp(scale_exp + factor + disease_effect[i]));
        if (y_bef <= Median_pos) {
          y_aft = exp(logistic_rng(((mu + Xobs_bg_early_pos[i,] * alpha)), 1.0/beta));
          yhat_factor_total[i] = y_aft + y_bef;
        } else {
          y_aft = exp(logistic_rng(((mu + Xobs_bg_late_pos[i,] * alpha)), 1.0/beta));
          yhat_factor_total[i] = y_aft + y_bef;
        }
      }
      for (i in 1:Ncen) {
        y_bef = weibull_rng(shape_exp, exp(scale_exp + disease_effect[Nobs + i]));
        if (y_bef <= Median_neg) {
          y_aft = exp(logistic_rng(((mu + Xcen_bg_early_neg[i,] * alpha)), 1.0/beta));
          yhat_total[Nobs + i] = y_aft + y_bef;
        } else {
          y_aft = exp(logistic_rng(((mu + Xcen_bg_late_neg[i,] * alpha)), 1.0/beta));
          yhat_total[Nobs + i] = y_aft + y_bef;
        }
        y_bef =  weibull_rng(shape_exp, exp(scale_exp + factor + disease_effect[Nobs + i]));
        if (y_bef <= Median_pos) {
          y_aft = exp(logistic_rng(((mu + Xcen_bg_early_pos[i,] * alpha)), 1.0/beta));
          yhat_factor_total[Nobs + i] = y_aft + y_bef;
        } else {
          y_aft = exp(logistic_rng(((mu + Xcen_bg_late_pos[i,] * alpha)), 1.0/beta));
          yhat_factor_total[Nobs + i] = y_aft + y_bef;
        }
      }
    }
    
    
