data {
  int<lower=0> n_ids;
  int<lower=0> n_dyads;
  int<lower=0> n_beh;
  array[n_beh + 1] int<lower=0> behav_types; // 1: counts (poiss), 2: counts (proportion, binomial)
  array[n_dyads, n_beh] int<lower=0> interactions;
  matrix[n_dyads, n_beh] interactions_cont;
  array[n_dyads, 2] int<lower=0> dyads_navi;
  matrix[n_dyads, n_beh] obseff;
  array[n_dyads, n_beh] int<lower=0> obseff_int;

  // priors for intercepts
  matrix[n_beh, 2] prior_matrix;

  // info about indexing and number of optional shape parameter(s)
  int<lower=0> gamma_shape_n;
  array[n_beh] int<lower=0> gamma_shape_pos;
  int<lower=0> beta_shape_n;
  array[n_beh] int<lower=0> beta_shape_pos;

  // additional predictors
  vector[n_ids] indi_cat_pred; //<lower=0,upper=1>
  vector[n_ids] indi_covariate_pred;
  vector[n_dyads] dyad_cat_pred; //<lower=0,upper=1>
  vector[n_dyads] dyad_covariate_pred;
  // and flags
  int do_indi_cat;
  int do_indi_covariate;
  int do_dyad_cat;
  int do_dyad_covariate;
}

parameters {
  real<lower=0> indi_soc_sd; // SD parameter
  vector[n_ids] indi_soc_vals_z; // blups on z-scale
  real<lower=0> dyad_soc_sd; // SD parameter
  vector[n_dyads] dyad_soc_vals_z; // blups on z-scale
  vector[n_beh] beh_intercepts; // intercepts for behaviours

  vector<lower=0>[gamma_shape_n] shapes;
  vector<lower=0>[beta_shape_n] shapes_beta;

  // real indi_intercept;
  real indi_cat_eff;
  real indi_covariate_eff;
  // real dyad_intercept;
  real dyad_cat_eff;
  real dyad_covariate_eff;
}

transformed parameters {
  vector[n_ids] indi_soc_vals;  // actual blups
  vector[n_dyads] dyad_soc_vals;  // actual blups
  vector[n_ids] indi_aux = rep_vector(0.0, n_ids);
  vector[n_dyads] dyad_aux = rep_vector(0.0, n_dyads);
  vector[n_dyads] scaled_indi_sums;  // sums of dyadic values, scaled

  if (do_indi_cat == 1) {
    indi_aux = indi_aux + indi_cat_eff * indi_cat_pred;
  }
  if (do_indi_covariate == 1) {
    indi_aux = indi_aux + indi_covariate_eff * indi_covariate_pred;
  }
  indi_soc_vals = (indi_soc_sd * indi_soc_vals_z) + indi_aux;


  if (do_dyad_cat == 1) {
    dyad_aux = dyad_aux + dyad_cat_eff * dyad_cat_pred;
  }
  if (do_dyad_covariate == 1) {
    dyad_aux = dyad_aux + dyad_covariate_eff * dyad_covariate_pred;
  }
  dyad_soc_vals = (dyad_soc_sd * dyad_soc_vals_z) + dyad_aux;

  scaled_indi_sums = sqrt(0.5) * (indi_soc_vals[dyads_navi[, 1]] + indi_soc_vals[dyads_navi[, 2]]);
}

model {
  // linear predictor
  vector[n_dyads] lp = scaled_indi_sums + dyad_soc_vals;


  for (i in 1:n_beh) {
    if (behav_types[i] == 1) {
      interactions[, i] ~ poisson(exp(lp + beh_intercepts[i] + log(obseff[, i])));
    }
    if (behav_types[i] == 2) {
      interactions[, i] ~ binomial_logit(obseff_int[, i], lp + beh_intercepts[i]);
    }
    if (behav_types[i] == 3) {
      interactions_cont[, i] ~ gamma(shapes[gamma_shape_pos[i]],
                                     shapes[gamma_shape_pos[i]] * exp(-(lp + beh_intercepts[i] + log(obseff[, i]))))  ; //
    }
    if (behav_types[i] == 4) {
      interactions_cont[, i] ~ beta(inv_logit(lp + beh_intercepts[i]) * shapes_beta[beta_shape_pos[i]],
                                     (1 - inv_logit(lp + beh_intercepts[i])) *  shapes_beta[beta_shape_pos[i]] )  ; //
    }
  }

  for (i in 1:n_beh) {
    // declare priors for intercepts based on observed data
    // has been done outside Stan (and is passed as data)
    beh_intercepts[i] ~ student_t(3, prior_matrix[i, 1], prior_matrix[i, 2]);

    // extra priors for shape/dispersion parameters
    if (behav_types[i] == 3) {
      shapes[gamma_shape_pos[i]] ~ gamma(0.01, 0.01);
    }
    if (behav_types[i] == 4) {
      shapes_beta[beta_shape_pos[i]] ~ gamma(0.1, 0.1);
    }
  }
  indi_soc_vals_z ~ normal(0, 1);
  dyad_soc_vals_z ~ normal(0, 1);
  // indi_soc_sd ~ student_t(3, 0, 1);
  // dyad_soc_sd ~ student_t(3, 0, 1);
  indi_soc_sd ~ exponential(2);
  dyad_soc_sd ~ exponential(2);
  indi_cat_eff ~ normal(0, 0.5);
  indi_covariate_eff ~ normal(0, 0.5);
  dyad_cat_eff ~ normal(0, 0.5);
  dyad_covariate_eff ~ normal(0, 0.5);

  // indi_intercept ~ normal(0, 1);
  // dyad_intercept ~ normal(0, 1);

}

generated quantities {
  array[n_dyads, n_beh] int interactions_pred;
  matrix[n_dyads, n_beh] interactions_pred_cont;
  matrix[n_dyads, n_beh] log_lik; // pointwise log likelihoods

  {
    vector[n_dyads] LP = scaled_indi_sums + dyad_soc_vals; // linear predictor without behaviour-specific intercept
    for (i in 1:n_beh) {
      if (behav_types[i] == 1) {
        interactions_pred[, i] = poisson_log_rng(LP + beh_intercepts[i] + log(obseff[, i]));
        for (jj in 1:n_dyads) { // loglik
          log_lik[jj, i] = poisson_log_lpmf(interactions[jj, i] | LP[jj] + beh_intercepts[i] + log(obseff[jj, i]));
        }
      }
      if (behav_types[i] == 2) {
        interactions_pred[, i] = binomial_rng(obseff_int[, i], inv_logit(LP + beh_intercepts[i]));
        for (jj in 1:n_dyads) { // loglik
          log_lik[jj, i] = binomial_logit_lpmf(interactions[jj, i] | obseff_int[jj, i], LP[jj] + beh_intercepts[i]);
        }
      }
      if (behav_types[i] == 3) {
        interactions_pred_cont[, i] = to_vector(gamma_rng(shapes[gamma_shape_pos[i]],
                                      shapes[gamma_shape_pos[i]] *
                                        exp(-(LP + beh_intercepts[i] + log(obseff[, i])))));
        for (jj in 1:n_dyads) { // loglik
          log_lik[jj, i] = gamma_lpdf(interactions_cont[jj, i] | shapes[gamma_shape_pos[i]],
                                       shapes[gamma_shape_pos[i]] * exp(-(LP[jj] + beh_intercepts[i] + log(obseff[jj, i]))));
        }
      }
      if (behav_types[i] == 4) {
        interactions_pred_cont[, i] = to_vector(beta_rng(
                                          shapes_beta[beta_shape_pos[i]] * (    inv_logit(LP + beh_intercepts[i])),
                                          shapes_beta[beta_shape_pos[i]] * (1 - inv_logit(LP + beh_intercepts[i]))
                                          ));
        for (jj in 1:n_dyads) { // loglik
          log_lik[jj, i] = beta_lpdf(interactions_cont[jj, i] |
                                     inv_logit(LP[jj] + beh_intercepts[i])  * shapes_beta[beta_shape_pos[i]],
                                (1 - inv_logit(LP[jj] + beh_intercepts[i])) * shapes_beta[beta_shape_pos[i]]);

        }
      }
    }
  }
}


