// interactions can be any combination of counts,
//   discrete proportions (count and trials), continuous proportions (beta),
//   positive continuous (gamma)

// we require actually a second observation effort matrix that is
//   of type 'integer', otherwise stan will complain

// similarly, there are two matrices for the actual interactions,
//   one with integers (for pois and binomial) and
//   one with continuous numbers (for beta and gamma models)

functions {
  real lin2prob(real x, real obseff) {
    return(1 - exp(-exp(x) * obseff));
  }

  real prob2lin(real x, real obseff) {
    return(log((-log(1 - x))/obseff));
  }

  real gamma0_lpdf(real y, real alpha, real beta, real hu, real obseff) {
    // real bern_mu = 1.0 - exp(-(exp(hu) * obseff)); //
    real bern_mu = lin2prob(hu, obseff);
    if (y == 0) {
      return bernoulli_lpmf(0 | bern_mu);
    } else {
      return bernoulli_lpmf(1 | bern_mu) + gamma_lpdf(y | alpha, beta);
    }
  }

  real gamma0_rng(real alpha, real beta, real hu, real obseff) {
    real out = 0.0;
    int aux = binomial_rng(1, lin2prob(hu, obseff));
    if (aux == 1) {
      out = gamma_rng(alpha, beta);
    }
    return out;
  }
}

data {
  int<lower=0> n_ids;
  int<lower=0> n_dyads;
  int<lower=0> n_beh;
  // int<lower=0> behav_types[n_beh + 1]; // 1: counts (poiss), 2: counts (proportion, binomial) // old way
  array[n_beh + 1] int<lower=0> behav_types; // 1: counts (poiss), 2: counts (proportion, binomial)
  array[n_dyads, n_beh] int<lower=0> interactions;
  matrix[n_dyads, n_beh] interactions_cont;
  array[n_dyads, 2] int<lower=0> dyads_navi;
  matrix[n_dyads, n_beh] obseff;
  array[n_dyads, n_beh] int<lower=0> obseff_int;

  // flag for prior predictive check
  int<lower=0,upper=1> prior_only;

  // priors for intercepts
  matrix[n_beh, 2] prior_matrix;
  matrix[n_beh, 2] prior_matrix2;

  // info about indexing and number of optional shape parameter(s)
  int<lower=0> gamma_shape_n;
  array[n_beh] int<lower=0> gamma_shape_pos;
  int<lower=0> beta_shape_n;
  array[n_beh] int<lower=0> beta_shape_pos;
}
parameters {
  real<lower=0> indi_soc_sd; // SD parameter
  vector[n_ids] indi_soc_vals_z; // blups on z-scale
  real<lower=0> dyad_soc_sd; // SD parameter
  vector[n_dyads] dyad_soc_vals_z; // blups on z-scale
  vector[n_beh] beh_intercepts; // intercepts for behaviours
  vector[gamma_shape_n] beh_intercepts_add; // second set of intercepts for behaviours

  vector<lower=0>[gamma_shape_n] shapes_gamma_raw; // unconstrained
  vector<lower=0>[beta_shape_n] shapes_beta;
}
transformed parameters {
  vector[n_ids] indi_soc_vals;  // actual blups
  vector[n_dyads] dyad_soc_vals;  // actual blups
  vector[n_dyads] scaled_indi_sums;  // sums of dyadic values, scaled
  vector<lower=0>[gamma_shape_n] shapes_gamma = exp(shapes_gamma_raw);
  indi_soc_vals = (indi_soc_sd * indi_soc_vals_z);
  dyad_soc_vals = (dyad_soc_sd * dyad_soc_vals_z);
  scaled_indi_sums = sqrt(0.5) * (indi_soc_vals[dyads_navi[, 1]] + indi_soc_vals[dyads_navi[, 2]]);
}
model {
  vector[n_dyads] lp = rep_vector(0.0, n_dyads);
  // linear predictor
  lp = lp + scaled_indi_sums + dyad_soc_vals;

  if (prior_only == 0) {
    for (i in 1:n_beh) {
      // first attempt at gamma0
      if (behav_types[i] == 5) {
        for (n in 1:n_dyads) {
          interactions_cont[n, i] ~ gamma0(
            // alpha (shape):
            shapes_gamma[gamma_shape_pos[i]],
            // beta (shape/mu):
            shapes_gamma[gamma_shape_pos[i]] / exp(lp[n] + beh_intercepts_add[gamma_shape_pos[i]] + log(obseff[n, i])),
            // bernoulli intercept:
            beh_intercepts[i],
            // obseff for bernoulli offset
            obseff[n, i]
            );
        }
      }
      if (behav_types[i] == 1) {
        interactions[, i] ~ poisson(exp(lp + beh_intercepts[i] + log(obseff[, i])));
      }
      if (behav_types[i] == 2) {
        interactions[, i] ~ binomial_logit(obseff_int[, i], lp + beh_intercepts[i]);
      }
      if (behav_types[i] == 3) {
        interactions_cont[, i] ~ gamma(shapes_gamma[gamma_shape_pos[i]],
                                       shapes_gamma[gamma_shape_pos[i]]  / exp(lp + beh_intercepts[i] + log(obseff[, i])))  ; //
      }
      if (behav_types[i] == 4) {
        interactions_cont[, i] ~ beta(inv_logit(lp + beh_intercepts[i]) * shapes_beta[beta_shape_pos[i]],
                                     (1 - inv_logit(lp + beh_intercepts[i])) *  shapes_beta[beta_shape_pos[i]] ); //
      }
    }
  }


  for (i in 1:n_beh) {
    // declare priors for intercepts based on observed data
    // has been done outside Stan (and is passed as data)
    beh_intercepts[i] ~ student_t(3, prior_matrix[i, 1], prior_matrix[i, 2]);

    // extra priors for shape/dispersion parameters
    if (behav_types[i] == 3) {
      // shapes[gamma_shape_pos[i]] ~ gamma(0.01, 0.01);
      shapes_gamma_raw[gamma_shape_pos[i]] ~ normal(0, 1);
    }
    if (behav_types[i] == 4) {
      shapes_beta[beta_shape_pos[i]] ~ gamma(0.1, 0.1);
    }
    if (behav_types[i] == 5) {
      // gamma prior comes from prior_matrix2
      beh_intercepts_add[gamma_shape_pos[i]] ~ student_t(3, prior_matrix2[i, 1], prior_matrix2[i, 2]); // gamma in mixture
      shapes_gamma_raw[gamma_shape_pos[i]] ~ normal(0, 1);
    }
  }
  indi_soc_vals_z ~ normal(0, 1);
  dyad_soc_vals_z ~ normal(0, 1);
  // indi_soc_sd ~ student_t(3, 0, 1);
  // dyad_soc_sd ~ student_t(3, 0, 1);
  indi_soc_sd ~ exponential(2);
  dyad_soc_sd ~ exponential(2);
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
        // interactions_pred_cont[, i] = to_vector(poisson_log_rng(LP + beh_intercepts[i] + log(obseff[, i])));
        for (jj in 1:n_dyads) { // loglik
          log_lik[jj, i] = poisson_log_lpmf(interactions[jj, i] | LP[jj] + beh_intercepts[i] + log(obseff[jj, i]));
        }
      }
      if (behav_types[i] == 2) {
        interactions_pred[, i] = binomial_rng(obseff_int[, i], inv_logit(LP + beh_intercepts[i]));
        // interactions_pred_cont[, i] = to_vector(binomial_rng(obseff_int[, i], inv_logit(LP + beh_intercepts[i])));
        for (jj in 1:n_dyads) { // loglik
          log_lik[jj, i] = binomial_logit_lpmf(interactions[jj, i] | obseff_int[jj, i], LP[jj] + beh_intercepts[i]);
        }
      }
      if (behav_types[i] == 3) {
        interactions_pred_cont[, i] = to_vector(gamma_rng(shapes_gamma[gamma_shape_pos[i]],
        shapes_gamma[gamma_shape_pos[i]] *
        exp(-(LP + beh_intercepts[i] + log(obseff[, i])))));
        for (jj in 1:n_dyads) { // loglik
          log_lik[jj, i] = gamma_lpdf(interactions_cont[jj, i] | shapes_gamma[gamma_shape_pos[i]],
          shapes_gamma[gamma_shape_pos[i]] * exp(-(LP[jj] + beh_intercepts[i] + log(obseff[jj, i]))));
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
      if (behav_types[i] == 5) {
        for (n in 1:n_dyads) {
          interactions_pred_cont[n, i] = gamma0_rng(
            // alpha (shape):
            shapes_gamma[gamma_shape_pos[i]],
            // beta (shape/mu):
            shapes_gamma[gamma_shape_pos[i]] / exp(LP[n] + beh_intercepts_add[gamma_shape_pos[i]] + log(obseff[n, i])),
            // bernoulli intercept:
            beh_intercepts[i],
            // obseff for bernoulli offset
            obseff[n, i]
            );
        }
      }
    }

  }

}


