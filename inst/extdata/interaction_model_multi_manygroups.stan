// model for multiple groups with modelled group-level pars
//   (behavior baseline, indi SD, dyad SD)
// hard constraint: one behavior only!!! and must be the same for all periods
// hard constraint 2: frequencies only (count or prop)
// only tested with count!

data {
  int<lower=1> n_pwise_ids; // unique period-ids
  int<lower=1> n_pwise_dyads;
  int<lower=0> n_periods;
  int<lower=1,upper=2> behav_types;
  int<lower=1,upper=1> n_beh;
  // 1: counts (poiss); 2: counts (proportion, binomial);
  array[n_pwise_dyads, 1] int<lower=0> interactions;
  matrix[n_pwise_dyads, 1] obseff; // for count
  array[n_pwise_dyads, 1] int<lower=0> obseff_int; // for prop

  array[n_pwise_dyads] int<lower=0> index_period; // for each dyad: which period does it occur in
  array[n_pwise_ids] int<lower=0> index_period_individual; // for each ID: which period does it occur in

  array[n_pwise_dyads] int<lower=0> index_dyad;
  array[n_pwise_dyads] int<lower=0> index_id1;
  array[n_pwise_dyads] int<lower=0> index_id2;

  array[n_periods] int n_dyads_perperiod;

  // priors for intercept
  matrix[n_periods, 2] prior_matrix;
}
transformed data {
}
parameters {
  vector<lower=0>[n_periods] indi_soc_sd;
  vector[n_pwise_ids] indi_soc_vals_z; // blups on z-scale

  vector<lower=0>[n_periods] dyad_soc_sd;
  vector[n_pwise_dyads] dyad_soc_vals_z; // blups on z-scale

  vector[n_periods] beh_intercepts; // intercepts for behaviours

}
transformed parameters {
  vector[n_pwise_ids] indi_soc_vals = rep_vector(0.0, n_pwise_ids); // blups on actual scale
  vector[n_pwise_dyads] dyad_soc_vals = rep_vector(0.0, n_pwise_dyads); // blups on actual scale

  dyad_soc_vals = dyad_soc_vals_z .* dyad_soc_sd[index_period];
  indi_soc_vals = indi_soc_vals_z .* indi_soc_sd[index_period_individual];

  vector[n_pwise_dyads] scaled_indi_sums = sqrt(0.5) * (indi_soc_vals[index_id1] + indi_soc_vals[index_id2]);  // sums of dyadic values, scaled
}
model {
  vector[n_pwise_dyads] lp = rep_vector(0.0, n_pwise_dyads);
  lp += scaled_indi_sums;
  lp += dyad_soc_vals;
  lp += beh_intercepts[index_period];

  if (behav_types == 1) {
    interactions[, 1] ~ poisson_log(lp + log(obseff[, 1]));
  }
  if (behav_types == 2) {
    interactions[, 1] ~ binomial_logit(obseff_int[, 1], lp);
  }

  for (i in 1:n_periods) {
    beh_intercepts[i] ~ student_t(3, prior_matrix[i, 1], prior_matrix[i, 2]);
    indi_soc_sd[i] ~ exponential(2);
    dyad_soc_sd[i] ~ exponential(2);
  }

  indi_soc_vals_z ~ std_normal();
  dyad_soc_vals_z ~ std_normal();
}

generated quantities {
  array[n_pwise_dyads, 1] int interactions_pred; // matrix for pp functions...
  vector[n_pwise_dyads] lp = rep_vector(0.0, n_pwise_dyads);
  lp += scaled_indi_sums;
  lp += dyad_soc_vals;
  lp += beh_intercepts[index_period];

  if (behav_types == 1) {
    interactions_pred[, 1] = poisson_log_rng(lp + log(obseff[, 1]));
  }
  if (behav_types == 2) {
    interactions_pred[, 1] = binomial_rng(obseff_int[, 1], lp);
  }
}
