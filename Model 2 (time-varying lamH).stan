//--- Time-varying dengue catalytic model ---//
// assumes time-varying FOI prior to data
// assumes complete immunity after 2nd infection
// assumes equal transmissability of 4 serotypes

data {

  int nA; // N age groups
  int nT; // N time points
  int max_age; //max age of case data
  array[nT, nA] int cases; // reported case data
  matrix[nT, nA] pop; // population data
  array[2, nA] int ageLims; // lower & upper bounds of age groups
  row_vector[max_age] age; // age as sequence from 0 to 118

}

transformed data {

  array[nT] int sum_cases; // observed total cases per year

  for (t in 1 : nT)
    sum_cases[t] = sum(cases[t,  : ]);
}

parameters {

  array[5] real<lower=0, upper=0.25> lam_H; // historic average FOI
  array[nT] real<lower=0, upper=0.25> lam_t; // time varying FOI
  vector<lower=0, upper=1>[2] report;  //report[1] is reporting rate of 2nd infections, report[2] is the same for 1st infections
  real<lower=0, upper = min(1./report)> chi; // relative reporting rate of children aged 1-15 yrs old

}

transformed parameters {
  
  row_vector<lower=0, upper=1>[max_age] susc0; // proportion susceptible a time 0
  row_vector<lower=0, upper=1>[max_age] mono0; // proportion monotypic at time 0
  row_vector<lower=0, upper=1>[max_age] multi0; // proportion multitypic at time 0
  matrix<lower=0, upper=1>[nT, max_age] susc; // proportion susceptible
  matrix<lower=0, upper=1>[nT, max_age] mono; // proportion monotypic
  matrix<lower=0, upper=1>[nT, max_age] multi; // proportion multitypic
  matrix<lower=0, upper=1>[nT, max_age] inc1; // incidence of primary infections
  matrix<lower=0, upper=1>[nT, max_age] inc2; // incidence of secondary infections
  array[nT] vector<lower=0>[nA] Ecases; // expected reported cases per year and age group
  vector<lower=0>[nT] Ecases_sum; // expected reported cases per year
  array[nT] vector<lower=0>[nA] prob_Ecases; // expected probability of reported cases per year and age group

  //--- immune profiles at beginning of time series (time 0)
  
  for(i in 1 : 119){
    if(i <= 21){
      susc0[i] = exp(-4 * lam_H[1] * age[i]);
      mono0[i] = 4 * exp(-3 * lam_H[1] * age[i]) * (1 - exp(-lam_H[1] * age[i]));
    } else if(i > 21 && i <= 41){
      susc0[i] = exp(-4 * (20 * lam_H[1] + lam_H[2] * (age[i] - 20)));
      mono0[i] = 4 * exp(-3 * (20 * lam_H[1] + lam_H[2] * (age[i] - 20))) * 
      (1 - exp(-(20 * lam_H[1] + lam_H[2] * (age[i] - 20))));
    } else if (i > 41 && i <= 61){
      susc0[i] = exp(-4 * (20 * (lam_H[1] + lam_H[2]) + lam_H[3] * (age[i] - 40)));
      mono0[i] = 4 * exp(-3 * (20 * (lam_H[1] + lam_H[2]) + lam_H[3] * (age[i] - 40))) * 
      (1 - exp(-(20 * (lam_H[1] + lam_H[2]) + lam_H[3] * (age[i] - 40))));
    } else if(i > 61 && i <= 81){
      susc0[i] = exp(-4 * (20 * (lam_H[1] + lam_H[2] + lam_H[3]) + lam_H[4] * (age[i] - 60)));
      mono0[i] = 4 * exp(-3 * (20 * (lam_H[1] + lam_H[2] + lam_H[3]) + lam_H[4] * (age[i] - 60))) *
      (1 - exp(-(20 * (lam_H[1] + lam_H[2] + lam_H[3]) + lam_H[4] * (age[i] - 60))));
    } else{
      susc0[i] = exp(-4 * (20 * (lam_H[1] + lam_H[2] + lam_H[3] + lam_H[4]) + lam_H[5] * (age[i] - 80)));
      mono0[i] = 4 * exp(-3 * (20 * (lam_H[1] + lam_H[2] + lam_H[3] + lam_H[4]) + lam_H[5] * (age[i] - 80))) *
      (1 - exp(-(20 * (lam_H[1] + lam_H[2] + lam_H[3] + lam_H[4]) + lam_H[5] * (age[i] - 80))));
    }
    multi0[i] = 1 - (susc0[i] + mono0[i]);
  }

  //--- infants (assumes no infections in <1 year olds), we have observed cases also for below 1 year

  // for(t in 1:nT) susc[t,1] = 1;
  // for(t in 1:nT) mono[t,1] = 0;
  // for(t in 1:nT) multi[t,1] = 0;
  //--- immune profiles at time 1 
  
  susc[1, 1] = 1;

  mono[1, 1] = 0;

  multi[1, 1] = 0;

  susc[1, 2:max_age] = susc0[1:(max_age-1)] - 4 * lam_t[1] * susc0[1:(max_age-1)];

  mono[1, 2:max_age] = mono0[1:(max_age-1)] + 4 * lam_t[1] * susc0[1:(max_age-1)]
                     - 3 * lam_t[1] * mono0[1:(max_age-1)];

  multi[1, 2:max_age] = multi0[1:(max_age-1)] + 3 * lam_t[1] * mono0[1:(max_age-1)];

  inc1[1, ] = 4 * lam_t[1] * susc0;

  inc2[1, ] = 3 * lam_t[1] * mono0;

  //--- immune profiles at subsequent time steps 

  for (t in 2 : nT) {
    susc[t, 1] = 1;

    mono[t, 1] = 0;

    multi[t, 1] = 0;
    
    susc[t, 2:max_age] = susc[t - 1, 1:(max_age-1)]
                       - 4 * lam_t[t] * susc[t - 1, 1:(max_age-1)];

    mono[t, 2:max_age] = mono[t - 1, 1:(max_age-1)]
                       + 4 * lam_t[t] * susc[t - 1, 1:(max_age-1)]
                       - 3 * lam_t[t] * mono[t - 1, 1:(max_age-1)];

    multi[t, 2:max_age] = multi[t - 1, 1:(max_age-1)]
                        + 3 * lam_t[t] * mono[t - 1, 1:(max_age-1)];

    inc1[t, ] = 4 * lam_t[t] * susc[t - 1, ];

    inc2[t, ] = 3 * lam_t[t] * mono[t - 1, ];

  }

  //--- expected reported cases for each year and age group

  for (t in 1 : nT) 

    for (a in 1 : nA) {

      if (a <= 3){    // NEED TO CHANGE IF INCLUDING 0 YEAR OLDS IN LIKELIHOOD
      Ecases[t, a] = chi*(report[1]
                     * mean(inc2[t, ageLims[1, a] : ageLims[2, a]])
                        + report[2]
                          * mean(inc1[t, ageLims[1, a] : ageLims[2, a]]))
                     * pop[t, a];
      }
      else{
        Ecases[t, a] = ((report[1]
                     * mean(inc2[t, ageLims[1, a] : ageLims[2, a]]))
                        + (report[2]
                          * mean(inc1[t, ageLims[1, a] : ageLims[2, a]])))
                     * pop[t, a];
      }
    }

  // Ensure the log likelihood not 0

  for (t in 1 : nT) 

    for (a in 1 : nA) {

      if (Ecases[t, a] == 0) {
        Ecases[t, a] = 0.0001;

      } else {
        Ecases[t, a] = Ecases[t, a];
      }

    }

  //---  expected reported cases total per year

  for (t in 1 : nT) {
    Ecases_sum[t] = sum(Ecases[t,  : ]);
  }

  //--- Expected probability of reported cases per year and age group

  for (t in 1 : nT) {
    prob_Ecases[t,  : ] = Ecases[t,  : ] / sum(Ecases[t,  : ]);
  }

}

model {

  //--- priors
  lam_H ~ normal(0, 1);
  lam_t ~ normal(0, 1);
  report[1] ~ normal(0,1);
  report[2] ~ normal(0,1);
  chi ~ normal(1, 1); 

  //--- likelihood 
  // poisson likelihood for total per year

  for (t in 1 : nT) target += poisson_lpmf(sum_cases[t] | Ecases_sum[t]);

  // multinomial likelihood for age and year

  for (t in 1 : nT) target += multinomial_lpmf(cases[t,  : ] | prob_Ecases[t,  : ]);

}

generated quantities {
 
  vector[nT] log_lik;
  vector[nT] log_lik1;
  vector[nT] log_lik2;
 
  for (t in 1 : nT)     log_lik1[t] = poisson_lpmf(sum_cases[t] | Ecases_sum[t]);
 
  for (t in 1 : nT)     log_lik2[t] = multinomial_lpmf(cases[t,  : ] | prob_Ecases[t,  : ]);
 
  for (t in 1 : nT)     log_lik[t] = log_lik1[t] + log_lik2[t];
 
}
