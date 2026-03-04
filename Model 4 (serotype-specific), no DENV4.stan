//--- Time-varying dengue catalytic model ---//
// assumes constant endemic FOI prior to data
// assumes complete immunity after 2nd infection
// assumes equal transmissability of 4 serotypes

data {

  int nA; // N age groups
  int nT; // N time points
  int max_age; //max age of case data
  array[nT, nA, 3] int cases; // reported case data
  matrix[nT, nA] pop; // population data
  array[2, nA] int ageLims; // lower & upper bounds of age groups
  row_vector[max_age] age; // age as sequence 

}

transformed data {

  array[3, nT] int sum_cases; // observed total cases per year
  
  for (t in 1 : nT)
  for(s in 1 : 3){
    sum_cases[s, t] = sum(cases[t,  : , s]);
    }
}

parameters {

  array[5] simplex[3] lam_H; // historic average serotype FOI proportions
  array[nT] simplex[3] lam_t; // time varying yeraly serotype FOI proportions 
  vector<lower=0, upper=1>[2] report; //primary and secondary reporting rates now combined into one vector
  real<lower=0, upper = min(1./report)> chi; // relative reporting rate of children aged 1-15 yrs old
  array[5] real<lower=0, upper=1> sigma_H; //total historic FOI (multiplier for lam_H as using dirichlet distribution)
  array[nT] real<lower=0, upper=1> sigma_t; ; //yearly total FOI (multiplier for lam_t as using dirichlet distribution)

}

transformed parameters {

  row_vector<lower=0, upper=1>[max_age] susc0; // proportion susceptible a time 0
  array[3] row_vector<lower=0, upper=1>[max_age] mono0; // proportion monotypic at time 0
  row_vector<lower=0, upper=1>[max_age] multi0; // proportion multitypic at time 0
  matrix<lower=0, upper=1>[nT, max_age] susc; // proportion susceptible
  array[3] matrix<lower=0, upper=1>[nT, max_age] mono; // proportion monotypic
  matrix<lower=0, upper=1>[nT, max_age] multi; // proportion multitypic
  array[3] matrix<lower=0, upper=1>[nT, max_age] inc1; // incidence of primary infections
  array[3] matrix<lower=0, upper=1>[nT, max_age] inc2; // incidence of secondary infections
  array[nT, 3] vector<lower=0>[nA] Ecases; // expected reported cases per year and age group
  array[3] vector<lower=0>[nT] Ecases_sum; // expected reported cases per year
  array[nT, 3] vector<lower=0>[nA] prob_Ecases; // expected probability of reported cases per year and age group

  //--- immune profiles at beginning of time series (time 0)
  for(i in 1 : max_age){
  if(i <= 21){
    susc0[i] = exp(-sum(sigma_H[1]*lam_H[1,]) * age[i]);
    for(s in 1:3){
      mono0[s, i] = exp(-(sum(sigma_H[1]*lam_H[1,]) - sigma_H[1]*lam_H[1,s]) * age[i]) .* 
        (1 - exp(-sigma_H[1]*lam_H[1,s] * age[i]));
    }
  } else if(i > 21 && i <= 41){
    susc0[i] = exp(-(sum(sigma_H[1]*lam_H[1,]) * 20 + sum(sigma_H[2]*lam_H[2,]) * (age[i] - 20)));
    for(s in 1:3){
      mono0[s, i] = exp(-((sum(sigma_H[1]*lam_H[1,]) - sigma_H[1]*lam_H[1,s]) * 20 + 
                            (sum(sigma_H[2]*lam_H[2,]) - sigma_H[2]*lam_H[2,s]) * (age[i] - 20))) .*
        (1 - exp(-(sigma_H[1]*lam_H[1,s] * 20 + sigma_H[2]*lam_H[2,s] * (age[i] - 20))));
    }
  } else if (i > 41 && i <= 61){
    susc0[i] = exp(-(sum(sigma_H[1]*lam_H[1,]) * 20 + sum(sigma_H[2]*lam_H[2,]) * 20 + 
                    sum(sigma_H[3]*lam_H[3,]) * (age[i] - 40)));
    for(s in 1:3){
      mono0[s, i] = exp(-((sum(sigma_H[1]*lam_H[1,]) - sigma_H[1]*lam_H[1,s]) * 20 + 
                            (sum(sigma_H[2]*lam_H[2,]) - sigma_H[2]*lam_H[2,s]) * 20 + 
                            (sum(sigma_H[3]*lam_H[3,]) - sigma_H[3]*lam_H[3,s]) * (age[i] - 40))) .*
        (1 - exp(-(sigma_H[1]*lam_H[1,s] * 20 + sigma_H[2]*lam_H[2,s] * 20 + sigma_H[3]*lam_H[3,s] * (age[i] - 40))));
    }
  } else if(i > 61 && i <= 81){
    susc0[i] = exp(-(sum(sigma_H[1]*lam_H[1,]) * 20 + sum(sigma_H[2]*lam_H[2,]) * 20 + 
                    sum(sigma_H[3]*lam_H[3,]) * 20 + sum(sigma_H[4]*lam_H[4,]) * (age[i] - 60)));
    for(s in 1:3){
      mono0[s, i] = exp(-((sum(sigma_H[1]*lam_H[1,]) - sigma_H[1]*lam_H[1,s]) * 20 + 
                            (sum(sigma_H[2]*lam_H[2,]) - sigma_H[2]*lam_H[2,s]) * 20 + 
                            (sum(sigma_H[3]*lam_H[3,]) - sigma_H[3]*lam_H[3,s]) * 20 + 
                            (sum(sigma_H[4]*lam_H[4,]) - sigma_H[4]*lam_H[4,s]) * (age[i] - 60))) .*
        (1 - exp(-(sigma_H[1]*lam_H[1,s] * 20 + sigma_H[2]*lam_H[2,s] * 20 + sigma_H[3]*lam_H[3,s] * 20 + 
                     sigma_H[4]*lam_H[4,s] * (age[i] - 60))));
    }
  } else{
    susc0[i] = exp(-(sum(sigma_H[1]*lam_H[1,]) * 20 + sum(sigma_H[2]*lam_H[2,]) * 20 + 
                    sum(sigma_H[3]*lam_H[3,]) * 20 + sum(sigma_H[4]*lam_H[4,]) * 20 +
                    sum(sigma_H[5]*lam_H[5,]) * (age[i] - 80)));
    for(s in 1:3){
      mono0[s, i] = exp(-((sum(sigma_H[1]*lam_H[1,]) - sigma_H[1]*lam_H[1,s]) * 20 + 
                            (sum(sigma_H[2]*lam_H[2,]) - sigma_H[2]*lam_H[2,s]) * 20 + 
                            (sum(sigma_H[3]*lam_H[3,]) - sigma_H[3]*lam_H[3,s]) * 20 + 
                            (sum(sigma_H[4]*lam_H[4,]) - sigma_H[4]*lam_H[4,s]) * 20 + 
                            (sum(sigma_H[5]*lam_H[5,]) - sigma_H[5]*lam_H[5,s]) * (age[i] - 80))) .*
        (1 - exp(-(sigma_H[1]*lam_H[1,s] * 20 + sigma_H[2]*lam_H[2,s] * 20 + sigma_H[3]*lam_H[3,s] * 20 + 
                     sigma_H[4]*lam_H[4,s] * 20 + sigma_H[5]*lam_H[5,s] * (age[i] - 80))));
    }
  }
  multi0[i] = 1 - susc0[i] - sum(mono0[, i]);
}

  //--- infants (assumes no infections in <1 year olds), we have observed cases also for below 1 year

  // for(t in 1:nT) susc[t,1] = 1;
  // for(t in 1:nT) mono[t,1] = 0;
  // for(t in 1:nT) multi[t,1] = 0;
  //--- immune profiles at time 1 
  
  susc[1, 1] = 1;
  
  multi[1, 1] = 0;
  
  susc[1, 2:max_age] = susc0[1:(max_age-1)] - sum(sigma_t[1]*lam_t[1]) * susc0[1:(max_age-1)];
  
  
  for(s in 1:3){
  mono[s, 1, 1] = 0;

  mono[s, 1, 2:max_age] = mono0[s, 1:(max_age-1)] + sigma_t[1]*lam_t[1][s] * susc0[1:(max_age-1)]
                     - (sum(sigma_t[1]*lam_t[1]) - sigma_t[1]*lam_t[1][s]) * mono0[s, 1:(max_age-1)];

  inc1[s, 1, ] = sigma_t[1]*lam_t[1][s] * susc0;
  
  for(a in 1:max_age){
  inc2[s, 1, a] = sigma_t[1]*lam_t[1][s] * (sum(mono0[, a]) - mono0[s,a]);
  
  }}

  for(a in 2:max_age){
    multi[1, a] = multi0[(a-1)] + sum(inc2[, 1, (a-1)]);
  }
  
  
  
  

  //--- immune profiles at subsequent time steps 

if(nT >= 2){
  for (t in 2 : nT) {
    susc[t, 1] = 1;

    multi[t, 1] = 0;
    
    susc[t, 2:max_age] = susc[t - 1, 1:(max_age-1)]
                       - sum(sigma_t[t]*lam_t[t]) * susc[t - 1, 1:(max_age-1)];
                       
    for(s in 1:3){    
    mono[s, t, 1] = 0;

    mono[s, t, 2:max_age] = mono[s, t - 1, 1:(max_age-1)]
                       + sigma_t[t]*lam_t[t][s] * susc[t - 1, 1:(max_age-1)]
                       - (sum(sigma_t[t]*lam_t[t]) - sigma_t[t]*lam_t[t][s]) * mono[s, t - 1, 1:(max_age-1)];

    inc1[s, t, ] = sigma_t[t]*lam_t[t][s] * susc[t - 1, ];
    
    for(a in 1:max_age){
    inc2[s, t, a] = sigma_t[t]*lam_t[t][s] * (sum(mono[, t - 1, a]) - mono[s, t - 1, a]);
    }}
                       
    for(a in 2:max_age){
    multi[t, a] = multi[t - 1, (a - 1)] + sum(inc2[, t, (a - 1)]);
    }                   

  }
}

  //--- expected reported cases for each year and age group

  for (t in 1 : nT) 

    for (a in 1 : nA) {
      
      for(s in 1 : 3){
        
        if (a <= 3){ 

      Ecases[t, s, a] = chi*(report[1]
                     * mean(inc2[s, t, ageLims[1, a] : ageLims[2, a]])
                        + report[2]
                          * mean(inc1[s, t, ageLims[1, a] : ageLims[2, a]]))
                     * pop[t, a];
    }
    else{
      
       Ecases[t, s, a] = (report[1]
                     * mean(inc2[s, t, ageLims[1, a] : ageLims[2, a]])
                        + report[2]
                          * mean(inc1[s, t, ageLims[1, a] : ageLims[2, a]]))
                     * pop[t, a];
    }
      }
    }

  // Apply your conditional logic to ensure that the log likelihood never works with 0 but close (0.0001)

  for (t in 1 : nT) 

    for (a in 1 : nA) {
      
      for(s in 1:3){

      if (Ecases[t, s, a] == 0) {
        Ecases[t, s, a] = 0.0001;

      } else {
        Ecases[t, s, a] = Ecases[t, s, a];
      }

    }
    }

  //---  expected reported cases total per year

  for (t in 1 : nT) {
    for(s in 1:3){
    Ecases_sum[s, t] = sum(Ecases[t, s, : ]);
  }
  }

  //--- Expected probability of reported cases per year and age group

  for (t in 1 : nT) {
    for(s in 1:3){
    prob_Ecases[t, s, : ] = Ecases[t, s, : ] / sum(Ecases[t, s, : ]);
  }
  }

}

model {

  //--- priors
  row_vector[3] alpha = [1, 1, 1]; //so dirichlet prior symmetric
  for (t in 1:nT) {
    lam_t[t] ~ dirichlet(alpha);
  }
  for (i in 1:5) {
  lam_H[i,] ~ dirichlet(alpha);
  }
  report ~ normal(0, 1); 
  sigma_H ~ uniform(0,1);
  sigma_t ~ uniform(0,1);
  chi ~ normal(1, 1); 
  //--- likelihood 
  // poisson likelihood for total per year

  for(s in 1 : 3) for (t in 1 : nT) target += poisson_lpmf(sum_cases[s, t] | Ecases_sum[s, t]);

  // multinomial likelihood for age and year

  for(s in 1 : 3) for (t in 1 : nT) target += multinomial_lpmf(cases[t,  : , s] | prob_Ecases[t, s, : ]);

}

generated quantities {
 
  array[3] vector[nT] log_lik;
  array[3] vector[nT] log_lik1;
  array[3] vector[nT] log_lik2;
 
  for(s in 1 : 3) for (t in 1 : nT)     log_lik1[s, t] = poisson_lpmf(sum_cases[s, t] | Ecases_sum[s, t]);
 
  for(s in 1 : 3) for (t in 1 : nT)     log_lik2[s, t] = multinomial_lpmf(cases[t,  : , s] | prob_Ecases[t, s, : ]);
 
  for(s in 1 : 3)  for (t in 1 : nT)     log_lik[s, t] = log_lik1[s, t] + log_lik2[s, t];
 
}
