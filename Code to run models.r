library(tidyverse)
library(rstan)

## Running the Model 1, from Imai et al. 2016##

#Case data

mex_cases_age_state_no0 <- readRDS("Mexico dengue data/all_case_data.rds")

mex_cases_mat_no0 <- as.matrix(mex_cases_age_state_no0[, 3:18])
mex_states_order <- unique(mex_cases_age_state_no0[, 2])

#Pop data

mex_pop_age_state_no0 <- readRDS("Mexico dengue data/pop_data.rds") %>%
  mutate(NOM_ENT = ifelse(NOM_ENT == "ciudad de m'exico", "distrito federal", NOM_ENT)) %>%
  arrange(NOM_ENT, Year)

mex_pop_mat_no0 <- as.matrix(mex_pop_age_state_no0[, 3:18])

#Run model

mod_imai <- rstan::stan_model(file = "Model from Imai et al. 2016 (Model 1).stan")

nT <- 8 #Number of years running model for
nA <- ncol(mex_cases_mat_no0) #Number of age groups

age_min_imai <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75) #Minimum age in age groups
age_max_imai <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 119) #Maximum age in age groups

fit_imai <- vector(mode = "list", length = 32)

for (i in seq(1, 32)[-c(2, 6, 29)]) {
  data <- list(
    nA = nA,
    nT = nT,
    cases = mex_cases_mat_no0[((i - 1) * 8 + 1):(i * 8), ],
    pop = mex_pop_mat_no0[((i - 1) * 8 + 1):(i * 8), ],
    ageLims = rbind(age_min_imai, age_max_imai)
  )
  
  
  fit_imai[[i]] <- rstan::sampling(
    mod_imai,
    data = data,
    chains = 3,
    #number of chains
    iter = 6000,
    #number of iterations for each chain
    warmup = 1000,
    #number of these iterations which are warmup
    cores = 3,
    #run chains in parallel
    refresh = 100
  )
}

chains_imai <- vector(mode = "list", length = 32)
fit_summary_imai <- vector(mode = "list", length = 32)
diagnostics_imai <- vector(mode = "list", length = 32)

for (i in seq(1, 32)[-c(2, 6, 29)]) {
  chains_imai[[i]] <- rstan::extract(fit_imai[[i]])
  fit_summary_imai[[i]] <- rstan::summary(fit_imai[[i]])
  diagnostics_imai[[i]] <- cbind(fit_summary_imai[[i]]$summary[, "Rhat"], fit_summary_imai[[i]]$summary[, "n_eff"])
}

#Extract parameter posteriors

rho_imai <- vector(mode = "list", length = 32)
gamma1_imai <- vector(mode = "list", length = 32)
lam_imai <- vector(mode = "list", length = 32)
pars_imai <- vector(mode = "list", length = 32)

for (i in seq(1, 31)[-c(2, 6, 29)]) {
  rho_imai[[i]] <-
    quantile(chains_imai[[i]]$rho, c(0.5, 0.025, 0.975))
  gamma1_imai[[i]] <-
    quantile(chains_imai[[i]]$gamma1, c(0.5, 0.025, 0.975))
  lam_imai[[i]] <-
    quantile(chains_imai[[i]]$lam, c(0.5, 0.025, 0.975))
  pars_imai[[i]] <- rbind(rho_imai[[i]], gamma1_imai[[i]], lam_imai[[i]])
  pars_imai[[i]] <- as.data.frame(pars_imai[[i]])
  pars_imai[[i]]$pars <- c("rho", "gamma", "lam")
  pars_imai[[i]]$ENTIDAD_RES <- rep(as.character(mex_states_order[i, ]), 3)
  colnames(pars_imai[[i]]) <- c("med", "ciL", "ciU", "pars", "ENTIDAD_RES")
}

pars_df_imai <- bind_rows(pars_imai)

## Running Model 5, based on O'Driscoll et al. 2019##

#Case data

mex_cases_age_state_no0 <- readRDS("Mexico dengue data/all_case_data.rds")

mex_cases_mat_no0 <- as.matrix(mex_cases_age_state_no0[, 3:18])
mex_states_order <- unique(mex_cases_age_state_no0[, 2])

#Pop data

mex_pop_age_state_no0 <- mex_pop_age_state_no0 <- readRDS("Mexico dengue data/pop_data.rds") %>%
  mutate(NOM_ENT = ifelse(NOM_ENT == "ciudad de m'exico", "distrito federal", NOM_ENT)) %>%
  arrange(NOM_ENT, Year)

mex_pop_mat_no0 <- as.matrix(mex_pop_age_state_no0[, 3:18])

#Run model function

run_model <- function(state, mex_cases_mat, mex_pop_mat) {
  mod <- rstan::stan_model(file = "Model 2 (time-varying lamH).stan")
  
  nT <- 8  #Number of years running model for
  nA <- ncol(mex_cases_mat) #Number of age groups
  
  age_min <- c(2, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76) #Minimum age in age groups
  age_max <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 119) #Maximum age in age groups
  age <- seq(0, 118) #Sequence of total number of age years
  
  data <- list(
    nA = nA,
    nT = nT,
    cases = mex_cases_mat[((state - 1) * 8 + 1):(state * 8), ],
    pop = mex_pop_mat[((state - 1) * 8 + 1):(state * 8), ],
    age = age,
    ageLims = rbind(age_min, age_max),
    max_age = 119
  )
  
  fit <- rstan::sampling(
    mod,
    data = data,
    chains = 3,
    #number of chains
    iter = 6000,
    #number of iterations for each chain
    warmup = 1000,
    #number of these iterations which are warmup
    cores = 3,
    #run chains in parallel
    refresh = 100
  )
  
  saveRDS(fit, file  = paste0("fit_output_m2_", state, ".rds"))
  
  chains <- rstan::extract(fit)
  saveRDS(chains, file  = paste0("chains_m2_", state, ".rds"))
  
  fit_summary <- rstan::summary(fit)
  saveRDS(fit_summary, file  = paste0("fit_summary_m2_", state, ".rds"))
  
  rhat_values <- fit_summary$summary[, "Rhat"]
  ess_bulk_values <- fit_summary$summary[, "n_eff"]
  
  diagnostics <- cbind(rhat_values, ess_bulk_values)
  saveRDS(diagnostics, file  = paste0("diagnostics_m2_", state, ".rds"))
  
}

#eg:

run_model(17, mex_cases_mat_no0, mex_pop_mat_no0) #State 17 (alphabetically) is Morelos, check using mex_states_order

chains_m2 <- vector(mode = "list", length = 32)

#for(i in seq(1, 31)[-c(2, 6, 29)]){ #this is all states
for (i in 17:17) {
  chains_m2[[i]] <- readRDS(paste0("chains_m2_", i, ".rds"))
}

#Extract parameter posteriors

rho_m2 <- vector(mode = "list", length = 32)
gamma_m2 <- vector(mode = "list", length = 32)
rho_young_m2 <- vector(mode = "list", length = 32)
gamma_young_m2 <- vector(mode = "list", length = 32)
chi_m2 <- vector(mode = "list", length = 32)
lam_H_m2 <- vector(mode = "list", length = 32)
pars_m2 <- vector(mode = "list", length = 32)

#for (i in seq(1, 31)[-c(2, 6, 29)]) { #this is all states
for (i in 17:17) {
  rho_m2[[i]] <-
    quantile(chains_m2[[i]]$report[, 1], c(0.5, 0.025, 0.975))
  gamma_m2[[i]] <-
    quantile(chains_m2[[i]]$report[, 2], c(0.5, 0.025, 0.975))
  chi_m2[[i]] <-
    quantile(chains_m2[[i]]$chi, c(0.5, 0.025, 0.975))
  lam_H_1_m2[[i]] <-
    quantile(chains_m2[[i]]$lam_H[, 1], c(0.5, 0.025, 0.975))
  lam_H_2_m2[[i]] <-
    quantile(chains_m2[[i]]$lam_H[, 2], c(0.5, 0.025, 0.975))
  lam_H_3_m2[[i]] <-
    quantile(chains_m2[[i]]$lam_H[, 3], c(0.5, 0.025, 0.975))
  lam_H_4_m2[[i]] <-
    quantile(chains_m2[[i]]$lam_H[, 4], c(0.5, 0.025, 0.975))
  lam_H_5_m2[[i]] <-
    quantile(chains_m2[[i]]$lam_H[, 5], c(0.5, 0.025, 0.975))
  rho_young_m2[[i]] <-
    quantile(chains_m2[[i]]$report[, 1] * chains_m2[[i]]$chi, c(0.5, 0.025, 0.975))
  gamma_young_m2[[i]] <-
    quantile(chains_m2[[i]]$report[, 2] * chains_m2[[i]]$chi, c(0.5, 0.025, 0.975))
  pars_m2[[i]] <- rbind(
    rho_m2[[i]],
    gamma_m2[[i]],
    chi_m2[[i]],
    lam_H_1_m2[[i]],
    lam_H_2_m2[[i]],
    lam_H_3_m2[[i]],
    lam_H_4_m2[[i]],
    lam_H_5_m2[[i]],
    rho_young_m2[[i]],
    gamma_young_m2[[i]]
  )
  pars_m2[[i]] <- as.data.frame(pars_m2[[i]])
  pars_m2[[i]]$pars <- c(
    "rho",
    "gamma",
    "chi",
    "lam_H_1",
    "lam_H_2",
    "lam_H_3",
    "lam_H_4",
    "lam_H_5",
    "rho_young",
    "gamma_young"
  )
  pars_m2[[i]]$ENTIDAD_RES <- rep(as.character(mex_states_order[i, ]), 10)
  colnames(pars_m2[[i]]) <- c("med", "ciL", "ciU", "pars", "ENTIDAD_RES")
}

pars_df_m2 <- bind_rows(pars_m2)

lam_m2 <- vector(mode = "list", length = 32)

#for (i in seq(1, 31)[-c(2, 6, 9, 29)]) { #this is all states
for (i in 17:17) {
  lam_m2[[i]] <- data.frame(
    t = seq(1, data$nT),
    med = NA,
    ciL = NA,
    ciU = NA,
    ENTIDAD_RES = mex_states_order[i, ]
  )
  
  for (t in 1:data$nT)
    lam_m2[[i]][t, 2:4] <-
      quantile(chains_m2[[i]]$lam_t[, t], c(0.5, 0.025, 0.975))
  
  lam_m2[[i]]$type <-
    c("lam_T1",
      "lam_T2",
      "lam_T3",
      "lam_T4",
      "lam_T5",
      "lam_T6",
      "lam_T7",
      "lam_T8")
}

lam_df_m2 <- bind_rows(lam_m2)


## Running Model 6, the serotype-specific model##

#Attribute cases to each serotype

mex_cases_age_state_no0 <- readRDS("Mexico dengue data/all_case_data.rds")
mex_states_order <- unique(mex_cases_age_state_no0[, 2])

prop_sero_1 <- readRDS("Mexico dengue data/proportion_of_serotype.rds") %>%
  filter(serotype == "DENV1") %>%
  select(!serotype)

prop_sero_2 <- readRDS("Mexico dengue data/proportion_of_serotype.rds") %>%
  filter(serotype == "DENV2") %>%
  select(!serotype)

prop_sero_3 <- readRDS("Mexico dengue data/proportion_of_serotype.rds") %>%
  filter(serotype == "DENV3") %>%
  select(!serotype)

prop_sero_4 <- readRDS("Mexico dengue data/proportion_of_serotype.rds") %>%
  filter(serotype == "DENV4") %>%
  select(!serotype)

hosp_cases_no0 <- readRDS("Mexico dengue data/hospitalised_case_data.rds")

non_hosp_cases_no0 <- readRDS("Mexico dengue data/non_hospitalised_case_data.rds")

hosp_cases_mat_no0 <- as.matrix(hosp_cases_no0[, 3:18])
non_hosp_cases_mat_no0 <- as.matrix(non_hosp_cases_no0[, 3:18])

hosp_cases_mat_no0[is.na(hosp_cases_mat_no0)] <- 0
non_hosp_cases_mat_no0[is.na(non_hosp_cases_mat_no0)] <- 0

cases_mat_no0_1 = diag(prop_sero_1$prop_hosp) %*% hosp_cases_mat_no0 + diag(prop_sero_1$prop_non_hosp) %*% non_hosp_cases_mat_no0
cases_mat_no0_2 = diag(prop_sero_2$prop_hosp) %*% hosp_cases_mat_no0 + diag(prop_sero_2$prop_non_hosp) %*% non_hosp_cases_mat_no0
cases_mat_no0_3 = diag(prop_sero_3$prop_hosp) %*% hosp_cases_mat_no0 + diag(prop_sero_3$prop_non_hosp) %*% non_hosp_cases_mat_no0
cases_mat_no0_4 = diag(prop_sero_4$prop_hosp) %*% hosp_cases_mat_no0 + diag(prop_sero_4$prop_non_hosp) %*% non_hosp_cases_mat_no0
cases_array_seros <- vector(mode = "list", length = 4)
cases_array_seros[[1]] <- cases_mat_no0_1
cases_array_seros[[2]] <- cases_mat_no0_2
cases_array_seros[[3]] <- cases_mat_no0_3
cases_array_seros[[4]] <- cases_mat_no0_4
cases <- array(unlist(cases_array_seros), dim = c(nT * 32, nA, 4))

#Pop data

mex_pop_age_state_no0 <- mex_pop_age_state_no0 <- readRDS("Mexico dengue data/Mexico dengue data/pop_data.rds") %>%
  mutate(NOM_ENT = ifelse(NOM_ENT == "ciudad de m'exico", "distrito federal", NOM_ENT)) %>%
  arrange(NOM_ENT, Year)

mex_pop_mat_no0 <- as.matrix(mex_pop_age_state_no0[, 3:18])

#Run model function

run_model_m4 <- function(state, cases, pop) {
  age_min <- c(2, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76)
  age_max <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 119)
  age <- seq(0, 118)
  
  
  if (state %in% c(5, 12, 20, 21, 23, 27, 30, 31)) {
    mod <- rstan::stan_model(file = paste0("Model 4 (serotype-specific), all serotypes.stan"))
  } else if (state %in% c(10, 11, 18, 24)) {
    mod <- rstan::stan_model(file = paste0("Model 4 (serotype-specific), no DENV3-4.stan"))
  } else{
    mod <- rstan::stan_model(file = paste0("Model 4 (serotype-specific), no DENV4.stan"))
  }
  data <- list(
    nA = 16,
    #numb of age groups
    cases = if (state %in% c(1, 9, 13, 22)) {
      array(round(cases[(state * 8), , 1:3]), dim = c(1, 16, 3))
    } else if (state == 3) {
      round(cases[((state * 8) - 1):(state * 8), , 1:3])
    } else if (state %in% c(23, 31)) {
      round(cases[((state * 8) - 1):(state * 8), , ])
    } else if (state == 10) {
      round(cases[((state - 1) * 8 + 1):((state - 1) * 8 + 2), , 1:2])
    } else if (state == 26) {
      round(cases[((state * 8) - 2):(state * 8), , 1:3])
    } else if (state == 18) {
      round(cases[((state * 8) - 5):((state * 8) - 2), , 1:2])
    } else if (state == 24) {
      round(cases[((state - 1) * 8 + 2):((state * 8) - 2), , 1:2])
    } else if (state == 4) {
      round(cases[((state * 8) - 4):(state * 8), , 1:3])
    } else if (state == 28) {
      round(cases[((state - 1) * 8 + 1):((state * 8) - 3), , 1:3])
    } else if (state %in% c(8, 25)) {
      round(cases[((state - 1) * 8 + 3):(state * 8), , 1:3])
    } else if (state == 11) {
      round(cases[((state - 1) * 8 + 1):((state * 8) - 2), , 1:2])
    } else if (state == 7) {
      round(cases[((state - 1) * 8 + 2):(state * 8), , 1:3])
    } else if (state %in% c(5, 12, 20, 21, 27, 30)) {
      round(cases[((state - 1) * 8 + 1):(state * 8), , ])
    } else{
      round(cases[((state - 1) * 8 + 1):(state * 8), , 1:3])
    },
    #case matrix
    pop = if (state %in% c(1, 9, 13, 22)) {
      array(pop[(state * 8), ], dim = c(1, 16))
    } else if (state %in% c(3, 23, 31)) {
      pop[((state * 8) - 1):(state * 8), ]
    } else if (state == 10) {
      pop[((state - 1) * 8 + 1):((state - 1) * 8 + 2), ]
    } else if (state == 26) {
      pop[((state * 8) - 2):(state * 8), ]
    } else if (state == 18) {
      pop[((state * 8) - 5):((state * 8) - 2), ]
    } else if (state == 24) {
      pop[((state - 1) * 8 + 2):((state * 8) - 2), ]
    } else if (state == 4) {
      pop[((state * 8) - 4):(state * 8), ]
    } else if (state == 28) {
      pop[((state - 1) * 8 + 1):((state * 8) - 3), ]
    } else if (state %in% c(8, 25)) {
      pop[((state - 1) * 8 + 3):(state * 8), ]
    } else if (state == 11) {
      pop[((state - 1) * 8 + 1):((state * 8) - 2), ]
    } else if (state == 7) {
      pop[((state - 1) * 8 + 2):(state * 8), ]
    } else{
      pop[((state - 1) * 8 + 1):(state * 8), ]
    },
    #population matrix
    age = age,
    ageLims = rbind(age_min, age_max),
    max_age = 119
  )
  data$nT = if (state %in% c(1, 9, 13, 22)) {
    1
  } else{
    dim(data$pop)[1]
  }
  
  if (model == "model_m6_c2_vh") {
    data$first_year = if (state %in% c(1, 9, 13, 22)) {
      8
    } else if (state %in% c(3, 23, 31)) {
      7
    } else if (state == 26) {
      6
    } else if (state == 4) {
      4
    } else if (state %in% c(8, 18, 25)) {
      3
    } else if (state %in% c(7, 24)) {
      2
    } else{
      1
    }
  }
  
  fit <- rstan::sampling(
    mod,
    data = data,
    chains = 3,
    #number of chains
    iter = 6000,
    #number of iterations for each chain
    warmup = 1000,
    #number of these iterations which are warmup
    cores = 3,
    #run chains in parallel
    refresh = 100,
    #control = list(adapt_delta = 0.95),
    seed = 2000 #seed 1
    #seed = 6000 #seed 2
    #seed = 7000 #seed 3
  )
  
  
  saveRDS(fit, file  = paste0("fit_output_m4_", state, ".rds"))
  
  chains <- rstan::extract(fit)
  saveRDS(chains, file  = paste0("chains_m4_", state, ".rds"))
  
  fit_summary <- rstan::summary(fit)
  saveRDS(fit_summary, file  = paste0("fit_summary_m4_", state, ".rds"))
  
  rhat_values <- fit_summary$summary[, "Rhat"]
  ess_bulk_values <- fit_summary$summary[, "n_eff"]
  
  diagnostics <- cbind(rhat_values, ess_bulk_values)
  saveRDS(diagnostics, file  = paste0("diagnostics_m4_", state, ".rds"))
  
}

#eg:

run_model_m4(5, cases, mex_pop_mat_no0)  #State 5 (alphabetically) is Chiapas, check using mex_states_order

chains_m4 <- vector(mode = "list", length = 32)

#for(i in seq(1, 31)[-c(2, 6, 10, 11, 15, 29, 31)]){ #this is all states
for (i in 5:5) {
  chains_m4[[i]] <- readRDS(paste0("chains_m4_", i, ".rds"))
}

#Extract parameter posteriors

rho_m4 <- vector(mode = "list", length = 32)
gamma_m4 <- vector(mode = "list", length = 32)
chi_m4 <- vector(mode = "list", length = 32)
lam_H_1_1_m4 <- vector(mode = "list", length = 32)
lam_H_2_1_m4 <- vector(mode = "list", length = 32)
lam_H_3_1_m4 <- vector(mode = "list", length = 32)
lam_H_4_1_m4 <- vector(mode = "list", length = 32)
lam_H_1_2_m4 <- vector(mode = "list", length = 32)
lam_H_2_2_m4 <- vector(mode = "list", length = 32)
lam_H_3_2_m4 <- vector(mode = "list", length = 32)
lam_H_4_2_m4 <- vector(mode = "list", length = 32)
lam_H_1_3_m4 <- vector(mode = "list", length = 32)
lam_H_2_3_m4 <- vector(mode = "list", length = 32)
lam_H_3_3_m4 <- vector(mode = "list", length = 32)
lam_H_4_3_m4 <- vector(mode = "list", length = 32)
lam_H_1_4_m4 <- vector(mode = "list", length = 32)
lam_H_2_4_m4 <- vector(mode = "list", length = 32)
lam_H_3_4_m4 <- vector(mode = "list", length = 32)
lam_H_4_4_m4 <- vector(mode = "list", length = 32)
lam_H_1_5_m4 <- vector(mode = "list", length = 32)
lam_H_2_5_m4 <- vector(mode = "list", length = 32)
lam_H_3_5_m4 <- vector(mode = "list", length = 32)
lam_H_4_5_m4 <- vector(mode = "list", length = 32)
pars_m4 <- vector(mode = "list", length = 32)

for (i in 5:5) {
#for(i in seq(1, 31)[-c(2, 6, 9, 10, 29)]){ this is all states
  
  rho_m4[[i]] <-
    quantile(chains_m4[[i]]$report[, 1], c(0.5, 0.025, 0.975))
  gamma_m4[[i]] <-
    quantile(chains_m4[[i]]$report[, 2], c(0.5, 0.025, 0.975))
  chi_m4[[i]] <-
    quantile(chains_m4[[i]]$chi, c(0.5, 0.025, 0.975))
  lam_H_1_1_m4[[i]] <-
    quantile(chains_m4[[i]]$lam_H[, 1, 1] * chains_m4[[i]]$sigma_H[, 1], c(0.5, 0.025, 0.975))
  lam_H_2_1_m4[[i]] <-
    quantile(chains_m4[[i]]$lam_H[, 1, 2] * chains_m4[[i]]$sigma_H[, 1], c(0.5, 0.025, 0.975))
  lam_H_1_2_m4[[i]] <-
    quantile(chains_m4[[i]]$lam_H[, 2, 1] * chains_m4[[i]]$sigma_H[, 2], c(0.5, 0.025, 0.975))
  lam_H_2_2_m4[[i]] <-
    quantile(chains_m4[[i]]$lam_H[, 2, 2] * chains_m4[[i]]$sigma_H[, 2], c(0.5, 0.025, 0.975))
  lam_H_1_3_m4[[i]] <-
    quantile(chains_m4[[i]]$lam_H[, 3, 1] * chains_m4[[i]]$sigma_H[, 3], c(0.5, 0.025, 0.975))
  lam_H_2_3_m4[[i]] <-
    quantile(chains_m4[[i]]$lam_H[, 3, 2] * chains_m4[[i]]$sigma_H[, 3], c(0.5, 0.025, 0.975))
  lam_H_1_4_m4[[i]] <-
    quantile(chains_m4[[i]]$lam_H[, 4, 1] * chains_m4[[i]]$sigma_H[, 4], c(0.5, 0.025, 0.975))
  lam_H_2_4_m4[[i]] <-
    quantile(chains_m4[[i]]$lam_H[, 4, 2] * chains_m4[[i]]$sigma_H[, 4], c(0.5, 0.025, 0.975))
  lam_H_1_5_m4[[i]] <-
    quantile(chains_m4[[i]]$lam_H[, 5, 1] * chains_m4[[i]]$sigma_H[, 5], c(0.5, 0.025, 0.975))
  lam_H_2_5_m4[[i]] <-
    quantile(chains_m4[[i]]$lam_H[, 5, 2] * chains_m4[[i]]$sigma_H[, 5], c(0.5, 0.025, 0.975))
  
  
  if (i %in% c(11, 18, 24)) {
    lam_H_3_1_m4[[i]] <- rep(NA, 3)
    lam_H_4_1_m4[[i]] <- rep(NA, 3)
    lam_H_3_2_m4[[i]] <- rep(NA, 3)
    lam_H_4_2_m4[[i]] <- rep(NA, 3)
    lam_H_3_3_m4[[i]] <- rep(NA, 3)
    lam_H_4_3_m4[[i]] <- rep(NA, 3)
    lam_H_3_4_m4[[i]] <- rep(NA, 3)
    lam_H_4_4_m4[[i]] <- rep(NA, 3)
    lam_H_3_5_m4[[i]] <- rep(NA, 3)
    lam_H_4_5_m4[[i]] <- rep(NA, 3)
    
  }
  else{
    lam_H_3_1_m4[[i]] <-
      quantile(chains_m4[[i]]$lam_H[, 1, 3] * chains_m4[[i]]$sigma_H[, 1],
               c(0.5, 0.025, 0.975))
    lam_H_3_2_m4[[i]] <-
      quantile(chains_m4[[i]]$lam_H[, 2, 3] * chains_m4[[i]]$sigma_H[, 2],
               c(0.5, 0.025, 0.975))
    lam_H_3_3_m4[[i]] <-
      quantile(chains_m4[[i]]$lam_H[, 3, 3] * chains_m4[[i]]$sigma_H[, 3],
               c(0.5, 0.025, 0.975))
    lam_H_3_4_m4[[i]] <-
      quantile(chains_m4[[i]]$lam_H[, 4, 3] * chains_m4[[i]]$sigma_H[, 4],
               c(0.5, 0.025, 0.975))
    lam_H_3_5_m4[[i]] <-
      quantile(chains_m4[[i]]$lam_H[, 5, 3] * chains_m4[[i]]$sigma_H[, 5],
               c(0.5, 0.025, 0.975))
    if (i %in% c(5, 12, 20, 21, 23, 27, 30, 31)) {
      lam_H_4_1_m4[[i]] <-
        quantile(chains_m4[[i]]$lam_H[, 1, 4] * chains_m4[[i]]$sigma_H[, 1],
                 c(0.5, 0.025, 0.975))
      lam_H_4_2_m4[[i]] <-
        quantile(chains_m4[[i]]$lam_H[, 2, 4] * chains_m4[[i]]$sigma_H[, 2],
                 c(0.5, 0.025, 0.975))
      lam_H_4_3_m4[[i]] <-
        quantile(chains_m4[[i]]$lam_H[, 3, 4] * chains_m4[[i]]$sigma_H[, 3],
                 c(0.5, 0.025, 0.975))
      lam_H_4_4_m4[[i]] <-
        quantile(chains_m4[[i]]$lam_H[, 4, 4] * chains_m4[[i]]$sigma_H[, 5],
                 c(0.5, 0.025, 0.975))
      lam_H_4_5_m4[[i]] <-
        quantile(chains_m4[[i]]$lam_H[, 5, 4] * chains_m4[[i]]$sigma_H[, 5],
                 c(0.5, 0.025, 0.975))
    } else{
      lam_H_4_1_m4[[i]] <- rep(NA, 3)
      lam_H_4_2_m4[[i]] <- rep(NA, 3)
      lam_H_4_3_m4[[i]] <- rep(NA, 3)
      lam_H_4_4_m4[[i]] <- rep(NA, 3)
      lam_H_4_5_m4[[i]] <- rep(NA, 3)
    }
  }
  
  pars_m4[[i]] <- rbind(
    rho_m4[[i]],
    gamma_m4[[i]],
    chi_m4[[i]],
    lam_H_1_1_m4[[i]],
    lam_H_2_1_m4[[i]],
    lam_H_3_1_m4[[i]],
    lam_H_4_1_m4[[i]],
    lam_H_1_2_m4[[i]],
    lam_H_2_2_m4[[i]],
    lam_H_3_2_m4[[i]],
    lam_H_4_2_m4[[i]],
    lam_H_1_3_m4[[i]],
    lam_H_2_3_m4[[i]],
    lam_H_3_3_m4[[i]],
    lam_H_4_3_m4[[i]],
    lam_H_1_4_m4[[i]],
    lam_H_2_4_m4[[i]],
    lam_H_3_4_m4[[i]],
    lam_H_4_4_m4[[i]],
    lam_H_1_5_m4[[i]],
    lam_H_2_5_m4[[i]],
    lam_H_3_5_m4[[i]],
    lam_H_4_5_m4[[i]]
  )
  pars_m4[[i]] <- as.data.frame(pars_m4[[i]])
  pars_m4[[i]]$pars <- c(
    "rho",
    "gamma",
    "chi",
    "lam_H_1_1",
    "lam_H_2_1",
    "lam_H_3_1",
    "lam_H_4_1",
    "lam_H_1_2",
    "lam_H_2_2",
    "lam_H_3_2",
    "lam_H_4_2",
    "lam_H_1_3",
    "lam_H_2_3",
    "lam_H_3_3",
    "lam_H_4_3",
    "lam_H_1_4",
    "lam_H_2_4",
    "lam_H_3_4",
    "lam_H_4_4",
    "lam_H_1_5",
    "lam_H_2_5",
    "lam_H_3_5",
    "lam_H_4_5"
  )
  pars_m4[[i]]$ENTIDAD_RES <- rep(as.character(mex_states_order[i, ]), 23)
  colnames(pars_m4[[i]]) <- c("med", "ciL", "ciU", "pars", "ENTIDAD_RES")
  
}

pars_df_m4 <- bind_rows(pars_m4)



lam_1_m4 <- vector(mode = "list", length = 32)
lam_2_m4 <- vector(mode = "list", length = 32)
lam_3_m4 <- vector(mode = "list", length = 32)
lam_4_m4 <- vector(mode = "list", length = 32)
FOI_m4 <- vector(mode = "list", length = 32)

for (i in 5:5) {
#for(i in seq(1, 31)[-c(2, 6, 9, 10, 29)]){ this is all states
  
  dim(chains_m4[[i]]$lam_t)[2]
  
  
  lam_1_m4[[i]] <- data.frame(
    t = seq(1, end_t),
    med = NA,
    ciL = NA,
    ciU = NA,
    ENTIDAD_RES = mex_states_order[i, ]
  )
  lam_2_m4[[i]] <-  lam_1_m4[[i]]
  lam_3_m4[[i]] <-  lam_1_m4[[i]]
  lam_4_m4[[i]] <-  lam_1_m4[[i]]
  FOI_m4[[i]] <-  lam_1_m4[[i]]
  
  lam_1_m4[[i]][t, 2:4] <-
    quantile(chains_m4[[i]]$lam_t[, t, 1] * chains_m4[[i]]$sigma_t[, t], c(0.5, 0.025, 0.975))
  lam_2_m4[[i]][t, 2:4] <-
    quantile(chains_m4[[i]]$lam_t[, t, 2] * chains_m4[[i]]$sigma_t[, t], c(0.5, 0.025, 0.975))
  
  
  for (t in 1:end_t) {
    if (i %in% c(11, 18, 24)) {
      lam_3_m4[[i]][t, 2:4] <- rep(NA, 3)
      lam_4_m4[[i]][t, 2:4] <- rep(NA, 3)
      FOI_m4[[i]][t, 2:4] <-
        quantile((chains_m4[[i]]$lam_t[, t, 1] + chains_m4[[i]]$lam_t[, t, 2]) *
                   chains_m4[[i]]$sigma_t[, t],
                 c(0.5, 0.025, 0.975)
        )
    }
    else if (i %in% c(5, 12, 20, 21, 23, 27, 30, 31)) {
      lam_3_m4[[i]][t, 2:4] <-
        quantile(chains_m4[[i]]$lam_t[, t, 3] * chains_m4[[i]]$sigma_t[, t],
                 c(0.5, 0.025, 0.975))
      lam_4_m4[[i]][t, 2:4] <-
        quantile(chains_m4[[i]]$lam_t[, t, 4] * chains_m4[[i]]$sigma_t[, t],
                 c(0.5, 0.025, 0.975))
      FOI_m4[[i]][t, 2:4] <-
        quantile((
          chains_m4[[i]]$lam_t[, t, 1] + chains_m4[[i]]$lam_t[, t, 2] +
            chains_m4[[i]]$lam_t[, t, 3] + chains_m4[[i]]$lam_t[, t, 4]
        ) * chains_m4[[i]]$sigma_t[, t],
        c(0.5, 0.025, 0.975)
        )
    } else{
      lam_3_m4[[i]][t, 2:4] <-
        quantile(chains_m4[[i]]$lam_t[, t, 3] * chains_m4[[i]]$sigma_t[, t],
                 c(0.5, 0.025, 0.975))
      lam_4_m4[[i]][t, 2:4] <- rep(NA, 3)
      FOI_m4[[i]][t, 2:4] <-
        quantile((
          chains_m4[[i]]$lam_t[, t, 1] + chains_m4[[i]]$lam_t[, t, 2] +
            chains_m4[[i]]$lam_t[, t, 3]
        ) * chains_m4[[i]]$sigma_t[, t],
        c(0.5, 0.025, 0.975)
        )
    }
  }
  
  
  if (i %in% c(1, 9, 13, 22)) {
    lam_1_m4[[i]]$type <- c("T8")
  } else if (i %in% c(3, 23, 31)) {
    lam_1_m4[[i]]$type <- c("T7", "T8")
  } else if (i == 10) {
    lam_1_m4[[i]]$type <- c("T1", "T2")
  } else if (i == 26) {
    lam_1_m4[[i]]$type <- c("T6", "T7", "T8")
  } else if (i == 18) {
    lam_1_m4[[i]]$type <-
      c("T3", "T4", "T5", "T6")
  } else if (i == 24) {
    lam_1_m4[[i]]$type <-
      c("T2", "T3", "T4", "T5", "T6")
  } else if (i == 4) {
    lam_1_m4[[i]]$type <-
      c("T4", "T5", "T6", "T7", "T8")
  } else if (i == 28) {
    lam_1_m4[[i]]$type <-
      c("T1", "T2", "T3", "T4", "T5")
  } else if (i %in% c(8, 25)) {
    lam_1_m4[[i]]$type <-
      c("T3", "T4", "T5", "T6", "T7", "T8")
  } else if (i == 11) {
    lam_1_m4[[i]]$type <-
      c("T1", "T2", "T3", "T4", "T5", "T6")
  } else if (i == 7) {
    lam_1_m4[[i]]$type <-
      c("T2", "T3", "T4", "T5", "T6", "T7", "T8")
  } else{
    lam_1_m4[[i]]$type <-
      c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8")
  }
  
  lam_2_m4[[i]]$type <-  lam_1_m4[[i]]$type
  lam_3_m4[[i]]$type <-  lam_1_m4[[i]]$type
  lam_4_m4[[i]]$type <-  lam_1_m4[[i]]$type
  FOI_m4[[i]]$type <-  lam_1_m4[[i]]$type
}

lam_df_m4 <- bind_rows(
  bind_rows(lam_1_m4),
  bind_rows(lam_2_m4),
  bind_rows(lam_3_m4),
  bind_rows(lam_4_m4),
  .id = "id"
) %>% rename(serotype = id) 
