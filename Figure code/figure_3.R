library(tidyverse)
library(MetBrewer)
library(sf)
library(patchwork)

mex_geo_adm_lev1 <- st_read("gadm41_MEX.gpkg", layer = "ADM_ADM_1")
mex_geo_adm_lev0 <- st_read("gadm41_MEX.gpkg", layer = "ADM_ADM_0")

mex_dengue_data_2016_2023 <- readRDS("mex_dengue_data_2016_2023.rds")
mex_pop_growth_rates <- readRDS("mex_pop_growth_rates.rds")
pars_df_m1 <- readRDS("pars_df_m1.rds")

e_cases_m1 <- vector(mode = "list", length = 32)

for (i in seq(1, 31)[-c(2, 6, 29)]) {
  e_cases_m1[[i]] <-
    quantile(rowSums((chains_m1[[i]]$inc1 + chains_m1[[i]]$inc2)/
                       matrix(rep((age_max_m1 - age_min_m1), 15000), ncol = 16, byrow = TRUE) * 
                       matrix(rep(colSums(mex_pop_mat_no0[((i-1)*8+1):(i*8),]), 15000), ncol = 16, byrow = TRUE)), 
             c(0.5, 0.025, 0.975))
  
  e_cases_m1[[i]]$ENTIDAD_RES <- as.character(mex_states_order[i,])
  e_cases_m1[[i]]$total_cases <- sum(mex_cases_mat_no0[((i-1)*8+1):(i*8),])
}

e_cases_m1_df <- bind_rows(e_cases_m1) %>%
  mutate(sens = total_cases/`50%`,
         sens_ciL = total_cases/`2.5%`,
         sens_ciU = total_cases/`97.5%`)

map_RR_agg <- left_join(mex_geo_adm_lev1 %>%
                          mutate(NAME_1 = iconv(tolower(NAME_1),
                                                to = "ASCII//TRANSLIT")),
                        left_join(mex_dengue_data_2016_2023 %>%
                                    filter(ESTATUS_CASO != 3) %>%
                                    count(ESTATUS_CASO, ENTIDAD_RES) %>% 
                                    group_by(ENTIDAD_RES) %>%
                                    mutate(total_cases = sum(n)) %>% 
                                    ungroup() %>% 
                                    distinct(ENTIDAD_RES, total_cases), 
                                  mex_pop_growth_rates %>%
                                    filter(NOM_ENT != "estados unidos mexicanos" & Age == "TOT") %>%
                                    mutate(total_pop = y_2016 + y_2017 + y_2018 + y_2019 + y_2020 + y_2021 + y_2022 + y_2023, 
                                           NOM_ENT = gsub("'", "", NOM_ENT)) %>%
                                    select(NOM_ENT, total_pop) %>%
                                    mutate(NOM_ENT = ifelse(NOM_ENT == "ciudad de mexico", "distrito federal", 
                                                            ifelse(NOM_ENT == "michoacan de ocampo", "michoacan", NOM_ENT))),
                                  by = join_by(ENTIDAD_RES == NOM_ENT)) %>%
                          mutate(RR_agg = total_cases/total_pop),
                        by = join_by(NAME_1 == ENTIDAD_RES))

map_age <- left_join(
  mex_geo_adm_lev1 %>%
    mutate(NAME_1 = iconv(tolower(NAME_1),
                          to = "ASCII//TRANSLIT")),
  mex_dengue_data_2016_2023 %>%
    filter(ESTATUS_CASO != 3) %>%
    group_by(ENTIDAD_RES) %>%
    summarise(mean = mean(EDAD_ANOS), median = median(EDAD_ANOS)),
  by = join_by(NAME_1 == ENTIDAD_RES)
)


map_params_m1 <- left_join(
  mex_geo_adm_lev1 %>%
    mutate(NAME_1 = iconv(tolower(NAME_1),
                          to = "ASCII//TRANSLIT")),
  pars_df_m1 %>% 
    pivot_wider(names_from = pars, values_from = c(med, ciL, ciU)),
  by = join_by(NAME_1 == ENTIDAD_RES)
) 

map_params_m1_sens <- left_join(
  mex_geo_adm_lev1 %>%
    mutate(NAME_1 = iconv(tolower(NAME_1),
                          to = "ASCII//TRANSLIT")),
  e_cases_m1_df,
  by = join_by(NAME_1 == ENTIDAD_RES)
)


f3_1 <- ggplot() +
  geom_sf(
    data = map_RR_agg,
    aes(fill = RR_agg),
    colour = "white",
    lwd = 0.3
  )  +
  geom_sf(data = mex_geo_adm_lev0,
          fill= NA,
          colour = "grey30",
          lwd = 0.3) +
  theme_classic() +
  scale_fill_gradientn(colors = met.brewer("OKeeffe2"), limits = c(0,0.004)) +
  theme(legend.position = "right",
        text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_rect(fill = NA, color = "white"),
  ) +
  labs(fill = "Annual \n incidence") + 
  guides(x = "none", y = "none")


f3_2 <- ggplot() +
  geom_sf(
    data = map_age,
    aes(fill = median),
    colour = "white",
    lwd = 0.3
  ) +
  geom_sf(data = mex_geo_adm_lev0,
          fill= NA,
          colour = "grey30",
          lwd = 0.3) +
  theme_classic() +
  scale_fill_gradientn(colors = met.brewer("OKeeffe2", direction = -1), limits = c(14, 35)) +
  theme(legend.position = "right",
        text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_rect(fill = NA, color = "white"),
  ) +
  labs(fill = "Median \n age") + 
  guides(x = "none", y = "none")

map_params_m1 <- left_join(
  mex_geo_adm_lev1 %>%
    mutate(NAME_1 = iconv(tolower(NAME_1),
                          to = "ASCII//TRANSLIT")),
  pars_df_m1 %>% 
    pivot_wider(names_from = pars, values_from = c(med, ciL, ciU)),
  by = join_by(NAME_1 == ENTIDAD_RES)
) 


f3_3 <- ggplot() +
  geom_sf(data = mex_geo_adm_lev1,
          fill = "grey50",
          colour = "white",
          lwd = 0.3) +
  geom_sf(
    data = map_params_m1 %>% filter(NAME_1 != "zacatecas" & NAME_1 != "tlaxcala" & NAME_1 != "chihuahua" & 
                                        NAME_1 != "baja california" & NAME_1 != "distrito federal"),
    aes(fill = med_lam),
    colour = "white",
    lwd = 0.3
  ) +
  geom_sf(data = mex_geo_adm_lev0,
          fill= NA,
          colour = "grey30",
          lwd = 0.3) +
  theme_classic() +
  scale_fill_gradientn(colors = met.brewer("OKeeffe2"), limits = c(0.0075,0.02)) +
  theme(legend.position = "right",
        text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_rect(fill = NA, color = "white"),
        legend.title = element_text(size = 10),
  ) + 
  labs(fill = "FOI") +
  guides(x = "none", y = "none")


f3_4 <- ggplot() +
  geom_sf(data = mex_geo_adm_lev1,
          fill = "grey50",
          colour = "white",
          lwd = 0.3) +
  geom_sf(
    data = map_params_m1_sens %>% filter(NAME_1 != "zacatecas" & NAME_1 != "tlaxcala" & NAME_1 != "chihuahua" & 
                                             NAME_1 != "baja california" & NAME_1 != "distrito federal"),
    aes(fill = sens),
    colour = "white",
    lwd = 0.3
  ) +
  geom_sf(data = mex_geo_adm_lev0,
          fill= NA,
          colour = "grey30",
          lwd = 0.3) +
  theme_classic() +
  scale_fill_gradientn(colors = met.brewer("OKeeffe2"), limits = c(0,0.16), breaks = c(0, 0.04, 0.08, 0.12, 0.16)) +
  theme(legend.position = "right",
        text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_rect(fill = NA, color = "white"),
        legend.title = element_text(size = 10),
  ) + 
  labs(fill = "Sensitivity") +
  guides(x = "none", y = "none")


(f3_1 | f3_2) / (f3_3 | f3_4) +
  plot_annotation(
    tag_levels = 'a',
    theme = theme(plot.tag = element_text(size = 10)))

ggsave(
  filename = "figure_3.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 190,
  height = 120,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  create.dir = FALSE
)