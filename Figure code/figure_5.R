library(tidyverse)
library(MetBrewer)
library(sf)

#Download from GADM
mex_geo_adm_lev1 <- st_read("gadm41_MEX.gpkg", layer = "ADM_ADM_1")
mex_geo_adm_lev0 <- st_read("gadm41_MEX.gpkg", layer = "ADM_ADM_0")

lam_df_m4 <- readRDS("lam_df_m4.rds")
pars_df_m4 <- readRDS("pars_df_m4.rds")

int_m4 <- bind_rows(lam_df_m4 %>%
                      mutate(year = as.factor(as.numeric(gsub("T", "", type)) + 2015)), 
                    pars_df_m4 %>% 
                      filter(pars %in% c("lam_H_1_1", "lam_H_2_1", "lam_H_3_1", "lam_H_4_1")) %>%
                      mutate(pars = gsub("lam_H_", "", pars)) %>%
                      mutate(pars = ifelse(
                        pars %in% c(1, 2, 3, 4), paste0(pars, "_constant"), pars
                      )) %>%
                      separate_wider_delim(pars, delim = "_", names = c("serotype", "type")) %>%
                      mutate(type = paste0("lam_H_", type)) %>%
                      mutate(type = "TH", year = "pre-2016")) %>%
  mutate(serotype = paste0("DENV-", serotype)) %>%
  select(c(serotype, med, ENTIDAD_RES, year)) %>%
  pivot_wider(names_from = "year", values_from = med) %>%
  mutate(`2022` = ifelse(is.na(`2022`) & is.na(`2021`) & is.na(`2020`) & is.na(`2019`) &
                           is.na(`2018`) & is.na(`2017`) & is.na(`2016`), `pre-2016`, `2022`),
         `2021` = ifelse(is.na(`2021`) & is.na(`2020`) & is.na(`2019`) &
                           is.na(`2018`) & is.na(`2017`) & is.na(`2016`), `pre-2016`, `2021`),
         `2020` = ifelse(is.na(`2020`) & is.na(`2019`) & is.na(`2018`) & is.na(`2017`) & 
                           is.na(`2016`), `pre-2016`, `2020`),
         `2019` = ifelse(is.na(`2019`) & is.na(`2018`) & is.na(`2017`) & is.na(`2016`), 
                         `pre-2016`, `2019`),
         `2018` = ifelse(is.na(`2018`) & is.na(`2017`) & is.na(`2016`), `pre-2016`, `2018`),
         `2017` = ifelse(is.na(`2017`) & is.na(`2016`), `pre-2016`, `2017`),
         `2016` = ifelse(is.na(`2016`), `pre-2016`, `2016`)) %>%
  pivot_longer(cols = !c(serotype, ENTIDAD_RES), names_to = "year", values_to = "med") %>%
  left_join(lam_df_m4 %>%
              mutate(year = as.factor(as.numeric(gsub("T", "", type)) + 2015),
                     serotype = paste0("DENV-", serotype)) %>%
              select(ENTIDAD_RES, serotype, year, med) %>% 
              rename(med_old = med),
            by = join_by(ENTIDAD_RES, serotype, year)) %>%
  mutate(historic = ifelse(year == "pre-2016", TRUE, FALSE)) %>%
  mutate(historic = ifelse(is.na(med_old) & !is.na(med), TRUE, historic)) %>%
  mutate(med_quan = cut(med, c(0, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5))) %>%
  mutate(year = factor(year, levels = c("pre-2016", 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023)))

map_params_m4 <- left_join(
  mex_geo_adm_lev1 %>%
    mutate(NAME_1 = iconv(tolower(NAME_1),
                          to = "ASCII//TRANSLIT")), 
  int_m4 %>% filter(year != "pre-2016"),
  by = join_by(NAME_1 == ENTIDAD_RES)
) 

ggplot() +
  geom_sf(
    data = map_params_m4 %>%
      filter(!is.na(med), ! NAME_1 %in% c("zacatecas", "tlaxcala", "chihuahua","baja california", "distrito federal",
                                          "durango", "yucatan")),
    aes(fill = serotype,
        alpha = med_quan),
    colour = "white",
    lwd = 0.1
  ) +
  geom_sf(
    data = map_params_m4 %>% filter(is.na(med)) %>%
      filter(! NAME_1 %in% c("zacatecas", "tlaxcala", "chihuahua","baja california", "distrito federal",
                             "durango", "yucatan")),
    fill = "grey50",
    colour = "white",
    lwd = 0.1
  ) +
  geom_sf(
    data = mex_geo_adm_lev1 %>%
      mutate(NAME_1 = iconv(tolower(NAME_1),
                            to = "ASCII//TRANSLIT")) %>% 
      filter(NAME_1 %in% c("zacatecas", "tlaxcala", "chihuahua","baja california", "distrito federal",
                           "durango", "yucatan")),
    fill = "grey50",
    colour = "white",
    lwd = 0.1
  ) +
  geom_sf(data = mex_geo_adm_lev0,
          fill = NA,
          colour = "grey30",
          lwd = 0.2) +
  facet_grid(year~serotype, axes = "all") +
  theme_classic() +
  MetBrewer::scale_fill_met_d("Egypt") +
  theme(text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.y = element_text(size = 10, angle = 0),
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill = NA, color = "white"),
        legend.title = element_text(size = 10),
  ) +
  labs(fill = "Serotype",
       alpha = "FOI") + 
  guides(x = "none", 
         y = "none")


ggsave(
  filename = "figure_5.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 190,
  height = 230,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  create.dir = FALSE
)