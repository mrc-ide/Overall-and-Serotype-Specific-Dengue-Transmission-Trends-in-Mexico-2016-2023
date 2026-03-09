library(tidyverse)
library(MetBrewer)
library(sf)

mex_geo_adm_lev1 <- st_read("gadm41_MEX.gpkg", layer = "ADM_ADM_1")

pars_df_m2 <- readRDS(pars_df_m2, "pars_df_m2.rds")

map_params_m2 <- left_join(
  mex_geo_adm_lev1 %>%
    mutate(NAME_1 = iconv(tolower(NAME_1), to = "ASCII//TRANSLIT")),
  pars_df_m2 %>%
    pivot_wider(names_from = pars, values_from = c(med, ciL, ciU)),
  by = join_by(NAME_1 == ENTIDAD_RES)
) %>%
  left_join(
    lam_df_m2 %>% group_by(ENTIDAD_RES) %>% summarise(mean_med = mean(med)),
    by = join_by(NAME_1 == ENTIDAD_RES)
  ) %>%
  left_join(lam_df_m2, by = join_by(NAME_1 == ENTIDAD_RES))

ggplot() +
  geom_sf(
    data = mex_geo_adm_lev1,
    fill = "grey50",
    colour = "grey10",
    lwd = 0.3
  ) +
  geom_sf(
    data = map_params_m2 %>%
      mutate(year = t + 2015) %>%
      filter(
        NAME_1 != "zacatecas" &
          NAME_1 != "tlaxcala" &
          NAME_1 != "chihuahua" &
          NAME_1 != "baja california" &
          NAME_1 != "distrito federal"
      ),
    aes(fill = med),
    colour = "white",
    lwd = 0.3
  ) +
  facet_wrap( ~ year, nrow = 3, axes = "all") +
  geom_sf(
    data = mex_geo_adm_lev0,
    fill = NA,
    colour = "grey30",
    lwd = 0.3
  ) +
  theme_classic() +
  scale_fill_gradientn(colors = met.brewer("OKeeffe2"), limits = c(0, 0.25)) +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.background = element_rect(fill = NA, color = "white"),
    legend.title = element_text(size = 10),
  ) +
  labs(fill = "FOI")  +
  guides(x = "none", y = "none")


ggsave(
  filename = "figure_4.tiff",
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
