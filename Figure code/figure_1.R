library(tidyverse)
library(MetBrewer)
library(patchwork)
library(binom)

mex_dengue_data_2016_2023 <- readRDS("mex_dengue_data_2016_2023.rds")
mex_pop_growth_rates <- readRDS("mex_pop_growth_rates.rds")

mex_inc_ot <- bind_cols(
  mex_dengue_data_2016_2023 %>%
    filter(ESTATUS_CASO != 3) %>%
    filter(
      ENTIDAD_RES != "otros paises" &
        ENTIDAD_RES != "otros paises de latinoamerica" &
        ENTIDAD_RES != "estados unidos de america"
    ) %>%
    count(YEAR) %>%
    group_by(YEAR) %>%
    rename(cases = n) %>%
    arrange(YEAR),
  mex_pop_growth_rates %>%
    filter(NOM_ENT == "estados unidos mexicanos", Sex == "Total", Age == "TOT") %>%
    select(!c(NOM_ENT, Sex, growth, y_2015)) %>%
    pivot_longer(
      names_to = "Year",
      values_to = "Pop",
      cols = 2:9,
      names_prefix = "y_"
    ) %>%
    mutate(Year = as.numeric(Year)) %>%
    arrange(Year)
) %>% mutate(inc = cases / Pop)


mex_prop_hosp_ot <- left_join(
  mex_dengue_data_2016_2023 %>%
    filter(ESTATUS_CASO != 3) %>%
    filter(is.na(TIPO_PACIENTE) == FALSE) %>%
    filter(
      ENTIDAD_RES != "otros paises" &
        ENTIDAD_RES != "otros paises de latinoamerica" &
        ENTIDAD_RES != "estados unidos de america"
    ) %>%
    count(YEAR, TIPO_PACIENTE) %>%
    group_by(YEAR) %>%
    mutate(prop = n / sum(n)) %>%
    filter(TIPO_PACIENTE == 2) %>%
    select(!TIPO_PACIENTE),
  mex_dengue_data_2016_2023 %>%
    filter(ESTATUS_CASO != 3) %>%
    filter(is.na(TIPO_PACIENTE) == FALSE) %>%
    filter(
      ENTIDAD_RES != "otros paises" &
        ENTIDAD_RES != "otros paises de latinoamerica" &
        ENTIDAD_RES != "estados unidos de america"
    ) %>%
    count(YEAR),
  by = join_by(YEAR)
) %>%
  rename(hosp = n.x, tot = n.y) %>%
  mutate(
    prop.lb = binom.confint(hosp, tot, methods = "exact")$lower,
    prop.ub = binom.confint(hosp, tot, methods = "exact")$upper
  )



mex_CFR_ot <- left_join(
  mex_dengue_data_2016_2023 %>%
    filter(ESTATUS_CASO != 3) %>%
    filter(
      ENTIDAD_RES != "otros paises" &
        ENTIDAD_RES != "otros paises de latinoamerica" &
        ENTIDAD_RES != "estados unidos de america"
    ) %>%
    count(YEAR, DEFUNCION) %>%
    group_by(YEAR) %>%
    mutate(prop = n / sum(n)) %>%
    filter(DEFUNCION == 1) %>%
    select(!DEFUNCION),
  mex_dengue_data_2016_2023 %>%
    filter(ESTATUS_CASO != 3) %>%
    filter(
      ENTIDAD_RES != "otros paises" &
        ENTIDAD_RES != "otros paises de latinoamerica" &
        ENTIDAD_RES != "estados unidos de america"
    ) %>%
    count(YEAR),
  by = join_by(YEAR)
) %>%
  rename(dec = n.x, tot = n.y) %>%
  mutate(
    prop.lb = binom.confint(dec, tot, methods = "exact")$lower,
    prop.ub = binom.confint(dec, tot, methods = "exact")$upper
  )


sero_plot_mex <- mex_dengue_data_2016_2023  %>%
  filter(ESTATUS_CASO != 3) %>%
  filter(RESULTADO_PCR != 5) %>%
  filter(
    ENTIDAD_RES != "otros paises" &
      ENTIDAD_RES != "estados unidos de america" &
      ENTIDAD_RES != "otros paises de latinoamerica"
  ) %>%
  count(RESULTADO_PCR, YEAR) %>%
  complete(RESULTADO_PCR, nesting(YEAR), fill = list(n = 0)) %>%
  complete(YEAR, nesting(RESULTADO_PCR), fill = list(n = 0)) %>%
  group_by(YEAR) %>%
  mutate(prop = n / sum(n)) %>%
  mutate(YEAR = as.numeric(YEAR))

mex_inc <- bind_cols(
  mex_dengue_data_2016_2023 %>%
    filter(ESTATUS_CASO != 3) %>%
    filter(
      ENTIDAD_RES != "otros paises" &
        ENTIDAD_RES != "otros paises de latinoamerica" &
        ENTIDAD_RES != "estados unidos de america"
    ) %>%
    count(YEAR, AGE_BAND) %>%
    rename(cases = n) %>%
    mutate(AGE_BAND = factor(
      AGE_BAND,
      levels = c(
        "[0,5)",
        "[5,10)",
        "[10,15)",
        "[15,20)",
        "[20,25)",
        "[25,30)",
        "[30,35)",
        "[35,40)",
        "[40,45)",
        "[45,50)",
        "[50,55)",
        "[55,60)",
        "[60,65)",
        "[65,70)",
        "[70,75)",
        "[75,120)"
      )
    )) %>%
    arrange(YEAR, AGE_BAND),
  mex_pop_growth_rates %>%
    filter(NOM_ENT == "estados unidos mexicanos", Sex == "Total", Age != "TOT") %>%
    select(!c(NOM_ENT, Sex, growth, y_2015)) %>%
    pivot_longer(
      names_to = "Year",
      values_to = "Pop",
      cols = 2:9,
      names_prefix = "y_"
    ) %>%
    mutate(Year = as.numeric(Year)) %>%
    arrange(Year)
) %>%
  mutate(inc = cases / Pop,
         year_2 = ifelse(
           YEAR == 2016,
           "Age Group \n \n 2016",
           ifelse(
             YEAR == 2017,
             "Age Group \n \n 2017",
             ifelse(
               YEAR == 2018,
               "Age Group \n \n 2018",
               ifelse(
                 YEAR == 2019,
                 "Age Group \n \n 2019",
                 ifelse(
                   YEAR == 2020,
                   "Age Group \n \n 2020",
                   ifelse(
                     YEAR == 2021,
                     "Age Group \n \n 2021",
                     ifelse(YEAR == 2022, "Age Group \n \n 2022", "Age Group \n \n 2023")
                   )
                 )
               )
             )
           )
         ))

mex_prop_hosp_ns <- mex_dengue_data_2016_2023 %>%
  filter(ESTATUS_CASO != 3) %>%
  filter(is.na(TIPO_PACIENTE) == FALSE) %>%
  filter(
    ENTIDAD_RES != "otros paises" &
      ENTIDAD_RES != "otros paises de latinoamerica" &
      ENTIDAD_RES != "estados unidos de america"
  ) %>%
  count(YEAR, AGE_BAND, TIPO_PACIENTE) %>%
  group_by(YEAR, AGE_BAND) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  filter(TIPO_PACIENTE == 2) %>%
  mutate(
    AGE_BAND = factor(
      AGE_BAND,
      levels = c(
        "[0,5)",
        "[5,10)",
        "[10,15)",
        "[15,20)",
        "[20,25)",
        "[25,30)",
        "[30,35)",
        "[35,40)",
        "[40,45)",
        "[45,50)",
        "[50,55)",
        "[55,60)",
        "[60,65)",
        "[65,70)",
        "[70,75)",
        "[75,120)"
      )
    ),
    year_2 = ifelse(
      YEAR == 2016,
      "Age Group \n \n 2016",
      ifelse(
        YEAR == 2017,
        "Age Group \n \n 2017",
        ifelse(
          YEAR == 2018,
          "Age Group \n \n 2018",
          ifelse(
            YEAR == 2019,
            "Age Group \n \n 2019",
            ifelse(
              YEAR == 2020,
              "Age Group \n \n 2020",
              ifelse(
                YEAR == 2021,
                "Age Group \n \n 2021",
                ifelse(YEAR == 2022, "Age Group \n \n 2022", "Age Group \n \n 2023")
              )
            )
          )
        )
      )
    )
  )

mex_CFR_ns <- mex_dengue_data_2016_2023 %>%
  filter(ESTATUS_CASO != 3) %>%
  filter(
    ENTIDAD_RES != "otros paises" &
      ENTIDAD_RES != "otros paises de latinoamerica" &
      ENTIDAD_RES != "estados unidos de america"
  ) %>%
  count(YEAR, AGE_BAND, DEFUNCION) %>%
  group_by(YEAR, AGE_BAND) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  filter(DEFUNCION == 1) %>%
  mutate(
    AGE_BAND = factor(
      AGE_BAND,
      levels = c(
        "[0,5)",
        "[5,10)",
        "[10,15)",
        "[15,20)",
        "[20,25)",
        "[25,30)",
        "[30,35)",
        "[35,40)",
        "[40,45)",
        "[45,50)",
        "[50,55)",
        "[55,60)",
        "[60,65)",
        "[65,70)",
        "[70,75)",
        "[75,120)"
      )
    ),
    year_2 = ifelse(
      YEAR == 2016,
      "Age Group \n \n 2016",
      ifelse(
        YEAR == 2017,
        "Age Group \n \n 2017",
        ifelse(
          YEAR == 2018,
          "Age Group \n \n 2018",
          ifelse(
            YEAR == 2019,
            "Age Group \n \n 2019",
            ifelse(
              YEAR == 2020,
              "Age Group \n \n 2020",
              ifelse(
                YEAR == 2021,
                "Age Group \n \n 2021",
                ifelse(YEAR == 2022, "Age Group \n \n 2022", "Age Group \n \n 2023")
              )
            )
          )
        )
      )
    )
  )


inc_short <- ggplot(aes(x = YEAR, y = inc), data = mex_inc_ot) +
  geom_line(aes(group = Age), colour = "#808fe1", size = 0.5) +
  geom_point(colour = "#808fe1", size = 0.7) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 10),
    strip.background = element_rect(fill = NA, color = "white"),
    panel.spacing = unit(0.1, "cm"),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  ylim(0, 0.002) +
  labs(y = "Incidence", x = "Year")

hosp_short <- ggplot(data = mex_prop_hosp_ot) +
  geom_point(aes(x = YEAR, y = prop),
             colour = "#97c684",
             size = 0.7) +
  geom_errorbar(
    aes(ymin = prop.lb, ymax = prop.ub, x = YEAR),
    colour = "#97c684",
    width = 0.3
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 10),
    strip.background = element_rect(fill = NA, color = "white"),
    panel.spacing = unit(0.1, "cm"),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  ylim(0, 0.22) +
  labs(y = "Proportion \n Hospitalised", x = "Year")

cfr_short <- ggplot(data = mex_CFR_ot) +
  geom_errorbar(
    aes(ymin = prop.lb, ymax = prop.ub, x = YEAR),
    colour = "#5c66a8",
    width = 0.3
  ) +
  geom_point(aes(x = YEAR, y = prop),
             colour = "#5c66a8",
             size = 0.7) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 10),
    strip.background = element_rect(fill = NA, color = "white"),
    panel.spacing = unit(0.1, "cm"),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  ylim(0, 0.0041) +
  labs(y = "CFR", x = "Year")

sero_short <- ggplot(aes(
  x = YEAR,
  y = prop,
  colour = as.factor(RESULTADO_PCR)
),
data = sero_plot_mex %>% mutate(RESULTADO_PCR = paste0("DENV", RESULTADO_PCR))) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = NA, color = "white"),
    panel.spacing = unit(5, "mm"),
    text = element_text(size = 10),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  ylim(0, 1) +
  labs(colour = "Serotype", y = "Proportion \n PCR Positive", x = "Year") +
  scale_x_continuous(breaks = c(2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023)) +
  scale_colour_manual(values = met.brewer("Egypt", 4))



inc_long <- ggplot(data = mex_inc) +
  geom_col(aes(y = inc, x = AGE_BAND),
           fill = "#808fe1",
           position = "dodge") +
  facet_grid(~ year_2, switch = "x") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = NA, color = "white"),
    text = element_text(size = 10),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    legend.background = element_blank(),
    panel.spacing = unit(0.1, "cm"),
    panel.border = element_blank()
  ) +
  scale_x_discrete(breaks = c("[0,5)", "[75,120)")) +
  labs(x = "Year", y = "Incidence") +
  ylim(0, 0.003)


ph_long <- ggplot(data = mex_prop_hosp_ns) +
  geom_col(aes(y = prop, x = AGE_BAND),
           fill = "#97c684",
           position = "dodge") +
  facet_grid(~ year_2, switch = "x") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = NA, color = "white"),
    text = element_text(size = 10),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    legend.position = "none",
    legend.background = element_blank(),
    panel.spacing = unit(0.1, "cm"),
    panel.border = element_blank()
  ) +
  scale_x_discrete(breaks = c("[0,5)", "[75,120)")) +
  labs(x = "Year", y = "Proportion Hospitalised") +
  ylim(0, 0.4)


cfr_long <- ggplot(data = mex_CFR_ns) +
  geom_col(aes(y = prop, x = AGE_BAND),
           fill = "#5c66a8",
           position = "dodge") +
  facet_grid(~ year_2, switch = "x") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = NA, color = "white"),
    text = element_text(size = 10),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    legend.position = "none",
    legend.background = element_blank(),
    panel.spacing = unit(0.1, "cm"),
    panel.border = element_blank()
  ) +
  scale_x_discrete(breaks = c("[0,5)", "[75,120)")) +
  labs(x = "Year", y = "CFR") +
  ylim(0, 0.04)


f1_layout  <- c(
  area(t = 1, l = 1),
  area(t = 1, l = 2),
  area(t = 1, l = 3),
  area(t = 1, l = 4, r = 5),
  area(t = 2, l = 1, r = 5),
  area(t = 3, l = 1, r = 5),
  area(t = 4, l = 1, r = 5)
)


free(inc_short) + free(hosp_short) + free(cfr_short) + free(sero_short) + inc_long + ph_long + cfr_long +
  plot_layout(design = f1_layout) +
  plot_annotation(tag_levels = 'a',
                  theme = theme(plot.tag = element_text(size = 10)))


ggsave(
  filename = "figure_1.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 190,
  height = 200,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  create.dir = FALSE
)