library(tidyverse)
library(MetBrewer)
library(patchwork)

mex_dengue_data_2016_2023 <- readRDS("mex_dengue_data_2016_2023.rds")
lam_df_m4 <- readRDS("lam_df_m4.rds")
pars_df_m4 <- readRDS("pars_df_m4.rds")
fit_m4 <- readRDS("fit_m4.rds")

oax_data_figure_6_p1 <- mex_dengue_data_2016_2023 %>%
  filter(ENTIDAD_RES == "oaxaca") %>%
  filter(ESTATUS_CASO != 3) %>%
  filter(!is.na(TIPO_PACIENTE)) %>%
  count(TIPO_PACIENTE, YEAR, AGE_BAND) %>%
  group_by(TIPO_PACIENTE) %>%
  rename(cases = n)

oax_data_figure_6_p1$AGE_BAND <- factor(
  oax_data_figure_6_p1$AGE_BAND,
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
)


oax_data_figure_6_p1 <- oax_data_figure_6_p1 %>% mutate(year_2 = ifelse(
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


f6_p1 <- ggplot(data = oax_data_figure_6_p1) +
  geom_col(aes(
    y = cases,
    x = AGE_BAND,
    fill = as.factor(TIPO_PACIENTE)
  ), position = "stack") +
  theme_classic() +
  theme(
    strip.placement = "outside",
    axis.text.x = element_blank(),
    legend.position = c(0.15, 0.75),
    legend.background = element_blank(),
    strip.background = element_rect(fill = NA, color = "white"),
    panel.spacing = unit(0.1, "cm"),
    panel.border = element_blank(),
    text = element_text(size = 10)
  ) +
  facet_wrap( ~ year_2, ncol = 8, strip.position = "bottom") +
  labs(x = "Year", y = "Cases") +
  scale_fill_manual(
    values = met.brewer("Derain", 2),
    name = "Hospitalisation Status",
    labels = c("Non-Hospitalised", "Hospitalised")
  ) +
  scale_x_discrete(breaks = c("[0,5)", "[75,120)")) +
  ylim(0, 3000)


oax_data_figure_6_p2 <- mex_dengue_data_2016_2023  %>%
  filter(ESTATUS_CASO != 3) %>%
  filter(RESULTADO_PCR != 5) %>%
  filter(!is.na(TIPO_PACIENTE)) %>%
  count(RESULTADO_PCR, YEAR, TIPO_PACIENTE, ENTIDAD_RES) %>%
  complete(RESULTADO_PCR,
           nesting(TIPO_PACIENTE, YEAR, ENTIDAD_RES),
           fill = list(n = 0)) %>%
  complete(YEAR,
           nesting(TIPO_PACIENTE, RESULTADO_PCR, ENTIDAD_RES),
           fill = list(n = 0)) %>%
  filter(ENTIDAD_RES == "oaxaca") %>%
  group_by(YEAR, TIPO_PACIENTE) %>%
  mutate(prop = n / sum(n)) %>%
  mutate(TIPO_PACIENTE = ifelse(
    TIPO_PACIENTE == 1,
    "Non-Hospitalised",
    ifelse(TIPO_PACIENTE == 2, "Hospitalised", NA)
  ))



f6_p2 <- ggplot(aes(
  x = as.numeric(YEAR),
  y = prop,
  colour = as.factor(RESULTADO_PCR)
), data = oax_data_figure_6_p2) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = NA, color = "white"),
    legend.position = "none",
    panel.spacing = unit(5, "mm"),
    text = element_text(size = 10)
  ) +
  facet_wrap( ~ TIPO_PACIENTE, nrow = 2, scales = "free_x") +
  ylim(0, 1) +
  labs(colour = "Serotype", y = "Proportion PCR Positive", x = "Year") +
  scale_x_continuous(breaks = c(2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023)) +
  scale_colour_manual(values = met.brewer("Egypt", 4))


oax_data_figure_6_p3 <- bind_rows(
  lam_df_m4 %>%
    filter(ENTIDAD_RES == "oaxaca"),
  pars_df_m4 %>%
    filter(
      pars != "rho" &
        pars != "gamma" & pars != "chi",
      ENTIDAD_RES == "oaxaca"
    ) %>%
    mutate(pars = gsub("lam_H_", "", pars)) %>%
    mutate(pars = ifelse(
      pars %in% c(1, 2, 3, 4), paste0(pars, "_constant"), pars
    )) %>%
    separate_wider_delim(pars, delim = "_", names = c("serotype", "type")) %>%
    mutate(type = paste0("lam_H_", type))
) %>%
  mutate(serotype = paste0("DENV-", serotype)) %>%
  mutate(year = ifelse(type == "T1", 2016, ifelse(
    type == "T2", 2017, ifelse(type == "T3", 2018, ifelse(
      type == "T4", 2019, ifelse(type == "T5", "2020", ifelse(
        type == "T6", 2021, ifelse(type == "T7", 2022, ifelse(
          type == "T8",
          2023,
          ifelse(
            type == "lam_H_1",
            "1995-2015",
            ifelse(
              type == "lam_H_2",
              "1975-1995",
              ifelse(
                type == "lam_H_3",
                "1955-1975",
                ifelse(
                  type == "lam_H_4",
                  "1935-1955",
                  ifelse(type == "lam_H_5", "1895-1935", "ERROR")
                )
              )
            )
          )
        ))
      ))
    ))
  ))) %>%
  mutate(year = factor(
    year,
    levels = c(
      "1895-1935",
      "1935-1955",
      "1955-1975",
      "1975-1995",
      "1995-2015",
      2016,
      2017,
      2018,
      2019,
      2020,
      2021,
      2022,
      2023
    )
  ))


f6_p3 <- ggplot(oax_data_figure_6_p3, aes(x = year, y = med, colour = serotype)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(year, ymin = ciL, ymax = ciU),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ylab("Serotype-specific FOI") +
  xlab('Year') +
  scale_y_continuous(
    breaks = c(0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07),
    limits = c(0, 0.07)
  ) +
  scale_colour_manual(values = met.brewer("Egypt", 4), name = "Serotype")


oax_data_figure_6_p4 <- as.data.frame(fit_m4[[20]]) %>%
  mutate(year_2 = ifelse(
    year == 1,
    "Age Group \n \n 2016",
    ifelse(
      year == 2,
      "Age Group \n \n 2017",
      ifelse(
        year == 3,
        "Age Group \n \n 2018",
        ifelse(
          year == 4,
          "Age Group \n \n 2019",
          ifelse(
            year == 5,
            "Age Group \n \n 2020",
            ifelse(
              year == 6,
              "Age Group \n \n 2021",
              ifelse(year == 7, "Age Group \n \n 2022", "Age Group \n \n 2023")
            )
          )
        )
      )
    )
  ))


f6_p4 <- ggplot(oax_data_figure_6_p4, aes(x = Age_Group , y = Cases)) +
  geom_point(aes(col = serotype), size = 0.5) +
  theme_classic() +
  geom_line(aes(
    x = Age_Group,
    y = pred,
    group = interaction(year, serotype),
    col = serotype
  )) +
  facet_wrap(~ year_2, ncol = 8, strip.position = "bottom") +
  geom_ribbon(aes(
    x = Age_Group,
    ymin = ciL,
    ymax = ciU,
    group = interaction(year, serotype),
    fill = serotype
  ),
  alpha = 0.2) +
  scale_fill_met_d(name = "Egypt") +
  scale_color_met_d(name = "Egypt") +
  labs(fill = "DENV Serotype", colour = "DENV Serotype", x = "Year") +
  ylim(0, 2500) +
  theme(
    strip.background = element_rect(fill = NA, color = "white"),
    text = element_text(size = 10),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    legend.position = c(0.1, 0.85),
    legend.background = element_blank(),
    panel.spacing = unit(0.1, "cm"),
    panel.border = element_blank()
  ) +
  scale_x_discrete(breaks = c("01-04 yrs", "75+ yrs")) +
  guides(fill = guide_legend(ncol = 2))


f6_layout  <- c(area(t = 1, l = 1, r = 2),
                area(t = 2, l = 1),
                area(t = 2, l = 2),
                area(t = 3, l = 1, r = 2))


f6_p1 + free(f6_p2) + free(f6_p3) + f6_p4 + plot_layout(design = f6_layout) +
  plot_annotation(
    title = "State - Oaxaca",
    tag_levels = 'a',
    theme = theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.tag = element_text(size = 10)
    )
  )



ggsave(
  filename = "figure_6.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 190,
  height = 240,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  create.dir = FALSE
)
