library(tidyverse)
library(MetBrewer)

FOI_m4 <- readRDS("FOI_m4.rds")
lam_df_m2 <- readRDS("lam_df_m2.rds")

FOI_df_m4 <- bind_rows(FOI_m4)

foi_comparison_m2_m4 <- bind_rows(
  lam_df_m2 %>%
    mutate(
      med = 4 * med,
      ciL = 4 * ciL,
      ciU = 4 * ciU
    ) %>%
    filter(
      ENTIDAD_RES %in% c(
        "chiapas",
        "guerrero",
        "jalisco",
        "michoacan",
        "morelos",
        "nuevo leon",
        "oaxaca",
        "tabasco",
        "veracruz",
        "coahuila",
        "colima",
        "sinaloa",
        "campeche",
        "puebla",
        "mexico",
        "guanajuato"
      )
    ),
  FOI_df_m4 %>% mutate(type = paste0("lam_", type)) %>%
    filter(
      ENTIDAD_RES %in% c(
        "chiapas",
        "guerrero",
        "jalisco",
        "michoacan",
        "morelos",
        "nuevo leon",
        "oaxaca",
        "tabasco",
        "veracruz",
        "coahuila",
        "colima",
        "sinaloa",
        "campeche",
        "puebla",
        "mexico",
        "guanajuato"
      )
    ),
  .id = "id"
) %>%
  mutate(id = ifelse(
    id == 1,
    "Non-serotype-specific (Model 2)",
    "Serotype-specific (Model 4)"
  )) %>%
  mutate(year = ifelse(type == "lam_T1", 2016, ifelse(
    type == "lam_T2", 2017, ifelse(type == "lam_T3", 2018, ifelse(
      type == "lam_T4", 2019, ifelse(type == "lam_T5", "2020", ifelse(
        type == "lam_T6", 2021, ifelse(type == "lam_T7", 2022, 2023)
      ))
    ))
  ))) %>%
  mutate(ENTIDAD_RES = ifelse(
    ENTIDAD_RES == "campeche",
    "Campeche",
    ifelse(
      ENTIDAD_RES == "chiapas",
      "Chiapas",
      ifelse(
        ENTIDAD_RES == "coahuila",
        "Coahuila",
        ifelse(
          ENTIDAD_RES == "colima",
          "Colima",
          ifelse(
            ENTIDAD_RES == "guanajuato",
            "Guanajuato",
            ifelse(
              ENTIDAD_RES == "guerrero",
              "Guerrero",
              ifelse(
                ENTIDAD_RES == "jalisco",
                "Jalisco",
                ifelse(
                  ENTIDAD_RES == "mexico",
                  "México",
                  ifelse(
                    ENTIDAD_RES == "michoacan",
                    "Michoacán",
                    ifelse(
                      ENTIDAD_RES == "morelos",
                      "Morelos",
                      ifelse(
                        ENTIDAD_RES == "oaxaca",
                        "Oaxaca",
                        ifelse(
                          ENTIDAD_RES == "nuevo leon",
                          "Nuevo Léon",
                          ifelse(
                            ENTIDAD_RES == "puebla",
                            "Puebla",
                            ifelse(
                              ENTIDAD_RES == "sinaloa",
                              "Sinaloa",
                              ifelse(
                                ENTIDAD_RES == "tabasco",
                                "Tabasco",
                                ifelse(ENTIDAD_RES == "veracruz", "Veracruz", NA)
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  ))


ggplot(foi_comparison_m2_m4) +
  geom_point(aes(x = year, y = med, colour = id), position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(
      year,
      ymin = ciL,
      ymax = ciU,
      colour = id
    ),
    width = 0,
    position = position_dodge(width = 0.5)
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 10),
    legend.position = "bottom",
    strip.background = element_rect(fill = NA, color = "white"),
    panel.spacing = unit(0.1, "cm"),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  labs(y = "Total FOI", x = "Year") +
  scale_colour_manual(values = met.brewer("Derain", 2), name = "Model", ) +
  facet_wrap( ~ ENTIDAD_RES, scales = "free")


ggsave(
  filename = "figure_7.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 190,
  height = 190,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  create.dir = FALSE
)