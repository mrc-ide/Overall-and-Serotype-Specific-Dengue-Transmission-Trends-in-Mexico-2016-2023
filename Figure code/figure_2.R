library(tidyverse)
library(MetBrewer)

mex_dengue_data_2016_2023 <- readRDS("mex_dengue_data_2016_2023.rds")

mex_dengue_serotypes_2016_2023_by_state_long <- mex_dengue_data_2016_2023  %>%
  filter(ESTATUS_CASO != 3,
         RESULTADO_PCR != 5,
         ENTIDAD_RES != "otros paises" & ENTIDAD_RES != "estados unidos de america" & ENTIDAD_RES != "otros paises de latinoamerica") %>%
  count(RESULTADO_PCR, YEAR, ENTIDAD_RES) %>%
  complete(RESULTADO_PCR, nesting(YEAR, ENTIDAD_RES), fill = list(n = 0)) %>%
  complete(YEAR, nesting(RESULTADO_PCR, ENTIDAD_RES), fill = list(n = 0)) %>%
  group_by(YEAR, ENTIDAD_RES) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(YEAR = as.factor(YEAR),
         ENTIDAD_RES_2 = ifelse(ENTIDAD_RES == "campeche", "Campeche",
                                ifelse(ENTIDAD_RES == "chiapas", "Chiapas",
                                       ifelse(ENTIDAD_RES == "coahuila", "Coahuila",
                                              ifelse(ENTIDAD_RES == "colima", "Colima",
                                                     ifelse(ENTIDAD_RES == "guerrero", "Guerrero",
                                                            ifelse(ENTIDAD_RES == "jalisco", "Jalisco",
                                                                   ifelse(ENTIDAD_RES == "michoacan", "Michoacán",
                                                                          ifelse(ENTIDAD_RES == "morelos", "Morelos",
                                                                                 ifelse(ENTIDAD_RES == "oaxaca", "Oaxaca",
                                                                                        ifelse(ENTIDAD_RES == "nuevo leon", "Nuevo León",
                                                                                               ifelse(ENTIDAD_RES == "puebla", "Puebla",
                                                                                                      ifelse(ENTIDAD_RES == "sinaloa", "Sinaloa",
                                                                                                             ifelse(ENTIDAD_RES == "tabasco", "Tabasco",
                                                                                                                    ifelse(ENTIDAD_RES == "veracruz", "Veracruz",
                                                                                                                           ifelse(ENTIDAD_RES == "aguascalientes", "Aguascalientes",
                                                                                                                                  ifelse(ENTIDAD_RES == "baja california", "Baja California",
                                                                                                                                         ifelse(ENTIDAD_RES == "baja california sur", "Baja California \n Sur",
                                                                                                                                                ifelse(ENTIDAD_RES == "durango", "Durango",
                                                                                                                                                       ifelse(ENTIDAD_RES == "guanajuato", "Guanajuato",
                                                                                                                                                              ifelse(ENTIDAD_RES == "hidalgo", "Hidalgo",
                                                                                                                                                                     ifelse(ENTIDAD_RES == "nayarit", "Nayarit",
                                                                                                                                                                            ifelse(ENTIDAD_RES == "sonora", "Sonora",
                                                                                                                                                                                   ifelse(ENTIDAD_RES == "mexico", "México",
                                                                                                                                                                                          ifelse(ENTIDAD_RES == "queretaro", "Querétaro",
                                                                                                                                                                                                 ifelse(ENTIDAD_RES == "san luis potosi", "San Luis Potosí",
                                                                                                                                                                                                        ifelse(ENTIDAD_RES == "distrito federal", "Ciudad de México",
                                                                                                                                                                                                               ifelse(ENTIDAD_RES == "chihuahua", "Chihuahua",
                                                                                                                                                                                                                      ifelse(ENTIDAD_RES == "quintana roo", "Quintana Roo",
                                                                                                                                                                                                                             ifelse(ENTIDAD_RES == "tlaxcala", "Tlaxcala",
                                                                                                                                                                                                                                    ifelse(ENTIDAD_RES == "yucatan", "Yucatán",
                                                                                                                                                                                                                                           ifelse(ENTIDAD_RES == "tamaulipas", "Tamaulipas",
                                                                                                                                                                                                                                                  ifelse(ENTIDAD_RES == "zacatecas", "Zacatecas", 
                                                                                                                                                                                                                                                         NA)))))))))))))))))))))))))))))))))



ggplot(aes(x = YEAR, y = prop, group = as.factor(RESULTADO_PCR), colour = as.factor(RESULTADO_PCR)), 
       data = mex_dengue_serotypes_2016_2023_by_state_long %>%
         filter(!is.na(prop)))+ 
  geom_line(linewidth = 0.7) +
  geom_point(size = 0.7) +
  facet_wrap( ~ ENTIDAD_RES_2, scales = "free_x", axes = "all", axis.labels = "all_x") +
  theme_classic() +
  theme(text = element_text(size=10),
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(0.1,"cm"),
        panel.border = element_blank(),
        legend.position = c(0.6, 0.01),
        legend.direction = "horizontal",
        axis.text.x=element_text(angle=90,hjust=1)) +
  ylim(0, 1) +
  labs(colour = "DENV Serotype",
       y = "Proportion",
       x = "Year") +
  scale_colour_manual(values = met.brewer("Egypt", 4)) 


ggsave(
  filename = "figure_2.tiff",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 190,
  height = 210,
  units = "mm",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  create.dir = FALSE
)