# Course population biology
# Acacia suaveolens ####

# Markus Bauer
# 2023-10-12

### Packages ###
library(here)
library(tidyverse)
library(renv)
library(ggbeeswarm)

### Start ###
rm(list = ls())
renv::status()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Load data ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#Silvertown & Charlesworth (2001) Tab. 1.1 ISBN: 978-0-632-04991-2

data <- read_csv(
  here("data", "data_raw_acacia_suaveolens.csv"),
  col_names = TRUE,
  na = c("", "NA", "na"),
  col_types =
    cols(
      .default = "?"
    )
)

### * Functions ####
theme_adjusted <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text = element_text(angle = 0, hjust = 0.5, size = 9,
                             color = "black"),
    axis.title = element_text(angle = 0, hjust = 0.5, size = 9,
                              color = "black"),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 9),
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Analyses ###################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Data exploration ##########################################################


### a Graphs of raw data -------------------------------------------------------

data %>%
  slice(-1) %>%
  ggplot(aes(x = age, y = number)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ log(x)) +
  theme_adjusted() +
  labs(y = "Number of Individuals [#]", x = "Age in years")

data %>%
  slice(-1) %>%
  ggplot(aes(x = age, y = seeds)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ log(x)) +
  theme_adjusted() +
  labs(y = "Seeds per individual [#]", x = "Age in years")
