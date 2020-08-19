rm(list = ls())
library(tidyverse)

df_hydro <- readRDS("df_hydro")

mod <- glm(formula = score ~ . , data = df_hydro)
summary(mod)
