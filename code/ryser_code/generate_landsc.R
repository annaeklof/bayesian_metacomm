
require(tidyverse)
source("functions.R")

no_patches <- seq(16, 70, 2)
for (n in no_patches) {
  lsc <- generate_coordinates(n, 2) %>% as_tibble %>% rename(x=V1, y=V2)
  write_csv(lsc, paste0("../../data/ryser_data/landscapes/lsc_", n, ".csv"))
}

