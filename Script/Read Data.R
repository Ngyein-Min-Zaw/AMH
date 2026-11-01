url <- "https://github.com/Ngyein-Min-Zaw/AMH/raw/refs/heads/main/Data/UV%20(04.08.2006)%20-%20completed.xlsx"
download.file(url,
     destfile = "./Data/downloaded_file.xlsx",
     mode = "wb"
)

library(readxl)
UV <- read_excel("./Data/downloaded_file.xlsx")

library(rsm)
library(tidyverse)
library(extrafont)
library(readxl)
library(ggtext)

## Load Data
RSM <- read_excel("./Data/(26.06.2026) absorbance data.xlsx")
RSM2 <- RSM %>% mutate(Temperature = str_remove(Temp, "°C") %>% as.numeric())

mfont <- "Times New Roman"
