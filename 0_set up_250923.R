setwd("/Users/joonkim/Thesis")

# Load the libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readxl)

DATA <- read_excel("WUENIC-2024-data.xlsx", sheet = "country_level")
DATA2 <- DATA
DATA2$category <- NA
DATA2$category[DATA$mics_2024=="Former-Gavi eligible countries"] <- "FGC"
DATA2$category[DATA$mics_2024 %in% c("Never-Gavi IDA-eligible economies", 
                                     "Never-Gavi eligible lower middle-income countries")] <- "NGC"
DATA2$category[DATA$mics_2024=="Gavi57"] <- "Gavi57"
DATA2$category2 <- DATA2$category
DATA2$category2[DATA2$category=="NGC" | DATA2$category=="FGC" ] <- "FGC+NGC"
DATA2$category3 <- DATA2$category
DATA2$category3[DATA2$category== "Gavi57" & DATA2$wb_long_2019 == "LMIC"] <- "Gavi_LMIC"
DATA2$category3[DATA2$category3== "Gavi57"]<- "Gavi_LIC"

save(DATA2, file="DATA2.RData")

#Former Gavi (19 countries)
FormerGavi <- subset(DATA, DATA$mics_2024=="Former-Gavi eligible countries") 
unique(FormerGavi$cname)

# Never Gavi (26 countries)
NeverGavi <- subset(DATA2, DATA2$mics_2024 %in% c("Never-Gavi eligible lower middle-income countries", "Never-Gavi IDA-eligible economies"))
unique(NeverGavi$cname)

