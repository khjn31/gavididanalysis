setwd("/Users/joonkim/Thesis")

# Load the libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readxl)

load("DATA2.RData")

Key_columns <- DATA2 %>% select(cname, vaccine, year, coverage, nvax, target, mics_2024, gni_2019, category, category2, category3)

Filtered_MCV1<- subset(Key_columns, is.na(Key_columns$category)==FALSE &
                     Key_columns$vaccine=="mcv1" &
                     Key_columns$year >= 2010 & Key_columns$year <= 2023)

# Function to calculate the weighted average by target population
weighted_avg <- function(coverage, target) {
  return(sum(coverage * target, na.rm = TRUE) / sum(target, na.rm = TRUE))
}

# Calculate weighted averages by category, year using target population
category_average_MCV1 <- Filtered_MCV1 %>%
  group_by(category, year) %>%
  mutate(weighted_coverage = weighted_avg (coverage, target))

# Calculate weighted averages by category2, year using target population
category_average2_MCV1 <- Filtered_MCV1 %>%
  group_by(category2, year) %>%
  mutate(weighted_coverage = weighted_avg(coverage, target))

# Calculate weighted averages by category3, year using target population
category_average3_MCV1 <- Filtered_MCV1 %>%
  group_by(category3, year) %>%
  mutate(weighted_coverage = weighted_avg (coverage, target))

# Combine unique entries by category, year, and coverage
unique_MCV1 <- category_average_MCV1 %>% 
  distinct(cname, category, year, vaccine = "mcv1", gni_2019, weighted_coverage)

unique2_MCV1 <- category_average2_MCV1 %>%
  distinct(cname, category2, year, vaccine = "mcv1", gni_2019, weighted_coverage)

unique3_MCV1 <- category_average3_MCV1 %>%
  distinct(cname, category3, year, vaccine = "mcv1", gni_2019, weighted_coverage)

# Filter to keep only FGC, NGC
unique3_filter_MCV1 <- subset(unique3_MCV1, category3 %in% c("FGC", "NGC"))

ggplot() +
  geom_point(data=unique3_filter_MCV1, aes(x=year, y=weighted_coverage, group = category3, color=category3))+
  geom_line(data=unique3_filter_MCV1, aes(x=year, y=weighted_coverage, group = category3, color=category3)) +
  scale_color_manual(values=c("NGC" = "blue", "FGC" = "red", "Gavi_LMIC" = "#43A047")) +  # Set colors explicitly
  scale_y_continuous(breaks=seq(0,100,10), limits=c(0,100)) +
  theme_bw()+
  labs(x= "Year", y="Coverage (%)", colour = "Group")+
  theme(title=element_text(angle=0, face= "bold", colour= "black", size=11),
        axis.title.y=element_text(angle=90, face= "bold", colour= "black", size=11),
        axis.title.x=element_text(angle=0, face= "bold", colour= "black", size=11),
        axis.text.x = element_text(size=11, face= "bold", colour= "black" ),
        axis.text.y = element_text(size=11, face= "bold", colour= "black" ))

###############with numeric points########################
ggplot(data = unique3_filter_MCV1, 
       aes(x = year, y = weighted_coverage, group = category3, color = category3)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1) +
  
  # ➕ Add numeric labels above each point
  geom_text(aes(label = round(weighted_coverage, 1)),
            vjust = -0.8, size = 3.5, fontface = "bold", show.legend = FALSE) +
  
  scale_color_manual(values = c("NGC" = "blue", "FGC" = "red", "Gavi_LMIC" = "#43A047")) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
  theme_bw() +
  labs(x = "Year", y = "Coverage (%)", colour = "Group") +
  theme(
    title = element_text(angle = 0, face = "bold", colour = "black", size = 12),
    axis.title.y = element_text(angle = 90, face = "bold", colour = "black", size = 12),
    axis.title.x = element_text(angle = 0, face = "bold", colour = "black", size = 12),
    axis.text.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 12, face = "bold", colour = "black")
  )

# Filter to exclude 'Gavi_LIC' and keep only FGC, NGC, and Gavi_LMIC
unique3_filter_MCV1_Gavi_LMIC <- subset(unique3_MCV1, category3 %in% c("FGC", "NGC", "Gavi_LMIC"))

ggplot() +
  geom_point(data=unique3_filter_MCV1_Gavi_LMIC, aes(x=year, y=weighted_coverage, group = category3, color=category3))+
  geom_line(data=unique3_filter_MCV1_Gavi_LMIC, aes(x=year, y=weighted_coverage, group = category3, color=category3)) +
  scale_color_manual(values=c("NGC" = "blue", "FGC" = "red", "Gavi_LMIC" = "#43A047")) +  # Set colors explicitly
  scale_y_continuous(breaks=seq(65,95,5), limits=c(65,95)) +
  theme_bw()+
  labs(x= "Year", y="Coverage", colour = "Group")+
  theme(title=element_text(angle=0, face= "bold", colour= "black", size=12),
        axis.title.y=element_text(angle=90, face= "bold", colour= "black", size=12),
        axis.title.x=element_text(angle=0, face= "bold", colour= "black", size=12),
        axis.text.x = element_text(size=12, face= "bold", colour= "black" ),
        axis.text.y = element_text(size=12, face= "bold", colour= "black" ))

#########################################################################################
#######.MCV1## FGC vs NGC (Pre:2017-19 vs Post:20-21)_short term_PTA test_DID analysis  ########## 
#########################################################################################

# Filter country-level data for relevant years and categories
data_FGC_NGC_st_mcv1 <- unique_MCV1 %>%
  filter(category %in% c("FGC","NGC"), year >= 2017, year <= 2021) %>%
  mutate(
    cname          = as.factor(cname),
    treatment      = ifelse(category == "FGC", 1L, 0L),
    post_treatment = ifelse(year >= 2020, 1L, 0L),
    gni_scaled     = (gni_2019 - 4000) / 1000   # scale GNI
  )

#Test parallel trend assumption
# Filter data for the pre-treatment period (2015-2019) and relevant categories (FGC, NGC)
pre_period_FGC_NGC_mcv1 <- data_FGC_NGC_st_mcv1 %>%
  filter(year >= 2017 & year <= 2019 & (category == "FGC" | category == "NGC"))

# Create a treatment indicator (1 for FGC, 0 for NGC)
pre_period_FGC_NGC_mcv1 <- pre_period_FGC_NGC_mcv1 %>%
  mutate(treatment = ifelse(category == "FGC", 1, 0))

# Run the regression with an interaction between treatment and year to test the parallel trend assumption
parallel_trend_model_mcv1 <- lm(weighted_coverage ~ year * treatment + gni_scaled, data = pre_period_FGC_NGC_mcv1)

# Display summary of the parallel trend model
summary(parallel_trend_model_mcv1)

#DiD analysis#
# Run the DiD model with country and year fixed effects, controlling for gni_scaled
did_model_FGC_NGC_st_mcv1 <- lm(weighted_coverage ~ treatment * post_treatment + gni_scaled, data = data_FGC_NGC_st_mcv1)

# View the summary of the model
summary(did_model_FGC_NGC_st_mcv1)

###############################################################################
####### MCV1: FGC vs NGC (Pre:2017-19 vs Post:20-23)_long term_DID analysis w/o FEs ########## 
###############################################################################

# Filter country-level data for relevant years and categories
data_FGC_NGC_lt_mcv1 <- unique_MCV1 %>%
  filter(category %in% c("FGC","NGC"), year >= 2017, year <= 2023) %>%
  mutate(
    cname          = as.factor(cname),
    treatment      = ifelse(category == "FGC", 1L, 0L),
    post_treatment = ifelse(year >= 2020, 1L, 0L),
    gni_scaled     = (gni_2019 - 4000) / 1000   # scale GNI
  )

#DiD analysis#
# Run the DiD model without country and year fixed effects, controlling for gni_scaled
did_model_FGC_NGC_lt_mcv1 <- lm(weighted_coverage ~ treatment * post_treatment + gni_scaled, data = data_FGC_NGC_lt_mcv1)

# View the summary of the model
summary(did_model_FGC_NGC_lt_mcv1)

####TABLE GENERATION#######
# Generate the table for both short-term and long-term models with confidence intervals
stargazer(did_model_FGC_NGC_st_mcv1,  # Short-term model
          did_model_FGC_NGC_lt_mcv1,  # Intermediate-term model
          type = "text",         # Use "html" for HTML or "latex" output
          title = "Difference-in-Differences Analysis for MCV1: FGC vs NGC (Short Term and Intermediate Term)", 
          align = TRUE, 
          dep.var.labels = c("Weighted Coverage"),  # Shared dependent variable
          column.labels = c("(1) Short Term", "(2) Long Term"),  # Labels for short and long term
          covariate.labels = c("Treatment (FGC)", 
                               "Post Treatment", 
                               "Treatment × Post Treatment", 
                               "Constant"), 
          omit.stat = c("f", "ser"),  # Omits F-stat and Standard Error
          ci = TRUE,                 # Add confidence intervals
          ci.level = 0.95,           # Set the confidence interval level to 95%
          digits = 3)

########################MCV1 w/t FEs###################################
data_FGC_NGC_st_MCV1 <- Filtered_MCV1 %>%
  dplyr::filter(category %in% c("FGC","NGC"),
                year >= 2017, year <= 2021) %>%
  dplyr::mutate(
    cname          = as.factor(cname),
    treatment      = ifelse(category == "FGC", 1L, 0L),
    post_treatment = ifelse(year >= 2020, 1L, 0L),
    gni_scaled     = (gni_2019 - 4000) / 1000
  )

# FE-DiD: Weighted (child-population average)
did_fe_st_w_MCV1 <- fixest::feols(
  coverage ~ treatment:post_treatment | cname + year,
  data     = data_FGC_NGC_st_MCV1,
  weights  = ~ target,
  cluster  = ~ cname
)

# FE-DiD: Unweighted (country-average)
did_fe_st_u_MCV1 <- fixest::feols(
  coverage ~ treatment:post_treatment | cname + year,
  data     = data_FGC_NGC_st_MCV1,
  cluster  = ~ cname
)

# -------------------------------
# 2. Long-term (LT) panel: 2017–2023
# -------------------------------
data_FGC_NGC_lt_MCV1 <- Filtered_MCV1 %>%
  dplyr::filter(category %in% c("FGC","NGC"),
                year >= 2017, year <= 2023) %>%
  dplyr::mutate(
    cname          = as.factor(cname),
    treatment      = ifelse(category == "FGC", 1L, 0L),
    post_treatment = ifelse(year >= 2020, 1L, 0L),
    gni_scaled     = (gni_2019 - 4000) / 1000
  )

# FE-DiD: Weighted
did_fe_lt_w_MCV1 <- fixest::feols(
  coverage ~ treatment:post_treatment | cname + year,
  data     = data_FGC_NGC_lt_MCV1,
  weights  = ~ target,
  cluster  = ~ cname
)

# FE-DiD: Unweighted
did_fe_lt_u_MCV1 <- fixest::feols(
  coverage ~ treatment:post_treatment | cname + year,
  data     = data_FGC_NGC_lt_MCV1,
  cluster  = ~ cname
)

# -------------------------------
# 3. Table presentation
# -------------------------------
etable(
  list(
    "ST Weighted"   = did_fe_st_w_MCV1,
    "ST Unweighted" = did_fe_st_u_MCV1,
    "LT Weighted"   = did_fe_lt_w_MCV1,
    "LT Unweighted" = did_fe_lt_u_MCV1
  ),
  se       = "cluster",
  coefstat = "confint",  # show 95% CIs
  dict     = c("treatment:post_treatment" = "Treatment × Post (DiD)"),
  fitstat  = c("n", "r2", "wr2"),
  digits   = 2,
  title    = "MCV1: Difference-in-Differences Estimates (2017–2021 vs 2017–2023)\nWeighted vs Unweighted, Country & Year Fixed Effects"
)

