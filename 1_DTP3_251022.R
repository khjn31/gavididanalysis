setwd("/Users/joonkim/Thesis")

# Load the libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readxl)

load("DATA2.RData")

Key_columns <- DATA2 %>% select(cname, vaccine, year, coverage, nvax, target, mics_2024, gni_2019, category,category2, category3)

Filtered_DTP3 <- subset(Key_columns, is.na(Key_columns$category)==FALSE &
                          Key_columns$vaccine=="dtp3" &
                          Key_columns$year >= 2010 & Key_columns$year <= 2023)

# Function to calculate the weighted average by target population
weighted_avg <- function(coverage, target) {
 return(sum(coverage * target, na.rm = TRUE) / sum(target, na.rm = TRUE))
}

# Calculate weighted averages by category, year using target population
category_average_DTP3 <- Filtered_DTP3 %>%
 group_by(category, year) %>%
 mutate(weighted_coverage = weighted_avg (coverage, target))

# Calculate weighted averages by category2, year using target population
category_average2_DTP3 <- Filtered_DTP3 %>%
 group_by(category2, year) %>%
 mutate(weighted_coverage = weighted_avg(coverage, target))

# Calculate weighted averages by category3, year using target population
category_average3_DTP3 <- Filtered_DTP3 %>%
 group_by(category3, year) %>%
 mutate(weighted_coverage = weighted_avg (coverage, target))

# Combine unique entries by category, year, and coverage
unique_DTP3 <- category_average_DTP3 %>%
 distinct(cname, category, year, vaccine = "dtp3", gni_2019, weighted_coverage)
unique2_DTP3 <- category_average2_DTP3 %>%
 distinct(cname, category2, year, vaccine = "dtp3", gni_2019, weighted_coverage)
unique3_DTP3 <- category_average3_DTP3 %>%
 distinct(cname, category3, year, vaccine = "dtp3", gni_2019, weighted_coverage)

# Filter to exclude 'Gavi_LIC' and keep only FGC and NGC
unique_filter_DTP3 <- subset(unique_DTP3, category %in% c("FGC", "NGC"))

ggplot() +
  geom_point(data=unique_filter_DTP3, aes(x=year, y=weighted_coverage, group = category, color=category))+
  geom_line(data=unique_filter_DTP3, aes(x=year, y=weighted_coverage, group = category, color=category)) +
  scale_color_manual(values=c("NGC" = "blue", "FGC" = "red", "Gavi_LMIC" = "#43A047")) +  # Set colors explicitly
  scale_y_continuous(breaks=seq(0,100,10), limits=c(0,100)) +
  scale_x_continuous(breaks=seq(min(unique_filter_DTP3$year),
                                max(unique_filter_DTP3$year), 
                                by=2)) +   # tick marks every 2 years
  theme_bw()+
  labs(x= "Year", y="Coverage (%)", colour = "Group")+
  theme(title=element_text(angle=0, face= "bold", colour= "black", size=11),
        axis.title.y=element_text(angle=90, face= "bold", colour= "black", size=11),
        axis.title.x=element_text(angle=0, face= "bold", colour= "black", size=11),
        axis.text.x = element_text(size=11, face= "bold", colour= "black" ),
        axis.text.y = element_text(size=11, face= "bold", colour= "black" ))

##############same graph with numeric points######################
ggplot() +
  geom_point(data = unique_filter_DTP3,
             aes(x = year, y = weighted_coverage, group = category, color = category),
             size = 2.5) +
  geom_line(data = unique_filter_DTP3,
            aes(x = year, y = weighted_coverage, group = category, color = category),
            linewidth = 1) +
  geom_text(data = unique_filter_DTP3,
            aes(x = year, y = weighted_coverage, label = round(weighted_coverage, 1), color = category),
            vjust = -0.8, size = 3.5, fontface = "bold", show.legend = FALSE) +  # add value labels
  
  scale_color_manual(values = c("NGC" = "blue", "FGC" = "red", "Gavi_LMIC" = "#43A047")) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(min(unique_filter_DTP3$year),
                                  max(unique_filter_DTP3$year),
                                  by = 2)) +
  theme_bw() +
  labs(x = "Year", y = "Coverage (%)", colour = "Group") +
  theme(
    title = element_text(angle = 0, face = "bold", colour = "black", size = 12),
    axis.title.y = element_text(angle = 90, face = "bold", colour = "black", size = 12),
    axis.title.x = element_text(angle = 0, face = "bold", colour = "black", size = 12),
    axis.text.x = element_text(size = 10, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 10, face = "bold", colour = "black")
  )

# Filter country-level data for relevant years and categories
data_FGC_NGC_st <- unique_DTP3 %>%
  filter(category %in% c("FGC","NGC"), year >= 2017, year <= 2021) %>%
  mutate(
    cname          = as.factor(cname),
    treatment      = ifelse(category == "FGC", 1L, 0L),
    post_treatment = ifelse(year >= 2020, 1L, 0L),
    gni_scaled     = (gni_2019 - 4000) / 1000   # scale GNI
  )

#Test parallel trend assumption
# Filter data for the pre-treatment period (2015-2019) and relevant categories (FGC, NGC)
pre_period_FGC_NGC <- data_FGC_NGC_st %>%
  filter(year >= 2017 & year <= 2019 & (category == "FGC" | category == "NGC"))

# Create a treatment indicator (1 for FGC, 0 for NGC)
pre_period_FGC_NGC <- pre_period_FGC_NGC %>%
  mutate(treatment = ifelse(category == "FGC", 1, 0))

# Run the regression with an interaction between treatment and year to test the parallel trend assumption
parallel_trend_model <- lm(coverage ~ year * treatment + gni_scaled, data = pre_period_FGC_NGC)

# Display summary of the parallel trend model
summary(parallel_trend_model)

############################################################################################
##################### OLS DID MODEL WITHOUT FEs#############################################
############################################################################################

# Run the DiD model with country and year fixed effects, controlling for gni_2019
did_model_FGC_NGC_st <- lm(weighted_coverage ~ treatment * post_treatment + gni_scaled, data = data_FGC_NGC_st)
# View the summary of the model
summary(did_model_FGC_NGC_st)

# Run the DiD model with country without controlling for gni_2019
did_model_FGC_NGC_st_no_gni <- lm(weighted_coverage ~ treatment * post_treatment, data = data_FGC_NGC_st)
# View the summary of the model
summary(did_model_FGC_NGC_st_no_gni)

# Install and load the stargazer package if not already installed
if (!requireNamespace("stargazer", quietly = TRUE)) {
  install.packages("stargazer")
}
library(stargazer)

# Generate the table for the model with 95% CI
stargazer(did_model_FGC_NGC_st, 
          type = "text",         # Use "html" for HTML or "latex" for LaTeX output
          title = "Difference-in-Differences Analysis for DTP3: FGC vs NGC (Short Term)", 
          align = TRUE, 
          dep.var.labels = c("Weighted Coverage"), 
          covariate.labels = c("Treatment (FGC)", 
                               "Post Treatment (2020-21)",
                               "Treatment × Post Treatment"), 
          omit.stat = c("f", "ser"),  # Omits F-stat and Standard Error
          ci = TRUE,                 # Add confidence intervals
          ci.level = 0.95,           # Set the confidence interval level to 95%
          digits = 3)


# Generate the table for the model
stargazer(did_model_FGC_NGC_st_no_gni, 
          type = "text",         # Use "html" for HTML or "latex" for LaTeX output
          title = "Difference-in-Differences Analysis for DTP3: FGC vs NGC (Short Term)", 
          align = TRUE, 
          dep.var.labels = c("Weighted Coverage"), 
          covariate.labels = c("Treatment (FGC)", 
                               "Post Treatment (2020-21)",
                               "Treatment × Post Treatment"), 
          omit.stat = c("f", "ser"),  # Omits F-stat and Standard Error
          ci = TRUE,                 # Add confidence intervals
          ci.level = 0.95,           # Set the confidence interval level to 95%
          digits = 3)

###############################################################################
####### FGC vs NGC (Pre:2017-19 vs Post:20-23)_long term_OLS DID analysis ########## 
###############################################################################

# Filter country-level data for relevant years and categories
data_FGC_NGC_lt <- unique_DTP3 %>%
  filter(category %in% c("FGC","NGC"), year >= 2017, year <= 2023) %>%
  mutate(
    cname          = as.factor(cname),
    treatment      = ifelse(category == "FGC", 1L, 0L),
    post_treatment = ifelse(year >= 2020, 1L, 0L),
    gni_scaled     = (gni_2019 - 4000) / 1000   # scale GNI
  )

#DiD analysis#
# Run the DiD model with country and year fixed effects, controlling for gni_2019
did_model_FGC_NGC_lt <- lm(weighted_coverage ~ treatment * post_treatment + gni_scaled, data = data_FGC_NGC_lt)

# View the summary of the model
summary(did_model_FGC_NGC_lt)

#DiD analysis#
# Run the DiD model without gni_2019 controlled
did_model_FGC_NGC_lt_no_gni <- lm(weighted_coverage ~ treatment * post_treatment, data = data_FGC_NGC_lt)

# View the summary of the model
summary(did_model_FGC_NGC_lt_no_gni)

# Generate the table for the model
stargazer(did_model_FGC_NGC_lt_no_gni, 
          type = "text",         # Use "html" for HTML or "latex" for LaTeX output
          title = "Difference-in-Differences Analysis for DTP3: FGC vs NGC (Long Term)", 
          align = TRUE, 
          dep.var.labels = c("Weighted Coverage"), 
          covariate.labels = c("Treatment (FGC)", 
                               "Post Treatment (2020-23)", 
                               "Treatment × Post Treatment"), 
          omit.stat = c("f", "ser"),  # Omits F-stat and Standard Error
          ci = TRUE,                 # Add confidence intervals
          ci.level = 0.95,           # Set the confidence interval level to 95%
          digits = 3)

# Generate the table for both short-term and long-term models
stargazer(did_model_FGC_NGC_st,  # Short-term model
          did_model_FGC_NGC_lt,  # Long-term model
          type = "text",         # Use "html" for HTML or "latex" output
          title = "Difference-in-Differences Analysis for DTP3: FGC vs NGC (Short and Long Term)", 
          align = TRUE, 
          dep.var.labels = c("Weighted Coverage"),  # Shared dependent variable
          column.labels = c("(1)", "(2)"),          # Labels for short and long term
          covariate.labels = c("Treatment (FGC)", 
                               "Post Treatment", 
                               "Treatment × Post Treatment", 
                               "TConstant"), 
          omit.stat = c("f", "ser"),  # Omits F-stat and Standard Error
          ci = TRUE,                 # Add confidence intervals
          ci.level = 0.95,           # Set the confidence interval level to 95%
          digits = 3)

# Generate the table for both short-term and long-term models with confidence intervals
stargazer(did_model_FGC_NGC_st,  # Short-term model
          did_model_FGC_NGC_lt,  # Long-term model
          type = "text",         # Use "html" for HTML or "latex" output
          title = "Difference-in-Differences Analysis for DTP3: FGC vs NGC (Short Term and Long Term)", 
          align = TRUE, 
          dep.var.labels = c("Weighted Coverage"),  # Shared dependent variable
          column.labels = c("(1) Short Term", "(2) Long Term"),  # Labels for short and long term
          covariate.labels = c("Treatment (FGC)", 
                               "Post Treatment", 
                               "GNI per Capita (2019)", 
                               "Treatment × Post Treatment", 
                               "Constant"), 
          ci = TRUE,            # Adds confidence intervals
          ci.level = 0.95,      # 95% confidence interval
          omit.stat = c("f", "ser"),  # Omits F-stat and Standard Error
          digits = 3)

####################### ROBUSTNESS CHECKS & SENSITIVTY ANALAYSIS ###################################
############ DTP3 without Indonesia ##########################
# Filter the dataset to exclude Indonesia
Filtered_no_indonesia <- Key_columns %>%
  filter(!is.na(category), cname != "Indonesia",  # Exclude Indonesia
         vaccine == "dtp3", year >= 2010 & year <= 2023)  # Filter relevant criteria

# Function to calculate the weighted average by target population
weighted_avg <- function(coverage, target) {
  return(sum(coverage * target, na.rm = TRUE) / sum(target, na.rm = TRUE))
}

# Calculate weighted averages by category (without Indonesia)
category_average_no_indonesia <- Filtered_no_indonesia %>%
  group_by(category, year) %>%
  summarize(
    weighted_coverage = weighted_avg(coverage, target),  # Weighted average
    .groups = "drop"  # Ensure grouping is dropped after summarization
  )

# Add the weighted averages back to the main dataset (without Indonesia)
Filtered_no_indonesia <- Filtered_no_indonesia %>%
  left_join(category_average_no_indonesia, by = c("category", "year"))

# Intermediate-term DiD analysis without Indonesia
data_FGC_NGC_no_indonesia_st <- Filtered_no_indonesia %>%
  filter(category %in% c("FGC","NGC"), year >= 2017, year <= 2023) %>%
  mutate(
    cname          = as.factor(cname),
    treatment      = ifelse(category == "FGC", 1L, 0L),
    post_treatment = ifelse(year >= 2020, 1L, 0L),
  )

# Run the DiD model
did_model_FGC_NGC_no_indonesia_lt <- lm(
  weighted_coverage ~ treatment * post_treatment + gni_2019,
  data = data_FGC_NGC_no_indonesia_lt
)

# Display the model summary
summary(did_model_FGC_NGC_no_indonesia_lt)

library(stargazer)
stargazer(
  did_model_FGC_NGC_no_indonesia_lt,
  type = "text",
  title = "Difference-in-Differences Analysis for DTP3: FGC vs NGC (Excluding Indonesia)",
  dep.var.labels = c("Weighted Coverage"),
  covariate.labels = c("Treatment (FGC)", "Post Treatment", "GNI per Capita (2019)", "Treatment × Post Treatment"),
  omit.stat = c("f", "ser"),
  ci = TRUE,                 # Add confidence intervals
  ci.level = 0.95,           # Set the confidence interval level to 95%
  digits = 3
)

# Graph of FGC and NGC without Indonesia
unique_filter_no_indonesia <- subset(Filtered_no_indonesia, category %in% c("FGC", "NGC"))

ggplot() +
  geom_point(data=unique_filter_no_indonesia, aes(x=year, y=weighted_coverage, group = category, color=category))+
  geom_line(data=unique_filter_no_indonesia, aes(x=year, y=weighted_coverage, group = category, color=category)) +
  scale_y_continuous(breaks=seq(0,100,10), limits=c(0,100)) +
  theme_bw()+
  labs(x= "Year", y="Coverage", colour = "Group")+
  theme(title=element_text(angle=0, face= "bold", colour= "black", size=12),
        axis.title.y=element_text(angle=90, face= "bold", colour= "black", size=12),
        axis.title.x=element_text(angle=0, face= "bold", colour= "black", size=12),
        axis.text.x = element_text(size=12, face= "bold", colour= "black" ),
        axis.text.y = element_text(size=12, face= "bold", colour= "black" ))

################################################################
############# Sensitivity check without Vietnam ##############
##############################################################

# Filter the dataset to exclude Vietnam
Filtered_no_vietnam <- Key_columns %>%
  filter(!is.na(category), cname != "Viet Nam",  # Exclude Vietnam
         vaccine == "dtp3", year >= 2010 & year <= 2023)  # Filter relevant criteria

# Function to calculate the weighted average by target population
weighted_avg <- function(coverage, target) {
  return(sum(coverage * target, na.rm = TRUE) / sum(target, na.rm = TRUE))
}

# Calculate weighted averages by category (without Vietnam)
category_average_no_vietnam <- Filtered_no_vietnam %>%
  group_by(category, year) %>%
  summarize(
    weighted_coverage = weighted_avg(coverage, target),  # Weighted average
    .groups = "drop"  # Ensure grouping is dropped after summarization
  )

# Add the weighted averages back to the main dataset (without Vietnam)
Filtered_no_vietnam <- Filtered_no_vietnam %>%
  left_join(category_average_no_vietnam, by = c("category", "year"))

# Long-term analysis (2017-2023)
data_FGC_NGC_no_vietnam_lt <- Filtered_no_vietnam %>%
  filter(category %in% c("FGC", "NGC"), year >= 2017 & year <= 2023) %>%  # Filter for relevant years
  mutate(
    treatment = ifelse(category == "FGC", 1, 0),  # Treatment for FGC countries
    post_treatment = ifelse(year >= 2020, 1, 0)   # Post-2020 for post-treatment period
  )

# Run the DiD model
did_model_FGC_NGC_no_vietnam_lt <- lm(
  weighted_coverage ~ treatment * post_treatment + gni_2019,
  data = data_FGC_NGC_no_vietnam_lt
)

# Display the model summary
summary(did_model_FGC_NGC_no_vietnam_lt)

# Intermediate-term sensitivity analysis table
stargazer(
  did_model_FGC_NGC_no_vietnam_lt,
  type = "text",
  title = "Difference-in-Differences Analysis for DTP3: FGC vs NGC (Excluding Vietnam)",
  dep.var.labels = c("Weighted Coverage"),
  covariate.labels = c("Treatment (FGC)", "Post Treatment", "GNI per Capita (2019)", "Treatment × Post Treatment"),
  omit.stat = c("f", "ser"),
  ci = TRUE,                 # Add confidence intervals
  ci.level = 0.95,           # Set the confidence interval level to 95%
  digits = 3
)

##### GRAPH without Vietnam ########
# Filter to exclude Vietnam and keep only FGC and NGC
unique_filter_no_vietnam <- subset(Filtered_no_vietnam, category %in% c("FGC", "NGC"))

ggplot() +
  geom_point(data = unique_filter_no_vietnam, aes(x = year, y = weighted_coverage, group = category, color = category)) +
  geom_line(data = unique_filter_no_vietnam, aes(x = year, y = weighted_coverage, group = category, color = category)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
  theme_bw() +
  labs(x = "Year", y = "Coverage", colour = "Group") +
  theme(
    title = element_text(angle = 0, face = "bold", colour = "black", size = 12),
    axis.title.y = element_text(angle = 90, face = "bold", colour = "black", size = 12),
    axis.title.x = element_text(angle = 0, face = "bold", colour = "black", size = 12),
    axis.text.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 12, face = "bold", colour = "black")
  )


################ Sensitivity test without 2023 ##########################

# Filter country-level data for relevant years and categories
data_FGC_NGC_st_without2023 <- unique_DTP3 %>%
  filter(category %in% c("FGC","NGC"), year >= 2017, year <= 2022) %>%
  mutate(
    cname          = as.factor(cname),
    treatment      = ifelse(category == "FGC", 1L, 0L),
    post_treatment = ifelse(year >= 2020, 1L, 0L),
    gni_scaled     = (gni_2019 - 4000) / 1000   # scale GNI
  )

# Run the DiD model without 2023
did_model_FGC_NGC_st_without2023 <- lm(weighted_coverage ~ treatment * post_treatment + gni_scaled, data = data_FGC_NGC_st_without2023)
# View the summary of the model
summary(did_model_FGC_NGC_st_without2023)

# Intermediate-term sensitivity analysis table
stargazer(
  did_model_FGC_NGC_st_without2023,
  type = "text",
  title = "Difference-in-Differences Analysis for DTP3: FGC vs NGC (Excluding 2023)",
  dep.var.labels = c("Weighted Coverage"),
  covariate.labels = c("Treatment (FGC)", "Post Treatment", "GNI per Capita (2019)", "Treatment × Post Treatment"),
  omit.stat = c("f", "ser"),
  ci = TRUE,                 # Add confidence intervals
  ci.level = 0.95,           # Set the confidence interval level to 95%
  digits = 3
)

###################################################################
##Generate descriptive statistics for FGC and NGC groups (2010-2023)
####################################################################
descriptive_table_DTP3 <- Filtered_DTP3 %>%
  filter(category %in% c("FGC","NGC"), year >= 2010, year <= 2023) %>%
  group_by(category) %>%
  summarise(
    Mean = sum(coverage * target, na.rm = TRUE) / sum(target, na.rm = TRUE), # weighted mean
    SD   = sd(coverage, na.rm = TRUE),     # unweighted SD
    Observations = n(),                    # number of country-year observations
    .groups = "drop"
  )

descriptive_table_DTP3

# Convert to data frame
descriptive_table_df <- as.data.frame(descriptive_table_DTP3)

# Stargazer output
stargazer(
  descriptive_table_df,
  type = "text",
  summary = FALSE,
  title = "Weighted coverage for FGC and NGC (2010–2023)",
  rownames = FALSE,
  digits = 2
)

####### combined summary table ##########

# Define a function to calculate descriptive statistics for a given dataset and vaccine
calculate_summary <- function(data, vaccine_name) {
  data %>%
    filter(category %in% c("FGC", "NGC"), year >= 2010 & year <= 2023) %>%  # Filter by category and year
    group_by(category) %>%  # Group by category
    summarize(
      Vaccine = vaccine_name,                              # Add the vaccine name as a column
      Mean = mean(weighted_coverage, na.rm = TRUE),        # Calculate mean of weighted coverage
      SD = sd(weighted_coverage, na.rm = TRUE),            # Calculate standard deviation
      Observations = n()                                   # Count the number of observations
    )
}

# Calculate descriptive statistics for each vaccine
summary_DTP3 <- calculate_summary(unique_DTP3, "DTP3")
summary_DTP1 <- calculate_summary(unique_DTP1, "DTP1")
summary_MCV1 <- calculate_summary(unique_MCV1, "MCV1")

# Combine the results into a single table
combined_summary <- bind_rows(summary_DTP3, summary_DTP1, summary_MCV1)

# Convert to a data frame for stargazer
combined_summary_df <- as.data.frame(combined_summary)

stargazer(
  combined_summary_df,
  type = "text",  # Change to "html" or "latex" for other outputs
  summary = FALSE,  # Disable summary stats for stargazer
  title = "Summary Statistics for DTP3, DTP1, and MCV1 Coverage (2010–2023)",
  rownames = FALSE,
  digits = 2
)


########## GNI per capita summary statistics ################
# Rename dataset 
gni_source <- unique_DTP3   

# Calculate GNI per capita statistics for FGC and NGC
gni_stats <- gni_source %>%
  filter(category %in% c("FGC", "NGC"), year >= 2010 & year <= 2023) %>%  # Filter by category and year
  group_by(category) %>%  # Group by FGC and NGC
  summarize(
    Mean_GNI = mean(gni_2019, na.rm = TRUE),  # Calculate mean GNI per capita
    SD_GNI = sd(gni_2019, na.rm = TRUE),      # Calculate standard deviation
    Observations = n()                        # Count the number of observations
  )

# Convert to a data frame for stargazer
gni_stats_df <- as.data.frame(gni_stats)

# Generate the stargazer-style table
if (!requireNamespace("stargazer", quietly = TRUE)) {
  install.packages("stargazer")
}
library(stargazer)

# Create the table
stargazer(
  gni_stats_df,
  type = "text",  # Change to "html" or "latex" for other outputs
  summary = FALSE,  # Disable summary stats for stargazer
  title = "GNI Per Capita Statistics for FGC and NGC (2010–2023)",
  rownames = FALSE,
  digits = 2
)

########################## DiD FEs##################

data_FGC_NGC_st <- Filtered_DTP3 %>%
  filter(category %in% c("FGC","NGC"),
         year >= 2017, year <= 2021) %>%
  mutate(
    cname          = as.factor(cname),
    treatment      = ifelse(category == "FGC", 1L, 0L),
    post_treatment = ifelse(year >= 2020, 1L, 0L),
    gni_scaled     = (gni_2019 - 4000) / 1000
  )

# FE-DiD (weights applied in the model)
did_fe <- fixest::feols(
  coverage ~ treatment:post_treatment | cname + year,
  data = data_FGC_NGC_st,
  weights = ~ target,
  cluster = ~ cname
)

# show regression summary
summary(did_fe)

# use etable for a cleaner output
etable(
  did_fe,
  se       = "cluster",
  coefstat = "confint",    # show 95% CIs
  digits   = 3,
  dict     = c("treatment:post_treatment" = "Treatment × Post (DiD)"),
  fitstat  = c("n", "r2", "wr2"),
  title    = "Difference-in-Differences with Country and Year Fixed Effects"
)

# FE-DiD (no weights = unweighted country-year averages)
did_fe_unweighted <- fixest::feols(
  coverage ~ treatment:post_treatment | cname + year,
  data = data_FGC_NGC_st,
  cluster = ~ cname
)

# show results
etable(
  did_fe_unweighted,
  se       = "cluster",
  coefstat = "confint",
  digits   = 3,
  dict     = c("treatment:post_treatment" = "Treatment × Post (DiD)"),
  fitstat  = c("n", "r2", "wr2"),
  title    = "DiD (Unweighted, Country and Year Fixed Effects)"
)

###############################################################################
####### FGC vs NGC (Pre:2017-19 vs Post:20-23)_long term_DID analysis with FEs ########## 
###############################################################################

# Long-term (LT) panel: country–year micro data (NOT grouped)
data_FGC_NGC_lt <- Filtered_DTP3 %>%
  dplyr::filter(category %in% c("FGC","NGC"),
                year >= 2017, year <= 2023) %>%
  dplyr::mutate(
    cname          = as.factor(cname),
    treatment      = ifelse(category == "FGC", 1L, 0L),
    post_treatment = ifelse(year >= 2020, 1L, 0L),
    gni_scaled     = (gni_2019 - 4000) / 1000   # optional: for interpretability
  )

# FE-DiD, weighted (child-average effect)
did_fe_lt_w <- fixest::feols(
  coverage ~ treatment:post_treatment | cname + year,
  data     = data_FGC_NGC_lt,
  weights  = ~ target,
  cluster  = ~ cname
)

# FE-DiD, unweighted (country-average effect)
did_fe_lt_u <- fixest::feols(
  coverage ~ treatment:post_treatment | cname + year,
  data     = data_FGC_NGC_lt,
  cluster  = ~ cname
)

# Nicely formatted table
fixest::etable(
  list("LT weighted" = did_fe_lt_w, "LT unweighted" = did_fe_lt_u),
  se       = "cluster",
  coefstat = "confint",
  dict     = c("treatment:post_treatment" = "Treatment × Post (DiD)"),
  fitstat  = c("n","r2","wr2"),
  digits   = 3,
  title    = "DiD (2017–2023): Country & Year Fixed Effects"
)

###### Table presentation##############
etable(
  list("ST Weighted"   = did_fe,
       "ST Unweighted" = did_fe_unweighted,
       "LT Weighted"   = did_fe_lt_w,
       "LT Unweighted" = did_fe_lt_u),
  se       = "cluster",
  coefstat = "confint",   # 95% CIs
  dict     = c("treatment:post_treatment" = "Treatment × Post (DiD)"),
  fitstat  = c("n","r2","wr2"),
  digits   = 2,
  title    = "Difference-in-Differences Estimates (2017–2021 vs 2017–2023)\nWeighted vs Unweighted, Country & Year Fixed Effects"
)

################################################################### 
####################### END #######################################
###################################################################



