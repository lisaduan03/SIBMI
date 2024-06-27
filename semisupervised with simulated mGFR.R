# loading in packages
library(nhanesA)
library(survey)
library(ggplot2)
library(dplyr)
library(stringr)
library(gridExtra)

##### Loading in data ######
# load DEMO data for 1999-2000
nhanesTables(data_group='DEMO', year=1999)
demo_9900 <- nhanes("DEMO")
Demo_9900 <- nhanesTranslate('DEMO', names(demo_9900), data=demo_9900)
# load LAB data for 1999-2000
biopro_9900 <- nhanes("LAB18") # the name of biochemistry profile
cyst_9900 <- nhanes("SSCYST_A")
# selecting columns
Demo_9900select <- Demo_9900[c("SEQN", "RIDEXPRG", "RIAGENDR", "RIDAGEYR", "RIDRETH2", "SDMVPSU", "SDMVSTRA", "WTMEC4YR", "SDDSRVYR", "RIDSTATR")] # should I use ETH1 or 2?
# Correction is recommended: The Deming regression (adjusting for errors in measurement) for the correction is Standard Creatinine (Y) = 1.013*NHANES Creatinine (X) + 0.147 (r = 0.984).
biopro_9900select <- biopro_9900[c("SEQN", "LBXSCR")]
cyst_9900select <- cyst_9900[c("SEQN", "SSCYPC", "WTSCY4YR")] #  Surplus sera cystatin 99-02 weights
# merging on ID (SEQN)
merged_9900 <- merge(Demo_9900select, biopro_9900select, by = "SEQN", all = TRUE)
merged_9900 <- merge(merged_9900, cyst_9900select, by = "SEQN", all = TRUE)

# standardizing serum creatinine
merged_9900$Standard_Creatinine <- 1.013 * merged_9900$LBXSCR + 0.147 # adjustment
merged_9900 <- merged_9900[, !(names(merged_9900) %in% "LBXSCR")] # drop original LBXSCR
merged_9900 <- merged_9900 %>% 
  rename(LBXSCR = Standard_Creatinine)

head(merged_9900)

# repeat for 0102
nhanesTables(data_group='LAB', year=2001)
demo_0102 <- nhanes("DEMO_B")
Demo_0102 <- nhanesTranslate('DEMO_B', names(demo_0102), data=demo_0102)
nhanesTables(data_group='LAB', year=2001)
biopro_0102 <- nhanes("L40_B")
cyst_0102 <- nhanes("SSCYST_B")
Demo_0102select <- Demo_0102[c("SEQN", "RIDEXPRG", "RIAGENDR", "RIDAGEYR", "RIDRETH2", "SDMVPSU", "SDMVSTRA", "WTMEC4YR", "SDDSRVYR", "RIDSTATR")]
biopro_0102select <- biopro_0102[c("SEQN", "LBDSCR")]
cyst_0102select <- cyst_0102[c("SEQN", "SSCYPC", "WTSCY4YR")] #  Surplus sera cystatin 99-02 weights
merged_0102 <- merge(Demo_0102select, biopro_0102select, by = "SEQN", all = TRUE)
merged_0102 <- merge(merged_0102, cyst_0102select, by = "SEQN", all = TRUE)
merged_0102 <- merged_0102 %>% 
  rename(LBXSCR = LBDSCR)
head(merged_0102)


# concatenate and display
merged_0102 <- merged_0102[, names(merged_9900)] # reordering columns in 0102 to match
merged_data_combined <- rbind(merged_9900, merged_0102)
head(merged_data_combined)

# Get a summary
summary(merged_data_combined)


#### trying to replicate Diao et al., 2021 JAMA (reply) #####
# get missing counts
missing_count_scr <- sum(is.na(merged_data_combined$LBXSCR))
print(paste("Number of missing serum creatinine values:", missing_count_scr))
print(paste("Number of total people:", nrow(merged_data_combined)))
print(paste("Response rate:", 1 - (missing_count_scr/nrow(merged_data_combined))))

missing_count_cys <- sum(is.na(merged_data_combined$SSCYPC))
print(paste("Number of missing serum cys c values:", missing_count_cys))
print(paste("Number of total people:", nrow(merged_data_combined)))
print(paste("Response rate:", 1 - (missing_count_cys/nrow(merged_data_combined))))

# filtering for complete data, adults, non pregnant. 4262 vs 4434 in study, but lieklye due to 2022 data removal of SSCYPC
merged_data_sub <- subset(merged_data_combined, complete.cases(SSCYPC, LBXSCR) 
                          & !RIDEXPRG %in% c("Yes, positive lab pregnancy test or self-reported pregnant at exam")
                          & RIDAGEYR > 18)


### adding eGFR ####
calculate_eGFR <- function(creatinine, age, sex, race, with_race) {
  if (sex == 1) {  # Male
    if (creatinine <= 0.9) {
      eGFR <- 141 * (creatinine / 0.9) ^ -0.411 * (0.993 ^ age)
    } else {
      eGFR <- 141 * (creatinine / 0.9) ^ -1.209 * (0.993 ^ age)
    }
  } else {  # Female
    if (creatinine <= 0.7) {
      eGFR <- 144 * (creatinine / 0.7) ^ -0.329 * (0.993 ^ age)
    } else {
      eGFR <- 144 * (creatinine / 0.7) ^ -1.209 * (0.993 ^ age)
    }
  }
  # Adjust for Black race
  if (race == "Non-Hispanic Black" && with_race) {
    print("NHB!!!!!")
    eGFR <- eGFR * 1.159
  }
  return(eGFR)
}

### adding columns to the raw data and weighted
merged_data_sub <- merged_data_sub %>%
  mutate(
    eGFR = mapply(calculate_eGFR, creatinine = LBXSCR, age = RIDAGEYR, sex = RIAGENDR, race = RIDRETH2, with_race = TRUE),
    eGFR_no_race = mapply(calculate_eGFR, creatinine = LBXSCR, age = RIDAGEYR, sex = RIAGENDR, race = RIDRETH2, with_race = FALSE),
    change_in_eGFR = eGFR_no_race - eGFR
  )



### Survey Weights ####
## figure out the other weights!!!! WTSCY4YR
nhanesDesign <- svydesign(id = ~SDMVPSU,  # Primary Sampling Units (PSU)
                          strata  = ~SDMVSTRA, # Stratification used in the survey
                          weights = ~WTMEC4YR,  # according to CDC Weighting Module, use this for 1999-2002 
                          nest    = TRUE,      # Whether PSUs are nested within strata
                          data    = merged_data_sub)

# filtering. 4262 compared to Diao 4434 likely because in April 2022 some records removed from SSCYPC
#dfsub <- subset(nhanesDesign , RIDSTATR %in% "Both Interviewed and MEC examined") # this line useless, caught by next line
#dfsub <- subset(nhanesDesign, complete.cases(SSCYPC, LBXSCR) 
        #        & !RIDEXPRG %in% c("Yes, positive lab pregnancy test or self-reported pregnant at exam")
         #       & RIDAGEYR > 18)

### descriptive stats ### matches Diao et al., 2021 (JAMA reply)

unweighted_proportions_race <- merged_data_sub %>%
  count(RIDRETH2) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  select(RIDRETH2, proportion)
print(unweighted_proportions_race)

unweighted_proportions_gender <- merged_data_sub %>%
  count(RIAGENDR) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  select(RIAGENDR, proportion)
print(unweighted_proportions_gender)


median_age <- median(merged_data_sub$RIDAGEYR, na.rm = TRUE)
print(median_age)

svymean( ~ RIDRETH2 , nhanesDesign )

## doesn't match up below
merged_data_sub_black <- merged_data_sub %>%
  filter(RIDRETH2 == "Non-Hispanic Black")
mean_eGFR <- mean(merged_data_sub$eGFR, na.rm = TRUE)
mean_eGFR_no_race <- mean(merged_data_sub$eGFR_no_race, na.rm = TRUE)


dfsub_black = subset(nhanesDesign, RIDRETH2 == "Non-Hispanic Black")


svymean( ~ eGFR , dfsub_black )

svymean( ~ eGFR_no_race , dfsub_black )


weighted_median_change_in_eGFR <- svyquantile(~ change_in_eGFR, design = dfsub_black, quantiles = 0.5)
print(weighted_median_change_in_eGFR)



###### Recreating Figure 1 in Diao et al., 2020 JAMA Paper ######

# without survery weights
merged_data_sub_black <- merged_data_sub %>%
  filter(RIDRETH2 == "Non-Hispanic Black")
ggplot(merged_data_sub_black, aes(x = change_in_eGFR)) +
  geom_histogram(binwidth = 1, fill = "salmon", color = "black", alpha = 0.7) +
  labs(x = expression(paste("change in eGFR ", (mL/min/1.73 ~ m^2), " (no race - race)")),
       y = "Black adults in NHANES (2015-2018)",
       title = "Changes in Reported eGFR for Black Adults Following Removal of Race From eGFRcr") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold")
  )

# Histogram with survery weights
dfsub_black = subset(nhanesDesign, RIDRETH2 == "Non-Hispanic Black")
svyhist(~ change_in_eGFR, design = dfsub_black,
        main = "Changes in Reported eGFR for Black Adults Following Removal of Race From eGFRcr",
        xlab = expression(paste("change in eGFR ", (mL/min/1.73 ~ m^2), " (no race - race)")),
        ylab = "Black adults in NHANES (2015-2018)",
        col = "green",
        border = "black",
        cex.main = 0.8)


######## Table 1 ########
# function that takes in cutoff ranges and implication name, creates row in a table
create_summary_table <- function(cutoff_ranges, implications) {
  all_tables <- list()
  
  for (i in seq_along(cutoff_ranges)) {
    cutoff_range <- cutoff_ranges[[i]]
    implication <- implications[[i]]
    
    lower_cutoff <- cutoff_range[1]
    upper_cutoff <- cutoff_range[2]
    
    # Weighted: proportions with race
    prop_and_se_with_race <- svyby(~(eGFR >= lower_cutoff & eGFR < upper_cutoff), ~RIDRETH2, design = dfsub_black, svymean, na.rm = TRUE) # mean of many 1s and 0s = proportion 
    prop_with_race <- prop_and_se_with_race["eGFR >= lower_cutoff & eGFR < upper_cutoffTRUE"]
    se_with_race <- prop_and_se_with_race["se.eGFR >= lower_cutoff & eGFR < upper_cutoffTRUE"]
    
    # Weighted:  proportions without race
    prop_and_se_without_race <- svyby(~(eGFR_no_race >= lower_cutoff & eGFR_no_race < upper_cutoff), ~RIDRETH2, design = dfsub_black, svymean, na.rm = TRUE)
    prop_without_race <- prop_and_se_without_race[prop_and_se_without_race$RIDRETH2 == "Non-Hispanic Black", "eGFR_no_race >= lower_cutoff & eGFR_no_race < upper_cutoffTRUE"]
    se_without_race <- prop_and_se_without_race[prop_and_se_without_race$RIDRETH2 == "Non-Hispanic Black", "se.eGFR_no_race >= lower_cutoff & eGFR_no_race < upper_cutoffTRUE"]
    
    # Unweighted: counts with and without race
    count_with_race <- merged_data_sub_black %>%
      summarise(n = sum(eGFR >= lower_cutoff & eGFR < upper_cutoff, na.rm = TRUE))
    
    count_without_race <- merged_data_sub_black %>%
      summarise(n = sum(eGFR_no_race >= lower_cutoff & eGFR_no_race < upper_cutoff, na.rm = TRUE))
    
    # Weighted: Calculate absolute change in proportions
    absolute_change <- prop_without_race - prop_with_race
    
    # Weighted: formula SE of difference of two proportions
    se_absolute_change <- sqrt(se_with_race^2 + se_without_race^2)
    
    # CIs
    ci_low <- absolute_change - 1.96 * se_absolute_change
    ci_high <- absolute_change + 1.96 * se_absolute_change
    
    # Wrap long text text to fit in the table
    wrapped_implication <- str_wrap(implication, width = 40)  # Adjust width as needed
    
    # creating row in table 
    data_table <- data.frame(
      "Implication" = wrapped_implication,
      "eGFR range" = paste(">=", lower_cutoff, "& <", upper_cutoff),
      "Including Race" = paste0(count_with_race$n, " (", format(prop_with_race * 100, nsmall = 2), ")"),
      "Removing Race" = paste0(count_without_race$n, " (", format(prop_without_race * 100, nsmall = 2), ")"),
      "Absolute change, weighted % (95% CI)" = paste0(format(absolute_change * 100, nsmall = 2),
                                                      " (", format(ci_low * 100, nsmall = 2), ", ", format(ci_high * 100, nsmall = 2), ")")    )
    
    
    # Append current data table to list
    all_tables[[i]] <- data_table
  }
  
  # Combine all data tables into one
  combined_table <- do.call(rbind, all_tables)
  
  return(combined_table)
}

# Example usage
cutoff_ranges <- list(c(0, 60), c(0, 45), c(13, 50), c(0, 30), c(15, 29), c(0, 20))
implications <- c("CKD Diagnosis", "CKD stage 3b or higher and related drug recommendations", 
                  "Medical nutrition therapy covered", "Reclassification: CKD stage 4 or higher and related drug recommendations", 
                  "Kidney disease education covered", "Eligible for kidney transplant waiting list")

summary_table <- create_summary_table(cutoff_ranges, implications)
print(summary_table)

myTable <- tableGrob(summary_table)
library(grid)
grid.draw(myTable)

# controls number of digits
options(digits=2)

nhanes_srvyr_design <- as_survey(dfsub_black)
nhanes_srvyr_design %>%
  summarize(proportion = survey_mean(eGFR_no_race < 60, na.rm = TRUE))


