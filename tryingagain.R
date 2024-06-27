# loading in packages
library(nhanesA)
library(survey)
library(ggplot2)
library(dplyr)
library(stringr)
library(gridExtra)

##### Loading in data ######
# load DEMO data for 2015-2016
nhanesTables(data_group='DEMO', year=2015)
demo_1516 <- nhanes("DEMO_I")
Demo_1516 <- nhanesTranslate('DEMO_I', names(demo_1516), data=demo_1516)
# load LAB data for 2015-2016
nhanesTables(data_group='LAB', year=2015)
biopro_1516 <- nhanes("BIOPRO_I")
# selecting columns
Demo_1516select <- Demo_1516[c("SEQN", "RIDEXPRG", "RIAGENDR", "RIDAGEYR", "RIDRETH3", "SDMVPSU", "SDMVSTRA", "WTMEC2YR", "SDDSRVYR")]
biopro_1516select <- biopro_1516[c("SEQN", "LBXSCR")]
# merging on ID (SEQN)
merged_1516 <- merge(Demo_1516select, biopro_1516select, by = "SEQN", all = TRUE)
# display
head(merged_1516)

# repeat for 1718
nhanesTables(data_group='DEMO', year=2017)
demo_1718 <- nhanes("DEMO_J")
Demo_1718 <- nhanesTranslate('DEMO_J', names(demo_1718), data=demo_1718)
nhanesTables(data_group='LAB', year=2017)
biopro_1718 <- nhanes("BIOPRO_J")
Demo_1718select <- Demo_1718[c("SEQN", "RIDEXPRG", "RIAGENDR", "RIDAGEYR", "RIDRETH3", "SDMVPSU", "SDMVSTRA", "WTMEC2YR", "SDDSRVYR")]
biopro_1718select <- biopro_1718[c("SEQN", "LBXSCR")]
merged_1718 <- merge(Demo_1718select, biopro_1718select, by = "SEQN", all = TRUE)

# Concatenate the two datasets
merged_data_combined <- rbind(merged_1516, merged_1718)

# Display the combined data
head(merged_data_combined)

# Get a summary
summary(merged_data_combined)
# missing_indices <- which(is.na(merged_data_combined$LBXSCR))

# Mean imputation
# mean_cre <- mean(merged_data_combined$LBXSCR, na.rm = TRUE)
# merged_data_combined$LBXSCR[missing_indices] <- mean_cre

#### [lltting] ###
calculate_eGFR_no_race <- function(creatinine, age, sex) {
  if (is.na(creatinine) | is.na(age) | is.na(sex)) {
    return(NA)
  }
    
  if (sex == 1) {  # Male
    if (creatinine <= 0.9) {
      eGFR <- 141 * (creatinine / 0.9) ^ -0.411 * 0.993 ^ age
    } else {
      eGFR <- 141 * (creatinine / 0.9) ^ -1.209 * 0.993 ^ age
    }
  } else {  # Female
    if (creatinine <= 0.7) {
      eGFR <- 144 * (creatinine / 0.7) ^ -0.329 * 0.993 ^ age
    } else {
      eGFR <- 144 * (creatinine / 0.7) ^ -1.209 * 0.993 ^ age
    }
  }
  return(eGFR)
}

# Calculate eGFR
calculate_eGFR <- function(creatinine, age, sex, race) {
  if (is.na(creatinine) | is.na(age) | is.na(sex)) {
    return(NA)
  }
    
  if (sex == 1) {  # Male
    if (creatinine <= 0.9) {
      eGFR <- 141 * (creatinine / 0.9) ^ -0.411 * 0.993 ^ age
    } else {
      eGFR <- 141 * (creatinine / 0.9) ^ -1.209 * 0.993 ^ age
    }
  } else {  # Female
    if (creatinine <= 0.7) {
      eGFR <- 144 * (creatinine / 0.7) ^ -0.329 * 0.993 ^ age
    } else {
      eGFR <- 144 * (creatinine / 0.7) ^ -1.209 * 0.993 ^ age
    }
  }
  # Adjust for Black race
  if (race == "Non-Hispanic Black") {
    print("NHB!!!!!")
    eGFR <- eGFR * 1.159
  }
  return(eGFR)
}

merged_data_combined$eGFR <- mapply(calculate_eGFR, 
                                    creatinine = merged_data_combined$LBXSCR, 
                                    age = merged_data_combined$RIDAGEYR, 
                                    sex = merged_data_combined$RIAGENDR, 
                                    race = merged_data_combined$RIDRETH3)

# Apply eGFR calculation without race adjustment
merged_data_combined$eGFR_no_race <- mapply(calculate_eGFR_no_race, 
                                            creatinine = merged_data_combined$LBXSCR, 
                                            age = merged_data_combined$RIDAGEYR, 
                                            sex = merged_data_combined$RIAGENDR)

merged_data_combined <- merged_data_combined %>%
  mutate(change_in_eGFR = eGFR_no_race - eGFR)

# Calculate mean eGFR without race adjustment
mean_egfr_no_race <- merged_data_combined %>%
  group_by(RIDRETH3) %>%
  summarise(mean_eGFR_no_race = mean(eGFR_no_race, na.rm = TRUE))

# Plot without race adjustment
plot_no_race <- ggplot(merged_data_combined, aes(x = factor(RIDRETH3), y = eGFR_no_race)) +
  geom_boxplot() +
  geom_text(data = mean_egfr_no_race, aes(x = factor(RIDRETH3), y = mean_eGFR_no_race,
                                          label = sprintf("%.1f", mean_eGFR_no_race)),
            vjust = -0.5, size = 4, color = "black", fontface = "bold") +
  labs(x = "Race/Ethnicity", y = "eGFR", title = "Distribution of eGFR without Race Adjustment") +
  scale_x_discrete(labels = c("Mexican American", "Other Hispanic", "Non-Hispanic White",
                              "Non-Hispanic Black", "Other Race/Multiracial")) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate mean eGFR with race adjustment
mean_egfr <- merged_data_combined %>%
  group_by(RIDRETH3) %>%
  summarise(mean_eGFR = mean(eGFR, na.rm = TRUE))

# Plot with race adjustment
plot_with_race <- ggplot(merged_data_combined, aes(x = factor(RIDRETH3), y = eGFR)) +
  geom_boxplot() +
  geom_text(data = mean_egfr, aes(x = factor(RIDRETH3), y = mean_eGFR,
                                  label = sprintf("%.1f", mean_eGFR)),
            vjust = -0.5, size = 4, color = "black", fontface = "bold") +
  labs(x = "Race/Ethnicity", y = "eGFR", title = "Distribution of eGFR with Race Adjustment") +
  scale_x_discrete(labels = c("Mexican American", "Other Hispanic", "Non-Hispanic White",
                              "Non-Hispanic Black", "Other Race/Multiracial")) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Arrange plots side by side using gridExtra
grid.arrange(plot_no_race, plot_with_race, ncol = 2)


nhanesDesign <- svydesign(id = ~SDMVPSU,  # Primary Sampling Units (PSU)
                          strata  = ~SDMVSTRA, # Stratification used in the survey
                          weights = ~ifelse(SDDSRVYR %in% c(9, 10), 0.5 * WTMEC2YR, WTMEC2YR),  # according to CDC Weighting Module
                          nest    = TRUE,      # Whether PSUs are nested within strata
                          data    = merged_data_combined)

datasub <- merged_data_combined %>%
  filter(RIDRETH3 == "Non-Hispanic Black" & 
           !RIDEXPRG %in% c("Cannot ascertain if the participant is pregnant at exam",
                            "Yes, positive lab pregnancy test or self-reported pregnant at exam"))


###### Recreating Figure 1 in Diao et al., 2020 JAMA Paper ######
datasub <- merged_data_combined %>%
  filter(RIDRETH3 == "Non-Hispanic Black" & 
           !RIDEXPRG %in% c("Cannot ascertain if the participant is pregnant at exam",
                            "Yes, positive lab pregnancy test or self-reported pregnant at exam"))

missing_count <- sum(is.na(datasub$LBXSCR))
print(paste("Number of missing serum creatinine values:", missing_count))
print(paste("Number of total people:", nrow(datasub)))
print(paste("Response rate:", 1 - (missing_count/nrow(datasub))))

# create a subset of nhanesDesign for NHB and nonpregnant adults
dfsub <- subset(nhanesDesign, RIDRETH3 == "Non-Hispanic Black" & 
                  !RIDEXPRG %in% c("Cannot ascertain if the participant is pregnant at exam",
                                   "Yes, positive lab pregnancy test or self-reported pregnant at exam"))


# without survery weights
ggplot(datasub, aes(x = change_in_eGFR)) +
  geom_histogram(binwidth = 1, fill = "salmon", color = "black", alpha = 0.7) +
  labs(x = expression(paste("change in eGFR ", (mL/min/1.73 ~ m^2), " (no race - race)")),
       y = "Black adults in NHANES (2015-2018)",
       title = "Changes in Reported eGFR for Black Adults Following Removal of Race From eGFRcr") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold")
  )

# Histogram with survery weights
svyhist(~ change_in_eGFR, design = dfsub,
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
    prop_and_se_with_race <- svyby(~(eGFR >= lower_cutoff & eGFR < upper_cutoff), ~RIDRETH3, design = nhanesDesign, svymean, na.rm = TRUE) # mean of many 1s and 0s = proportion 
    prop_with_race <- prop_and_se_with_race[prop_and_se_with_race$RIDRETH3 == "Non-Hispanic Black", "eGFR >= lower_cutoff & eGFR < upper_cutoffTRUE"]
    se_with_race <- prop_and_se_with_race[prop_and_se_with_race$RIDRETH3 == "Non-Hispanic Black", "se.eGFR >= lower_cutoff & eGFR < upper_cutoffTRUE"]
    
    # Weighted:  proportions without race
    prop_and_se_without_race <- svyby(~(eGFR_no_race >= lower_cutoff & eGFR_no_race < upper_cutoff), ~RIDRETH3, design = nhanesDesign, svymean, na.rm = TRUE)
    prop_without_race <- prop_and_se_without_race[prop_and_se_without_race$RIDRETH3 == "Non-Hispanic Black", "eGFR_no_race >= lower_cutoff & eGFR_no_race < upper_cutoffTRUE"]
    se_without_race <- prop_and_se_without_race[prop_and_se_without_race$RIDRETH3 == "Non-Hispanic Black", "se.eGFR_no_race >= lower_cutoff & eGFR_no_race < upper_cutoffTRUE"]
    
    # Unweighted: counts with and without race
    count_with_race <- datasub %>%
      summarise(n = sum(eGFR >= lower_cutoff & eGFR < upper_cutoff, na.rm = TRUE))
    
    count_without_race <- datasub %>%
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
