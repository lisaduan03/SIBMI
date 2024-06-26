# Load DEMO data for 2015-2016
nhanesTables(data_group='DEMO', year=2015)
demo_1516 <- nhanes("DEMO_I")
Demo_1516 <- nhanesTranslate('DEMO_I', names(demo_1516), data=demo_1516)

# Load LAB data for 2015-2016
nhanesTables(data_group='LAB', year=2015)
cbc_1516 <- nhanes("CBC_I")
biopro_1516 <- nhanes("BIOPRO_I")

# Select relevant columns from DEMO data
Demo_1516select <- Demo_1516[c("SEQN", "RIDEXPRG", "RIAGENDR", "RIDAGEYR", "RIDRETH3")]

# Select relevant columns from LAB data
cbc_1516select <- cbc_1516[c("SEQN", "LBXNEPCT")]
biopro_1516select <- biopro_1516[c("SEQN", "LBXSCR")]

# Merge the data for 2015-2016
merged_1516 <- merge(Demo_1516select, cbc_1516select, by = "SEQN", all = TRUE)
merged_1516 <- merge(merged_1516, biopro_1516select, by = "SEQN", all = TRUE)

# Display the merged data
head(merged_1516)

# repeat for 1718
# Load DEMO data for 2017-2018
nhanesTables(data_group='DEMO', year=2017)
demo_1718 <- nhanes("DEMO_J")
Demo_1718 <- nhanesTranslate('DEMO_J', names(demo_1718), data=demo_1718)

# Load LAB data for 2017-2018
nhanesTables(data_group='LAB', year=2017)
cbc_1718 <- nhanes("CBC_J")
biopro_1718 <- nhanes("BIOPRO_J")

# Select relevant columns from DEMO data
Demo_1718select <- Demo_1718[c("SEQN", "RIDEXPRG", "RIAGENDR", "RIDAGEYR", "RIDRETH3")]

# Select relevant columns from LAB data
cbc_1718select <- cbc_1718[c("SEQN", "LBXNEPCT")]
biopro_1718select <- biopro_1718[c("SEQN", "LBXSCR")]

# Merge the data for 2017-2018
merged_1718 <- merge(Demo_1718select, cbc_1718select, by = "SEQN", all = TRUE)
merged_1718 <- merge(merged_1718, biopro_1718select, by = "SEQN", all = TRUE)

# Display the merged data
head(merged_1718)

# concatenated and display
# Concatenate the two datasets
merged_data_combined <- rbind(merged_1516, merged_1718)

# Display the combined data
head(merged_data_combined)

# Get a summary
summary(merged_data_combined)

# creating a histogram of LBXNEPCT
hist(merged.data$LBXNEPCT,
     main = "histogram of segmented neutrophils present")

with(merged.data, plot(LBXNEPCT ~ RIDAGEYR))

# now, linear regression
lm1 <- with(merged.data, lm(LBXNEPCT ~ RIDAGEYR))
abline(lm1, col="blue", lwd=3)
lm1

# calculate Pearson correlation coefficient
pcc <- with(merged.data, cor(LBXNEPCT, RIDAGEYR, use = "complete.obs"))
pcc

# Segregating by race
race_table <- table(merged.data$RIDRETH3)
barplot(race_table, 
        main = "Frequency Distribution of Race/Hispanic Origin",
        xlab = "Race/Hispanic Origin", 
        ylab = "Frequency",
        col = "skyblue",
        las = 2,  # Rotate the x-axis labels for better readability
        cex.names = 0.8)  # Adjust the size of the x-axis labels

ggplot(merged.data, aes(x = RIDRETH3, y = LBXNEPCT)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Segmented Neutrophils by Race",
       x = "Race", y = "Segmented Neutrophils (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Now, calculating eGFR


