# Load DEMO data for 2015-2016
nhanesTables(data_group='DEMO', year=2015)
demo_1516 <- nhanes("DEMO_I")
Demo_1516 <- nhanesTranslate('DEMO_I', names(demo_1516), data=demo_1516)

# Load LAB data for 2015-2016
nhanesTables(data_group='LAB', year=2015)
cbc_1516 <- nhanes("CBC_I")
biopro_1516 <- nhanes("BIOPRO_I")

# Select relevant columns from data
Demo_1516select <- Demo_1516[c("SEQN", "RIDEXPRG", "RIAGENDR", "RIDAGEYR", "RIDRETH3")]
cbc_1516select <- cbc_1516[c("SEQN", "LBXHCT", "LBXHGB", "LBXRBCSI")] # hematocrit, hemoglobin, and red blood cell count
biopro_1516select <- biopro_1516[c("SEQN",  "LBXSAPSI")]

# Merge the data for 2015-2016
merged_1516 <- merge(Demo_1516select, cbc_1516select, by = "SEQN", all = TRUE)

# Display the merged data
head(merged_1516)
merged_1516 <- merge(merged_1516, biopro_1516select, by = "SEQN", all = TRUE)


# repeat for 1718
# Load DEMO data for 2017-2018
nhanesTables(data_group='DEMO', year=2017)
demo_1718 <- nhanes("DEMO_J")
Demo_1718 <- nhanesTranslate('DEMO_J', names(demo_1718), data=demo_1718)

# Load LAB data for 2017-2018
nhanesTables(data_group='LAB', year=2017)
cbc_1718 <- nhanes("CBC_J")
biopro_1718 <- nhanes("BIOPRO_I")

# Select relevant columns from  data
Demo_1718select <- Demo_1718[c("SEQN", "RIDEXPRG", "RIAGENDR", "RIDAGEYR", "RIDRETH3")]
cbc_1718select <- cbc_1718[c("SEQN", "LBXHCT", "LBXHGB", "LBXRBCSI")]
biopro_1718select <- biopro_1718[c("SEQN",  "LBXSAPSI")]

# Merge the data for 2017-2018
merged_1718 <- merge(Demo_1718select, cbc_1718select, by = "SEQN", all = TRUE)
merged_1718 <- merge(merged_1718, biopro_1718select, by = "SEQN", all = TRUE)

# Display the merged data
head(merged_1718)

# concatenated and display
merged_data_combined <- rbind(merged_1516, merged_1718)

head(merged_data_combined)

# plotting

# Plotting LBXHCT vs Age
plot_hct <- ggplot(merged_data_combined, aes(x = RIDAGEYR, y = LBXHCT, color = RIAGENDR)) +
  geom_point() +
  labs(x = "Age", y = "LBXHCT", color = "Gender") +
  ggtitle("Scatterplot of LBXHCT vs Age") +
  theme_minimal()

# Plotting LBXHGB vs Age
plot_hgb <- ggplot(merged_data_combined, aes(x = RIDAGEYR, y = LBXHGB, color = RIAGENDR)) +
  geom_point() +
  labs(x = "Age", y = "LBXHGB", color = "Gender") +
  ggtitle("Scatterplot of LBXHGB vs Age") +
  theme_minimal()

# Plotting LBXRBCSI vs Age
plot_rbcsi <- ggplot(merged_data_combined, aes(x = RIDAGEYR, y = LBXRBCSI, color = RIAGENDR)) +
  geom_point() +
  labs(x = "Age", y = "LBXRBCSI", color = "Gender") +
  ggtitle("Scatterplot of LBXRBCSI vs Age") +
  theme_minimal()

plot_hct
plot_hgb
plot_rbcsi

# Plot alkaline phosphatase
plot_apsi <- ggplot(merged_data_combined, aes(x = RIDAGEYR, y = LBXSAPSI)) +
  geom_point() +
  labs(x = "Age", y = "LBXSAPSI") +
  ggtitle("Scatterplot of LBXSAPSI vs Age") +
  theme_minimal()

plot_apsi
