# ... (Sections 1-10 identical to 0.Data_Cleaning.R) ...

# ==================== 2. Feature name processing function ====================
identify_feature_name <- function(feature_name) {
  parts <- str_split(feature_name, "_", simplify = TRUE)
  if (length(parts) >= 3) {
    return(parts[length(parts)])
  }
  return(feature_name)
}

# ==================== 3. Calculate correlation between features and HI ====================
cat("\n========== Calculate correlation between features and HI ==========\n")

# First ensure the HI column exists
if (!"Heterogeneity_Index" %in% colnames(data)) {
  stop("Heterogeneity_Index column not found in data; please calculate HI first")
}

# Extract all features (excluding non-feature columns)
non_feature_cols <- c("Image", "RITH_Phenotype", "Heterogeneity_Index", "Cluster_Label")
all_features <- setdiff(colnames(data), non_feature_cols)

# Calculate Spearman correlation between each feature and HI
cor_results <- data.frame()

for (feat in all_features) {
  cor_val <- cor(data[[feat]], data$Heterogeneity_Index,
                 method = "spearman", use = "complete.obs")

  cor_results <- rbind(cor_results, data.frame(
    Feature = feat,
    Feature_Name = identify_feature_name(feat),
    Correlation = cor_val
  ))
}

# Sort by correlation
cor_results <- cor_results[order(-cor_results$Correlation), ]
rownames(cor_results) <- NULL

# Save correlation results
write_csv(cor_results, "All_Features_HI_Correlation.csv")
cat("Feature-HI correlation analysis saved to: All_Features_HI_Correlation.csv\n")

# ==================== 4. Select top 5 positive and top 5 negative features ====================
cat("\n========== Select top 5 positive and top 5 negative features ==========\n")

# Select top 5 positively correlated features
top_pos_features <- cor_results %>%
  filter(Correlation > 0) %>%
  arrange(desc(Correlation)) %>%
  head(5)

# Select top 5 negatively correlated features (smallest, i.e. most negative)
top_neg_features <- cor_results %>%
  filter(Correlation < 0) %>%
  arrange(Correlation) %>%
  head(5)

cat("\nTop 5 positively correlated features:\n")
print(top_pos_features)

cat("\nTop 5 negatively correlated features:\n")
print(top_neg_features)

# Combine selected features
selected_features <- rbind(top_pos_features, top_neg_features)
selected_feature_names <- selected_features$Feature_Name
selected_feature_codes <- selected_features$Feature

# ==================== 5. Prepare plotting data (using full feature names) ====================
cat("\n========== Prepare plotting data (using full feature names) ==========\n")

# Add complete information to selected_features
selected_features <- selected_features %>%
  mutate(
    Filter_Type = sapply(Feature, identify_filter_type),
    Feature_Class = sapply(Feature, identify_feature_class),
    Base_Feature_Name = sapply(Feature, identify_feature_name),
    Full_Feature_Name = Feature,
    Display_Full = ifelse(
      Filter_Type == "original",
      paste0(Base_Feature_Name, " (Original)"),
      paste0(Base_Feature_Name, " (", Filter_Type, ")")
    ),
    Display_Full = paste0(Display_Full, " [", Feature_Class, "]")
  )

# Check uniqueness
cat("Check uniqueness of full feature names:\n")
dup_check <- selected_features %>%
  group_by(Display_Full) %>%
  summarise(count = n()) %>%
  filter(count > 1)

if (nrow(dup_check) > 0) {
  cat("Duplicate display names found; add index:\n")
  print(dup_check)

  selected_features <- selected_features %>%
    group_by(Display_Full) %>%
    mutate(
      dup_num = row_number(),
      Display_Full = ifelse(dup_num > 1,
                            paste0(Display_Full, " #", dup_num),
                            Display_Full)
    ) %>%
    ungroup()
}

# Create plotting data frame
plot_data <- data.frame(
  HI = data$Heterogeneity_Index,
  RITH_Phenotype = data$RITH_Phenotype
)

# Add selected features
for (i in 1:nrow(selected_features)) {
  feat_code <- selected_features$Feature[i]
  plot_data[[feat_code]] <- data[[feat_code]]
}

cat("Plotting data prepared, containing", ncol(plot_data) - 2, "features\n")

# ==================== 6. Create three types of display names ====================
cat("\n========== Create three types of display names ==========\n")

selected_features <- selected_features %>%
  mutate(
    Display_Name1 = paste0(
      Display_Full,
      "\nr = ", round(Correlation, 3),
      ifelse(Correlation > 0, " (+) ", " (-) ")
    ),
    Display_Name2 = paste0(
      ifelse(Correlation > 0, "↗ ", "↘ "),
      Display_Full,
      "\n(r = ", round(Correlation, 3), ")"
    ),
    Display_Name3 = paste0(
      "Feature ", sprintf("%02d", row_number()), ": ",
      Base_Feature_Name, "\n",
      "Filter: ", Filter_Type,
      " | Class: ", Feature_Class,
      "\nr = ", round(Correlation, 3),
      ifelse(Correlation > 0, " (Positive)", " (Negative)")
    )
  )

# Print feature information
cat("\nSelected feature information:\n")
for (i in 1:nrow(selected_features)) {
  cat("\n", i, ". ", selected_features$Feature[i], "\n", sep = "")
  cat("   Short name: ", selected_features$Base_Feature_Name[i], "\n", sep = "")
  cat("   Filter: ", selected_features$Filter_Type[i], "\n", sep = "")
  cat("   Class: ", selected_features$Feature_Class[i], "\n", sep = "")
  cat("   Correlation: ", round(selected_features$Correlation[i], 3), "\n", sep = "")
}

# ==================== 7. Convert data to long format ====================
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = all_of(selected_features$Feature),
    names_to = "Feature_Code",
    values_to = "Feature_Value"
  ) %>%
  left_join(selected_features, by = c("Feature_Code" = "Feature"))

# ==================== 8. Function to remove extreme values ====================
cat("\n========== Remove extreme values ==========\n")

remove_outliers <- function(x, lower_quantile = 0.025, upper_quantile = 0.975) {
  lower_bound <- quantile(x, lower_quantile, na.rm = TRUE)
  upper_bound <- quantile(x, upper_quantile, na.rm = TRUE)
  x[x < lower_bound] <- NA
  x[x > upper_bound] <- NA
  return(x)
}

plot_data_clean <- plot_data_long %>%
  group_by(Feature_Code) %>%
  mutate(
    Feature_Value_Clean = remove_outliers(Feature_Value)
  ) %>%
  ungroup() %>%
  filter(!is.na(Feature_Value_Clean))

# ==================== 9. Method 1: basic scatter plot (full feature names) ====================
cat("\n========== Method 1: basic scatter plot ==========\n")

# Sort by correlation
plot_data_clean$Display_Name1 <- factor(
  plot_data_clean$Display_Name1,
  levels = unique(plot_data_clean$Display_Name1[order(-plot_data_clean$Correlation)])
)

p_method1 <- ggplot(plot_data_clean,
                    aes(x = Feature_Value_Clean, y = HI, color = RITH_Phenotype)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", color = "black", se = FALSE, size = 0.7) +
  facet_wrap(~ Display_Name1, scales = "free_x", ncol = 2) +
  scale_color_manual(values = c("Low RITH" = "#377EB8", "High RITH" = "#E41A1C")) +
  labs(
    title = "",
    subtitle = paste("Top", nrow(selected_features), "features correlated with Heterogeneity Index"),
    x = "Feature Value (95% central range)",
    y = "Heterogeneity Index (HI)",
    color = "RITH Phenotype"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11, margin = margin(b = 15)),
    strip.text = element_text(size = 9, face = "bold", lineheight = 1.1),
    strip.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines")
  )

# Save figure
ggsave("HI_Scatter_Method1_FullNames.pdf", p_method1,
       width = 16, height = 12, dpi = 300)
cat("Method 1 figure saved: HI_Scatter_Method1_FullNames.pdf\n")

# ==================== 12. Save feature information table (fixed version) ====================
cat("\n========== Save feature information ==========\n")

# First use mutate to create new columns, then select
feature_info_table <- selected_features %>%
  mutate(
    Direction = ifelse(Correlation > 0, "Positive", "Negative"),
    Correlation_Strength = case_when(
      abs(Correlation) >= 0.7 ~ "Very Strong",
      abs(Correlation) >= 0.5 ~ "Strong",
      abs(Correlation) >= 0.3 ~ "Moderate",
      TRUE ~ "Weak"
    )
  ) %>%
  select(
    Feature_Code = Feature,
    Base_Name = Base_Feature_Name,
    Filter_Type,
    Feature_Class,
    Correlation,
    Direction,
    Correlation_Strength
  ) %>%
  arrange(desc(Correlation))

write_csv(feature_info_table, "Selected_Features_Detailed_Info.csv")
cat("Detailed feature information saved: Selected_Features_Detailed_Info.csv\n")

# Display feature information
cat("\nDetailed list of selected features:\n")
print(feature_info_table)
