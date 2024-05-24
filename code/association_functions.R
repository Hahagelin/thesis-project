# Association analysis functions

# Function for continuous variables
DEA_func <- function(df, clin_df, variable_string) {
  # Merge clinical and gene expression data based on sample ID
  clin_df <- clin_df |> 
    select(all_of(variable_string), sample_id)
  
  merged_df <- df |>
    left_join(clin_df, by = "sample_id")
  
  # Prepare data for linear regression
  outcome <- merged_df[[variable_string]]
  predictors <- merged_df %>% 
    select(-all_of(variable_string), -sample_id) %>% 
    data.frame()
  
  # Initialize results data frame
  results <- data.frame(Protein = colnames(predictors), Coefficient = NA, PValue = NA, stringsAsFactors = FALSE)
  
  # Loop through each gene to perform linear regression
  for (i in seq_along(colnames(predictors))) {
    npx <- predictors[, i]
    model <- glm(outcome ~ npx)  # Linear regression (default family is gaussian)
    summary_model <- summary(model)
    
    # Extract coefficient and p-value for the gene
    predictor_name <- colnames(predictors)[i]
    results$Estimate[i] <- summary_model$coefficients["npx", "Estimate"]
    results$PValue[i] <- summary_model$coefficients["npx", "Pr(>|t|)"]
    
  }
  # Adjust p-values for multiple testing, for example, using Benjamini-Hochberg method
  results$AdjPValue <- p.adjust(results$PValue, method = "BH")
  
  # Return the results
  return(results)
}

# For binary variables
DEA_func_logistic <- function(df, clin_df, variable_string) {
  
  # Merge clinical and gene expression data based on sample ID
  clin_df <- clin_df |> 
    select(all_of(variable_string), sample_id)
  
  merged_df <- df |>
    left_join(clin_df, by = "sample_id")
  
  # Prepare data for logistic regression
  outcome <- as.factor(merged_df[[variable_string]]) # Ensure the outcome is a factor
  predictors <- merged_df %>% 
    select(-all_of(variable_string), -sample_id) %>% 
    data.frame()
  
  # Initialize results data frame
  results <- data.frame(Protein = colnames(predictors), Coefficient = NA, PValue = NA, stringsAsFactors = FALSE)
  
  # Loop through each gene to perform logistic regression
  for (i in seq_along(colnames(predictors))) {
    npx <- predictors[, i]
    model <- glm(outcome ~ npx, family = "binomial")
    summary_model <- summary(model)
    
    # Extract coefficient and p-value for the gene
    results$Coefficient[i] <- summary_model$coefficients["npx", "Estimate"]
    results$PValue[i] <- summary_model$coefficients["npx", "Pr(>|z|)"]
  }
  
  # Adjust p-values for multiple testing, for example, using Benjamini-Hochberg method
  results$AdjPValue <- p.adjust(results$PValue, method = "BH")
  
  # Return the results
  return(results)
}

# Function for ordinal variables
clm_maker <- function(df, clin_df, variable_string) {
  
  # Merge clinical and gene expression data based on sample ID
  clin_df <- clin_df |> 
    select(all_of(variable_string), sample_id)
  
  merged_df <- df |>
    left_join(clin_df, by = "sample_id")
  
  # Prepare data for logistic regression
  outcome <- as.factor(merged_df[[variable_string]]) # Ensure the outcome is a factor
  predictors <- merged_df %>% 
    select(-all_of(variable_string), -sample_id) %>% 
    data.frame()
  
  # Initialize results data frame
  results <- data.frame(Protein = colnames(predictors), Coefficient = NA, PValue = NA, stringsAsFactors = FALSE)
  
  # Loop through each gene to perform logistic regression
  for (i in seq_along(colnames(predictors))) {
    npx <- predictors[, i]
    model <- clm(outcome ~ npx)
    summary_model <- summary(model)
    
    # Extract coefficient and p-value for the gene
    results$Coefficient[i] <- summary_model$coefficients["npx", "Estimate"]
    results$PValue[i] <- summary_model$coefficients["npx", "Pr(>|z|)"]
  }
  
  # Adjust p-values for multiple testing, for example, using Benjamini-Hochberg method
  results$AdjPValue <- p.adjust(results$PValue, method = "BH")
  
  # Return the results
  return(results)
}

# Plotting DEA results in volcano plot
volcano <- function(results) {
  
  significance_threshold <- -log10(0.05)
  
  results <- results |> 
    mutate(rank = rank(`AdjPValue`)) |> 
    mutate(label = ifelse(rank <= 10, as.character(Protein), NA))
  
  volcano_plot <- ggplot(results, aes(x = Estimate, y = -log10(`AdjPValue`))) +
    geom_point(aes(color = ifelse(-log10(`AdjPValue`) > significance_threshold, "Significant", "Not significant")), alpha = 0.8) +
    geom_hline(yintercept = significance_threshold, linetype = "dashed") +
    scale_color_npg() +
    ggrepel::geom_text_repel(aes(label = label), box.padding = 0.1, point.padding = 1, 
                             segment.color = 'grey50', size = 3, max.iter = 10000) +
    labs(color = "Significance") +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(volcano_plot)
}

volcano_logistic <- function(results) {
  
  significance_threshold <- -log10(0.05)
  
  results <- results |> 
    mutate(rank = rank(`AdjPValue`)) |> 
    mutate(label = ifelse(rank <= 10, as.character(Protein), NA))
  
  volcano_plot <- ggplot(results, aes(x = Coefficient, y = -log10(`AdjPValue`))) +
    geom_point(aes(color = ifelse(-log10(`AdjPValue`) > significance_threshold, "Significant", "Not significant")), alpha = 0.8) +
    geom_hline(yintercept = significance_threshold, linetype = "dashed") +
    scale_color_npg() +
    ggrepel::geom_text_repel(aes(label = label), box.padding = 0.1, point.padding = 1, 
                             segment.color = 'grey50', size = 3, max.iter = 10000) +
    labs(color = "Significance") +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(volcano_plot)
}

# Multinomial

multinomial_analysis <- function(df, clin_df, variable_string) {
  
  # Merge clinical and gene expression data based on sample ID
  clin_df <- clin_df |> 
    select(all_of(variable_string), sample_id)
  
  merged_df <- df |>
    left_join(clin_df, by = "sample_id")
  
  # Prepare data for logistic regression
  outcome <- as.factor(merged_df[[variable_string]]) # Ensure the outcome is a factor
  predictors <- merged_df %>% 
    select(-all_of(variable_string), -sample_id) %>% 
    data.frame()
  
  # Initialize results data frame
  results <- data.frame(Protein = colnames(predictors), 
                        Coefficient_1 = NA, PValue_1 = NA, 
                        Coefficient_2 = NA, PValue_2 = NA, 
                        AdjPValue_1 = NA, AdjPValue_2 = NA,
                        stringsAsFactors = FALSE)
  
  # Loop through each gene to perform logistic regression
  for (i in seq_along(colnames(predictors))) {
    npx <- predictors[, i]
    model <- vglm(outcome ~ npx, family = "multinomial")
    
    summary_model <- summary(model) |> 
      coef()
    
    # Extract coefficient and p-value for the gene
    results$Coefficient_1[i] <- summary_model[["npx:1", "Estimate"]]
    results$PValue_1[i] <- summary_model[["npx:1", "Pr(>|z|)"]]
    
    results$Coefficient_2[i] <- summary_model[["npx:2", "Estimate"]]
    results$PValue_2[i] <- summary_model[["npx:2", "Pr(>|z|)"]]
  }
  
  # Adjust p-values for multiple testing, for example, using Benjamini-Hochberg method
  results$AdjPValue_1 <- p.adjust(results$PValue_1, method = "BH")
  results$AdjPValue_2 <- p.adjust(results$PValue_2, method = "BH")
  
  # Return the results
  return(results)
}


volcano_multi <- function(results, option) {
  
  significance_threshold <- -log10(0.05)
  
  results_1 <- results |> 
    mutate(rank = rank(`AdjPValue_1`)) |> 
    mutate(label = ifelse(rank <= 10, as.character(Protein), NA))
  
  volcano_plot_1 <- ggplot(results_1, aes(x = Coefficient_1, y = -log10(`AdjPValue_1`))) +
    geom_point(aes(color = ifelse(-log10(`AdjPValue_1`) > significance_threshold, "Significant", "Not significant")), alpha = 0.8) +
    geom_hline(yintercept = significance_threshold, linetype = "dashed") +
    scale_color_npg() +
    ggrepel::geom_text_repel(aes(label = label), box.padding = 0.1, point.padding = 1, 
                             segment.color = 'grey50', size = 3, max.iter = 10000) +
    labs(color = "Significance") +
    theme_minimal() +
    theme(legend.position = "none")
  
  results_2 <- results |> 
    mutate(rank = rank(`AdjPValue_2`)) |> 
    mutate(label = ifelse(rank <= 10, as.character(Protein), NA))
  
  volcano_plot_2 <- ggplot(results_2, aes(x = Coefficient_2, y = -log10(`AdjPValue_2`))) +
    geom_point(aes(color = ifelse(-log10(`AdjPValue_2`) > significance_threshold, "Significant", "Not significant")), alpha = 0.8) +
    geom_hline(yintercept = significance_threshold, linetype = "dashed") +
    scale_color_npg() +
    ggrepel::geom_text_repel(aes(label = label), box.padding = 0.1, point.padding = 1, 
                             segment.color = 'grey50', size = 3, max.iter = 10000) +
    labs(color = "Significance") +
    theme_minimal() +
    theme(legend.position = "none")
  
  if (option == "1") {
    return(volcano_plot_1)
  }
  else {
    return(volcano_plot_2)
  }
  
}
