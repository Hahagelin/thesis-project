

prepper <- function(npx, clinical, targets) {

  # Ensuring factor format of appropriate clinical variables
  clinical$primary_blood_draw <- as_factor(clinical$primary_blood_draw)
  clinical$x_BC <- as_factor(clinical$x_BC)
  clinical$menopause_status <- as_factor(clinical$menopause_status)
  clinical$mht_status <- as_factor(clinical$mht_status)
  clinical$smoking_status <- as_factor(clinical$smoking_status)
  clinical$cancer_stage <- as_factor(clinical$cancer_stage)
  
  selected_columns <- c("sample_id", targets)
  
  # Combining selected clinical variables with the NPX
  prepped_npx <- npx |> 
    left_join(clinical |> select(all_of(selected_columns)),  by = "sample_id") |> 
    select(-sample_id) |> 
    filter(if_all(all_of(targets), ~ !is.na(.)))
  
  return(prepped_npx)
}

classifier <- function(skane_df, sthlm_df, target, engine, option = "binary") {
  
  # Convert target to formula
  target_formula <- as.formula(paste(target, "~ ."))
  
  # Set seeds and prepare data splits with target as strata
  set.seed(123)
  skane_split <- initial_split(skane_df, prop = 3/4, strata = !!rlang::sym(target), breaks = 4)
  skane_train <- training(skane_split)
  skane_test <- testing(skane_split)
  
  set.seed(234)
  skane_folds <- vfold_cv(skane_train, v = 3, strata = !!rlang::sym(target))
  
  set.seed(456)
  sthlm_split <- initial_split(sthlm_df, prop = 3/4, strata = !!rlang::sym(target), breaks = 4)
  sthlm_train <- training(sthlm_split)
  sthlm_test <- testing(sthlm_split)
  
  set.seed(780)
  sthlm_folds <- vfold_cv(sthlm_train, v = 3, strata = !!rlang::sym(target))
  
  # Define recipes using the dynamic formula
  skane_recipe <- recipe(target_formula, data = skane_train) |> 
    step_zv(all_predictors()) |> 
    step_corr(threshold = 0.4, method = "spearman")
  
  sthlm_recipe <- recipe(target_formula, data = sthlm_train) |> 
    step_zv(all_predictors()) |> 
    step_corr(threshold = 0.4, method = "spearman")
  
  if (engine == "xgboost" || engine == "lightgbm") {
    parsnip_model <- boost_tree(
      trees = 500,           # Number of trees
      tree_depth = 3,        # Maximum depth of a tree
      min_n = 20,            # Minimum observations in a node
      loss_reduction = 0,  # Minimum loss reduction required for a split
      sample_size = 0.9,     # Proportion of samples used per tree
      mtry = 20,              # Number of variables considered for each split
      learn_rate = 0.05,
      stop_iter = 10
    ) |>  
      set_engine(engine) |> 
      set_mode("classification") 
  }
  
  if (engine == "randomForest") {
    parsnip_model <- rand_forest(
      trees = 500,
      mtry = NULL,
      min_n = 20
    ) |>  
      set_engine(engine) |> 
      set_mode("classification") 
  }
  
  if (engine == "glmnet" && option == "binary") {
    
    parsnip_model <- logistic_reg(
      penalty = 0,
      mixture = 0
    ) |>  
      set_engine(engine) |> 
      set_mode("classification") 
  } 
  
  if (engine == "glmnet" || engine == "nnet" && option == "multinomial") {
    
    parsnip_model <- multinom_reg(
      penalty = 0,
      mixture = 0
    ) |>  
      set_engine(engine) |> 
      set_mode("classification") 
  } 
  
  # Create workflows
  skane_wf <- workflow() |> 
    add_recipe(skane_recipe) |>  
    add_model(parsnip_model)
  
  sthlm_wf <- workflow() |> 
    add_recipe(sthlm_recipe) |>  
    add_model(parsnip_model)
  
  # Fit models using resampling
  skane_fit <- fit_resamples(
    skane_wf,
    resamples = skane_folds,
    metrics = metric_set(accuracy, roc_auc, sens, spec),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )
  
  sthlm_fit <- fit_resamples(
    sthlm_wf,
    resamples = sthlm_folds,
    metrics = metric_set(accuracy, roc_auc, sens, spec),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )

  # Final fitting on full training set and evaluation on test set
  skane_final <- last_fit(
    skane_wf,
    split = skane_split,
    metrics = metric_set(accuracy, roc_auc, sens, spec)
  )

  sthlm_final <- last_fit(
    sthlm_wf,
    split = sthlm_split,
    metrics = metric_set(accuracy, roc_auc, sens, spec)
  )
  
  # Extracting models for VIP
  skane_model_for_vip <- extract_fit_parsnip(skane_final)
  sthlm_model_for_vip <- extract_fit_parsnip(sthlm_final)
  
  # Creating performance overview element in the output
  performance_summary <- list(
    Skane = collect_metrics(skane_fit),
    Stockholm = collect_metrics(sthlm_fit),
    Skane_Final_Test = collect_metrics(skane_final),
    Stockholm_Final_Test = collect_metrics(sthlm_final)
  )
  
  prediction_summary <- list(
    Skane = collect_predictions(skane_fit),
    Stockholm = collect_predictions(sthlm_fit),
    Skane_Final_Test = collect_predictions(skane_final),
    Stockholm_Final_Test = collect_predictions(sthlm_final)
  )
    
  
  ## Tuning
  
  # XGBoost
  if (engine == "xgboost" | engine == "lightgbm") {
    model_spec_tune <- boost_tree(
      trees = 500,                   
      tree_depth = tune(),               
      min_n = tune(),                     
      loss_reduction = tune(),            
      mtry = tune(),                     
      learn_rate = 0.01                
    ) |>
      set_engine(engine) |>           
      set_mode("classification")
    
    grid <- grid_latin_hypercube(
      tree_depth(),
      min_n(),
      loss_reduction(),
      mtry = finalize(mtry(), skane_train),  # Limits mtry() to range from 1 to ncol(train_data)
      size = 20 
    )
  }
  
  # Random forest
  if (engine == "randomForest") {
    model_spec_tune <- rand_forest(
      trees = 500,           
      mtry = tune(),             
      min_n = tune()            
    ) |> 
      set_engine("randomForest") |>
      set_mode("classification")
    
    
    grid <- grid_latin_hypercube(
      mtry = finalize(mtry(), skane_train),  
      min_n(),
      size = 20
    )
  }
  
  # Glmnet
  if (engine == "glmnet") {
    model_spec_tune <- logistic_reg(  # Use multinom_reg for multinomial cases
      penalty = tune(),
      mixture = tune()
    ) |> 
      set_engine("glmnet") |>
      set_mode("classification")
    
    grid <- grid_latin_hypercube(
      penalty(),
      mixture(),
      size = 20
    )
  }
  
  # Workflow
  tune_wf <- workflow() |> 
    add_formula(x_BC ~.) %>%
    add_model(model_spec_tune)
  
  # Running the tuning
  set.seed(234)  # Ensure reproducibility
  skane_tuned <- tune_grid(
    tune_wf,
    resamples = skane_folds,
    grid = grid,
    control = control_grid(save_pred = TRUE, verbose = TRUE)
  )
  
  sthlm_tuned <- tune_grid(
    tune_wf,
    resamples = sthlm_folds,
    grid = grid,
    control = control_grid(save_pred = TRUE, verbose = TRUE)
  )
  
  # Retrieving parameters
  best_auc_skane <- select_best(skane_tuned, metric = "roc_auc")
  best_params_skane <- collect_metrics(skane_tuned) |> 
    filter(.metric == "roc_auc")

  
  best_auc_sthlm <- select_best(sthlm_tuned, metric = "roc_auc")
  best_params_sthlm <- collect_metrics(sthlm_tuned) |> 
    filter(.metric == "roc_auc")

  
  # Return a list containing all relevant model objects
  return(list(
    Skane = list(
      Fit = skane_fit,
      Final = skane_final,
      VIP_Model = skane_model_for_vip,
      Tuned_AUC = best_auc_skane,
      Best_Params = best_params_skane
    ),
    Stockholm = list(
      Fit = sthlm_fit,
      Final = sthlm_final,
      VIP_Model = sthlm_model_for_vip,
      Tuned_AUC = best_auc_sthlm,
      Best_Params = best_params_sthlm
    ),
    Performance = performance_summary,
    Predictions = prediction_summary
  ))
}

regressor <- function(skane_df, sthlm_df, target, engine) {
  
  # Convert target to formula
  target_formula <- as.formula(paste(target, "~ ."))
  
  # Set seeds and prepare data splits with target as strata
  set.seed(123)
  skane_split <- initial_split(skane_df, prop = 3/4, strata = !!rlang::sym(target), breaks = 4)
  skane_train <- training(skane_split)
  skane_test <- testing(skane_split)
  
  set.seed(234)
  skane_folds <- vfold_cv(skane_train, v = 5, strata = !!rlang::sym(target))
  
  set.seed(456)
  sthlm_split <- initial_split(sthlm_df, prop = 3/4, strata = !!rlang::sym(target), breaks = 4)
  sthlm_train <- training(sthlm_split)
  sthlm_test <- testing(sthlm_split)
  
  set.seed(780)
  sthlm_folds <- vfold_cv(sthlm_train, v = 5, strata = !!rlang::sym(target))
  
  # Define recipes using the dynamic formula
  skane_recipe <- recipe(target_formula, data = skane_train) |> 
    step_zv(all_predictors()) |> 
    step_corr(threshold = 0.4, method = "spearman")
  
  sthlm_recipe <- recipe(target_formula, data = sthlm_train) |> 
    step_zv(all_predictors()) |> 
    step_corr(threshold = 0.4, method = "spearman")
  
  if (engine == "xgboost" || engine == "lightgbm") {
    parsnip_model <- boost_tree(
      trees = 500,           # Number of trees
      tree_depth = 10,        # Maximum depth of a tree
      min_n = 20,            # Minimum observations in a node
      loss_reduction = 5,    # Minimum loss reduction required for a split
      sample_size = 0.9,     # Proportion of samples used per tree
      mtry = 20,             # Number of variables considered for each split
      learn_rate = 0.05,
      stop_iter = 10
    ) |>  
      set_engine(engine) |> 
      set_mode("regression") 
  }
  
  # Create workflows
  skane_wf <- workflow() |> 
    add_recipe(skane_recipe) |>  
    add_model(parsnip_model)
  
  sthlm_wf <- workflow() |> 
    add_recipe(sthlm_recipe) |>  
    add_model(parsnip_model)
  
  # Fit models using resampling
  skane_fit <- fit_resamples(
    skane_wf,
    resamples = skane_folds,
    metrics = metric_set(rmse, rsq, mae),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )
  
  sthlm_fit <- fit_resamples(
    sthlm_wf,
    resamples = sthlm_folds,
    metrics = metric_set(rmse, rsq, mae),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )
  
  # Final fitting on full training set and evaluation on test set
  skane_final <- last_fit(
    skane_wf,
    split = skane_split,
    metrics = metric_set(rmse, rsq, mae)
  )
  
  sthlm_final <- last_fit(
    sthlm_wf,
    split = sthlm_split,
    metrics = metric_set(rmse, rsq, mae)
  )
  
  # Extracting models for VIP
  skane_model_for_vip <- extract_fit_parsnip(skane_final)
  sthlm_model_for_vip <- extract_fit_parsnip(sthlm_final)
  
  # Creating performance overview element in the output
  performance_summary <- list(
    Skane = collect_metrics(skane_fit),
    Stockholm = collect_metrics(sthlm_fit),
    Skane_Final_Test = collect_metrics(skane_final),
    Stockholm_Final_Test = collect_metrics(sthlm_final)
  )
  
  # Return a list containing all relevant model objects
  return(list(
    Skane = list(
      Fit = skane_fit,
      Final = skane_final,
      VIP_Model = skane_model_for_vip
    ),
    Stockholm = list(
      Fit = sthlm_fit,
      Final = sthlm_final,
      VIP_Model = sthlm_model_for_vip
    ),
    Performance = performance_summary
  ))
}
