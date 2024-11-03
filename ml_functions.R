# Function that balances the groups (case-control).
balance_groups <- function(join_data,
                           variable,
                           case = 1,
                           seed = 123) {
  
  Variable <- rlang::sym(variable)
  set.seed(seed)
  
  case_data <- join_data |> filter(!!Variable == case)
  control_data <- join_data |> filter(!!Variable != case)
  
  case_sample_num <- nrow(case_data)
  
  group <- case_data
  
  control_data <- control_data |>
    dplyr::filter(!!Variable != case) |>
    dplyr::sample_n(size = case_sample_num, replace = TRUE)
  
  group <- rbind(group, control_data)
  
  return(group)
}


# Function that generates subtitles for the variable importance plots.
generate_subtitle <- function(features,
                              accuracy,
                              sensitivity,
                              specificity,
                              auc,
                              subtitle = c("accuracy",
                                           "sensitivity",
                                           "specificity",
                                           "auc",
                                           "features",
                                           "top-features")) {
  
  subtitle_parts <- c()
  
  if ("accuracy" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('accuracy = ', round(accuracy, 2), '    '))
  }
  
  if ("sensitivity" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('sensitivity = ', round(sensitivity, 2), '    '))
  }
  
  if ("specificity" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('specificity = ', round(specificity, 2), '    '))
  }
  
  if ("auc" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('AUC = ', round(auc, 2), '    '))
  }
  
  if (length(subtitle_parts) > 0) {
    subtitle_parts <- c(subtitle_parts, '\n')
  }
  
  if ("features" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('Features = ', nrow(features), '    '))
  }
  
  if ("top-features" %in% subtitle) {
    subtitle_parts <- c(subtitle_parts, paste0('top-features = ',
                                               nrow(features |> dplyr::filter(Scaled_Importance >= 50)),
                                               '    '))
  }
  
  subtitle <- paste(subtitle_parts, collapse = '')
  
  return(subtitle)
}


# Function that creates boxplots of the top n features.
plot_protein_boxplot <- function(join_data,
                                 variable,
                                 features,
                                 points = TRUE,
                                 xaxis_names = TRUE,
                                 palette = disease_palette,
                                 nfeatures = 9) {
  
  Variable <- rlang::sym(variable)
  
  top_features <- features |>
    arrange(desc(Scaled_Importance)) |>
    select(Variable) |>
    mutate(across(everything(), as.character)) |>
    head(nfeatures)
  proteins <- top_features[["Variable"]]
  
  long_data <- join_data |>
    select(!!Variable, any_of(proteins)) |>
    pivot_longer(cols = !any_of(c(variable)),
                 names_to = "Protein",
                 values_to = "NPX")
  
  long_data$Protein <- factor(long_data$Protein, levels = proteins, labels = proteins)
  long_data[[variable]] <- as.factor(long_data[[variable]])
  
  # Create boxplot
  boxplot <- long_data |>
    ggplot(aes(x = !!Variable, y = NPX)) +
    geom_boxplot() +
    geom_boxplot(data = filter(long_data, !!Variable == 1),
                 fill = palette[variable],
                 show.legend = FALSE)
  
  if (isTRUE(points)) {
    boxplot <- boxplot +
      geom_point(data = filter(long_data, !!Variable != 1),
                 position = position_jitter(width = 0.1),
                 color = 'grey',
                 alpha = 0.3) +
      geom_point(data = filter(long_data, !!Variable == 1),
                 aes(fill = !!Variable),
                 position = position_jitter(width = 0.1),
                 color = palette[variable],
                 alpha = 0.5,
                 show.legend = FALSE)
  }
  
  boxplot_panel <- boxplot +
    theme(legend.position = 'none') +
    xlab('') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90))
  
  if (isFALSE(xaxis_names)) {
    boxplot_panel <- boxplot_panel +
      theme(axis.text.x = element_blank())
  }
  
  boxplot_panel <- boxplot_panel +
    facet_wrap(~ Protein, scale="free_y")
  
  return(boxplot_panel)
}


# Function to extract model results: Features & AUC.
extract_model_res <- function(model, diseases) {
  features_all <- lapply(diseases, function(disease) {
    features <- model[[disease]]$features |> mutate(Disease = disease)
  })
  auc_all <- lapply(diseases, function(disease) {
    auc <- model[[disease]]$auc |> as_tibble() |> rename(AUC = value) |> mutate(Disease = disease)
  })
  return(list("features" = bind_rows(features_all),
              "auc" = bind_rows(auc_all)))
}


## LASSO
lasso <- function(join_data,
                  variable,
                  seed = 123,
                  cv_sets = 5,
                  cor_threshold = 0.9,
                  grid_size = 10,
                  subtitle = c("accuracy",
                               "sensitivity",
                               "specificity",
                               "auc",
                               "features",
                               "top-features"),
                  palette = disease_palette,
                  binary_cols = NULL) {
  Variable <- rlang::sym(variable)
  
  # Prepare sets
  set.seed(seed)
  data_split <- initial_split(join_data, prop = 0.8, strata = !!Variable)
  
  train_set <- training(data_split) |> 
    mutate(!!Variable := as.factor(!!Variable))
  test_set <- testing(data_split) |> 
    mutate(!!Variable := as.factor(!!Variable))
  
  train_set <- balance_groups(train_set,
                              variable,
                              case = 1,
                              seed)
  
  train_folds <- vfold_cv(train_set, v = cv_sets, strata = !!Variable)
  
  print("Starting training...")
  
  # Train model - hyperparameter optimization
  formula <- stats::as.formula(paste(variable, "~ ."))
  
  recipe <- recipe(formula, data = train_set) |> 
    update_role(DAid, new_role = "id") |> 
    step_dummy(all_nominal(), -all_outcomes(), -has_role("id")) |> 
    step_zv(all_predictors(), -has_role("id")) |> 
    step_normalize(all_numeric(), -all_outcomes(), -has_role("id"), -all_of(binary_cols)) |> 
    step_corr(all_numeric(), -all_outcomes(), -has_role("id"), -any_of(binary_cols), threshold = cor_threshold)
  
  model <- logistic_reg(penalty = tune(), mixture = 1) |>
    set_engine("glmnet")
  
  workflow <- workflow() |> 
    add_recipe(recipe) |> 
    add_model(model)
  
  grid <- workflow |> 
    extract_parameter_set_dials() |>
    grid_space_filling(size = grid_size, type = "latin_hypercube")
  
  roc_res <- metric_set(roc_auc)
  
  ctrl <- control_grid(save_pred = TRUE, parallel_over = "everything", verbose = TRUE)
  tune <- workflow |>
    tune_grid(train_folds, grid = grid, control = ctrl, metrics = roc_res)
  
  print("Selecting best model...")
  
  # Select best model - Final fit
  best <- tune |>
    select_best(metric = "roc_auc") |>
    select(-.config)
  
  final_wf <- finalize_workflow(workflow, best)
  
  final <- final_wf |>
    fit(train_set)
  
  print("Evaluating model...")
  
  # Evaluate model
  splits <- make_splits(train_set, test_set)
  
  preds <- last_fit(final_wf, splits, metrics = metric_set(roc_auc))
  
  res <- predict(final, new_data = test_set)
  
  res <- bind_cols(res, test_set |> select(!!Variable))
  
  accuracy <- res |> accuracy(!!Variable, .pred_class)
  sensitivity <- res |> sensitivity(!!Variable, .pred_class, event_level = "second")
  specificity <- res |> specificity(!!Variable, .pred_class, event_level = "second")
  auc <- preds |> collect_metrics()
  cm <- res |> conf_mat(!!Variable, .pred_class)
  roc <- preds |>
    collect_predictions(summarize = F) |>
    roc_curve(truth = !!Variable, .pred_0) |>
    ggplot(aes(x = 1 - specificity, y = sensitivity)) +
    geom_path(colour = palette[variable], linewidth = 2) +
    geom_abline(lty = 3) +
    coord_equal() +
    theme_hpa()
  
  # Feature importance
  features <- final |>
    extract_fit_parsnip() |>
    vi() |>
    mutate(Importance = abs(Importance),
           Variable = forcats::fct_reorder(Variable, Importance)) |>
    arrange(desc(Importance)) |>
    mutate(Scaled_Importance = scales::rescale(Importance, to = c(0, 100))) |>
    filter(Scaled_Importance > 0)
  
  subtitle_text <- generate_subtitle(features, 
                                     accuracy$.estimate, 
                                     sensitivity$.estimate, 
                                     specificity$.estimate, 
                                     auc$.estimate, 
                                     subtitle)
  
  var_imp_plot <- features |>
    ggplot(aes(x = Scaled_Importance, y = Variable)) +
    geom_col(aes(fill = ifelse(Scaled_Importance > 50, variable, NA))) +
    labs(y = NULL) +
    scale_x_continuous(breaks = c(0, 100), expand = c(0, 0)) +  # Keep x-axis tick labels at 0 and 100
    scale_fill_manual(values = palette, na.value = "grey50") +
    ggtitle(label = paste0(variable, ''),
            subtitle = subtitle_text) +
    xlab('Importance') +
    ylab('Features') +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  return(list("final_wf" = final_wf,
              "res" = res,
              "accuracy" = accuracy$.estimate,
              "sensitivity" = sensitivity$.estimate,
              "specificity" = specificity$.estimate,
              "auc" = auc$.estimate,
              "confusion_matrix" = cm,
              "roc_curve" = roc,
              "features" = features,
              "var_imp_plot" = var_imp_plot))
}


# Elastic Net
elasticnet <- function(join_data,
                       variable,
                       seed = 123,
                       cv_sets = 5,
                       cor_threshold = 0.9,
                       grid_size = 10,
                       subtitle = c("accuracy",
                                    "sensitivity",
                                    "specificity",
                                    "auc",
                                    "features",
                                    "top-features"),
                       palette = disease_palette,
                       binary_cols = NULL) {
  Variable <- rlang::sym(variable)
  
  # Prepare sets
  set.seed(seed)
  data_split <- initial_split(join_data, prop = 0.8, strata = !!Variable)
  
  train_set <- training(data_split) |> 
    mutate(!!Variable := as.factor(!!Variable))
  test_set <- testing(data_split) |> 
    mutate(!!Variable := as.factor(!!Variable))
  
  train_set <- balance_groups(train_set,
                              variable,
                              case = 1,
                              seed)
  
  train_folds <- vfold_cv(train_set, v = cv_sets, strata = !!Variable)
  
  print("Starting training...")
  
  # Train model - hyperparameter optimization
  formula <- stats::as.formula(paste(variable, "~ ."))
  
  recipe <- recipe(formula, data = train_set) |> 
    update_role(DAid, new_role = "id") |> 
    step_dummy(all_nominal(), -all_outcomes(), -has_role("id")) |> 
    step_zv(all_predictors(), -has_role("id")) |> 
    step_normalize(all_numeric(), -all_outcomes(), -has_role("id"), -all_of(binary_cols)) |> 
    step_corr(all_numeric(), -all_outcomes(), -has_role("id"), -all_of(binary_cols), threshold = cor_threshold)
  
  model <- logistic_reg(penalty = tune(), mixture = tune()) |>
    set_engine("glmnet")
  
  workflow <- workflow() |> 
    add_recipe(recipe) |> 
    add_model(model)
  
  grid <- workflow |> 
    extract_parameter_set_dials() |>
    grid_space_filling(size = grid_size, type = "latin_hypercube")
  
  roc_res <- metric_set(roc_auc)
  
  ctrl <- control_grid(save_pred = TRUE, parallel_over = "everything", verbose = TRUE)
  tune <- workflow |>
    tune_grid(train_folds, grid = grid, control = ctrl, metrics = roc_res)
  
  print("Selecting best model...")
  
  # Select best model - Final fit
  best <- tune |>
    select_best(metric = "roc_auc") |>
    select(-.config)
  
  final_wf <- finalize_workflow(workflow, best)
  
  final <- final_wf |>
    fit(train_set)
  
  print("Evaluating model...")
  
  # Evaluate model
  splits <- make_splits(train_set, test_set)
  
  preds <- last_fit(final_wf, splits, metrics = metric_set(roc_auc))
  
  res <- predict(final, new_data = test_set)
  
  res <- bind_cols(res, test_set |> select(!!Variable))
  
  accuracy <- res |> accuracy(!!Variable, .pred_class)
  sensitivity <- res |> sensitivity(!!Variable, .pred_class, event_level = "second")
  specificity <- res |> specificity(!!Variable, .pred_class, event_level = "second")
  auc <- preds |> collect_metrics()
  cm <- res |> conf_mat(!!Variable, .pred_class)
  roc <- preds |>
    collect_predictions(summarize = F) |>
    roc_curve(truth = !!Variable, .pred_0) |>
    ggplot(aes(x = 1 - specificity, y = sensitivity)) +
    geom_path(colour = palette[variable], linewidth = 2) +
    geom_abline(lty = 3) +
    coord_equal() +
    theme_hpa()
  
  # Feature importance
  features <- final |>
    extract_fit_parsnip() |>
    vi() |>
    mutate(Importance = abs(Importance),
           Variable = forcats::fct_reorder(Variable, Importance)) |>
    arrange(desc(Importance)) |>
    mutate(Scaled_Importance = scales::rescale(Importance, to = c(0, 100))) |>
    filter(Scaled_Importance > 0)
  
  subtitle_text <- generate_subtitle(features, 
                                     accuracy$.estimate, 
                                     sensitivity$.estimate, 
                                     specificity$.estimate, 
                                     auc$.estimate, 
                                     subtitle)
  
  var_imp_plot <- features |>
    ggplot(aes(x = Scaled_Importance, y = Variable)) +
    geom_col(aes(fill = ifelse(Scaled_Importance > 50, variable, NA))) +
    labs(y = NULL) +
    scale_x_continuous(breaks = c(0, 100), expand = c(0, 0)) +  # Keep x-axis tick labels at 0 and 100
    scale_fill_manual(values = palette, na.value = "grey50") +
    ggtitle(label = paste0(variable, ''),
            subtitle = subtitle_text) +
    xlab('Importance') +
    ylab('Features') +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  return(list("final_wf" = final_wf,
              "res" = res,
              "accuracy" = accuracy$.estimate,
              "sensitivity" = sensitivity$.estimate,
              "specificity" = specificity$.estimate,
              "auc" = auc$.estimate,
              "confusion_matrix" = cm,
              "roc_curve" = roc,
              "features" = features,
              "var_imp_plot" = var_imp_plot))
}


# Random Forest
rf <- function(join_data,
               variable,
               seed = 123,
               cv_sets = 5,
               cor_threshold = 0.9,
               grid_size = 10,
               subtitle = c("accuracy",
                            "sensitivity",
                            "specificity",
                            "auc",
                            "features",
                            "top-features"),
               palette = disease_palette,
               binary_cols = NULL) {
  Variable <- rlang::sym(variable)
  
  # Prepare sets
  set.seed(seed)
  data_split <- initial_split(join_data, prop = 0.8, strata = !!Variable)
  
  train_set <- training(data_split) |> 
    mutate(!!Variable := as.factor(!!Variable))
  test_set <- testing(data_split) |> 
    mutate(!!Variable := as.factor(!!Variable))
  
  train_set <- balance_groups(train_set,
                              variable,
                              case = 1,
                              seed)
  
  train_folds <- vfold_cv(train_set, v = cv_sets, strata = !!Variable)
  
  print("Starting training...")
  
  # Train model - hyperparameter optimization
  formula <- stats::as.formula(paste(variable, "~ ."))
  
  recipe <- recipe(formula, data = train_set) |> 
    update_role(DAid, new_role = "id") |> 
    step_dummy(all_nominal(), -all_outcomes(), -has_role("id")) |> 
    step_zv(all_predictors(), -has_role("id")) |> 
    step_normalize(all_numeric(), -all_outcomes(), -has_role("id"), -all_of(binary_cols)) |> 
    step_corr(all_numeric(), -all_outcomes(), -has_role("id"), -all_of(binary_cols), threshold = cor_threshold)
  
  model <- rand_forest(trees = 1000,
                       min_n = tune(),
                       mtry = tune()) |>
    parsnip::set_mode("classification") |>
    parsnip::set_engine("ranger", importance = "permutation")
  
  prepped_recipe <- prep(recipe, training = train_set)  # Prep the recipe
  baked_data <- bake(prepped_recipe, new_data = train_set)
  remaining_predictors <- colnames(baked_data)[!colnames(baked_data) %in% c(variable, "DAid")]
  n_remaining_predictors <- length(remaining_predictors)
  
  workflow <- workflow() |> 
    add_recipe(recipe) |> 
    add_model(model)
  
  grid <- grid_space_filling(min_n(), 
                             mtry(range = c(floor(sqrt(n_remaining_predictors)),
                                            (floor(n_remaining_predictors/3)))),
                             size = grid_size,
                             type = "latin_hypercube")
  
  roc_res <- metric_set(roc_auc)
  
  ctrl <- control_grid(save_pred = TRUE, parallel_over = "everything", verbose = TRUE)
  tune <- workflow |>
    tune_grid(train_folds, grid = grid, control = ctrl, metrics = roc_res)
  
  print("Selecting best model...")
  
  # Select best model - Final fit
  best <- tune |>
    select_best(metric = "roc_auc") |>
    select(-.config)
  
  final_wf <- finalize_workflow(workflow, best)
  
  final <- final_wf |>
    fit(train_set)
  
  print("Evaluating model...")
  
  # Evaluate model
  splits <- make_splits(train_set, test_set)
  
  preds <- last_fit(final_wf, splits, metrics = metric_set(roc_auc))
  
  res <- predict(final, new_data = test_set)
  
  res <- bind_cols(res, test_set |> select(!!Variable))
  
  accuracy <- res |> accuracy(!!Variable, .pred_class)
  sensitivity <- res |> sensitivity(!!Variable, .pred_class, event_level = "second")
  specificity <- res |> specificity(!!Variable, .pred_class, event_level = "second")
  auc <- preds |> collect_metrics()
  cm <- res |> conf_mat(!!Variable, .pred_class)
  roc <- preds |>
    collect_predictions(summarize = F) |>
    roc_curve(truth = !!Variable, .pred_0) |>
    ggplot(aes(x = 1 - specificity, y = sensitivity)) +
    geom_path(colour = palette[variable], linewidth = 2) +
    geom_abline(lty = 3) +
    coord_equal() +
    theme_hpa()
  
  # Feature importance
  features <- final |>
    extract_fit_parsnip() |>
    vi() |>
    mutate(Importance = abs(Importance),
           Variable = forcats::fct_reorder(Variable, Importance)) |>
    arrange(desc(Importance)) |>
    mutate(Scaled_Importance = scales::rescale(Importance, to = c(0, 100))) |>
    filter(Scaled_Importance > 0)
  
  subtitle_text <- generate_subtitle(features, 
                                     accuracy$.estimate, 
                                     sensitivity$.estimate, 
                                     specificity$.estimate, 
                                     auc$.estimate, 
                                     subtitle)
  
  var_imp_plot <- features |>
    ggplot(aes(x = Scaled_Importance, y = Variable)) +
    geom_col(aes(fill = ifelse(Scaled_Importance > 50, variable, NA))) +
    labs(y = NULL) +
    scale_x_continuous(breaks = c(0, 100), expand = c(0, 0)) +  # Keep x-axis tick labels at 0 and 100
    scale_fill_manual(values = palette, na.value = "grey50") +
    ggtitle(label = paste0(variable, ''),
            subtitle = subtitle_text) +
    xlab('Importance') +
    ylab('Features') +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  return(list("final_wf" = final_wf,
              "res" = res,
              "accuracy" = accuracy$.estimate,
              "sensitivity" = sensitivity$.estimate,
              "specificity" = specificity$.estimate,
              "auc" = auc$.estimate,
              "confusion_matrix" = cm,
              "roc_curve" = roc,
              "features" = features,
              "var_imp_plot" = var_imp_plot))
}


## LASSO Multiclass
lasso_multi <- function(join_data,
                        variable,
                        seed = 123,
                        cv_sets = 5,
                        cor_threshold = 0.9,
                        grid_size = 10,
                        subtitle = c("features", "top-features")) {
  Variable <- rlang::sym(variable)
  
  # Prepare sets
  set.seed(seed)
  data_split <- initial_split(join_data, prop = 0.8, strata = !!Variable)
  
  train_set <- training(data_split) |> 
    mutate(!!Variable := as.factor(!!Variable))
  test_set <- testing(data_split) |> 
    mutate(!!Variable := as.factor(!!Variable))
  
  train_folds <- vfold_cv(train_set, v = cv_sets, strata = !!Variable)
  
  print("Starting training...")
  
  # Train model - hyperparameter optimization
  formula <- stats::as.formula(paste(variable, "~ ."))
  
  recipe <- recipe(formula, data = train_set) |> 
    update_role(DAid, new_role = "id") |> 
    step_dummy(all_nominal(), -all_outcomes(), -has_role("id")) |> 
    step_zv(all_predictors(), -has_role("id")) |> 
    step_normalize(all_numeric(), -all_outcomes(), -has_role("id")) |> 
    step_corr(all_numeric(), -all_outcomes(), -has_role("id"), threshold = cor_threshold)
  
  model <- multinom_reg(penalty = tune(), mixture = 1) |>
    set_engine("glmnet")
  
  workflow <- workflow() |> 
    add_recipe(recipe) |> 
    add_model(model)
  
  grid <- workflow |> 
    extract_parameter_set_dials() |>
    grid_space_filling(size = grid_size, type = "latin_hypercube")
  
  roc_res <- metric_set(roc_auc)
  
  ctrl <- control_grid(save_pred = TRUE, parallel_over = "everything", verbose = TRUE)
  tune <- workflow |>
    tune_grid(train_folds, grid = grid, control = ctrl, metrics = roc_res)
  
  print("Selecting best model...")
  
  # Select best model - Final fit
  best <- tune |>
    select_best(metric = "roc_auc") |>
    select(-.config)
  
  final_wf <- finalize_workflow(workflow, best)
  
  final <- final_wf |>
    fit(train_set)
  
  print("Evaluating model...")
  
  # Evaluate model
  splits <- make_splits(train_set, test_set)
  
  preds <- last_fit(final_wf, splits, metrics = metric_set(roc_auc))
  
  class_predictions <- predict(final, new_data = test_set, type = "class")
  prob_predictions <- predict(final, new_data = test_set, type = "prob")
  
  res <- bind_cols(test_set |> select(!!Variable), class_predictions, prob_predictions)
  
  cm <- res |> conf_mat(!!Variable, .pred_class)

  pred_cols <- grep("^\\.pred_", names(res |> select(-.pred_class)), value = TRUE)

  roc_data <- roc_curve(res, truth = !!Variable, !!!rlang::syms(pred_cols))
  roc <- autoplot(roc_data)
  
  features <- final |>
    extract_fit_parsnip() |>
    vi() |>
    mutate(Importance = abs(Importance),
           Variable = forcats::fct_reorder(Variable, Importance)) |>
    arrange(desc(Importance)) |>
    mutate(Scaled_Importance = scales::rescale(Importance, to = c(0, 100))) |>
    filter(Scaled_Importance > 0)
  
  subtitle_text <- generate_subtitle(features = features, subtitle = subtitle)
  
  var_imp_plot <- features |>
    ggplot(aes(x = Scaled_Importance, y = Variable)) +
    geom_col(aes(fill = ifelse(Scaled_Importance > 50, variable, NA))) +
    labs(y = NULL) +
    scale_x_continuous(breaks = c(0, 100), expand = c(0, 0)) +  # Keep x-axis tick labels at 0 and 100
    ggtitle(label = paste0(variable, " - Multiclassification"),
            subtitle = subtitle_text) +
    xlab('Importance') +
    ylab('Features') +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  # ROC
  final_predictions <- prob_predictions |> 
    mutate(DAid = test_set$DAid) |> 
    relocate(DAid)
  
  dat <- test_set |> 
    select(DAid, !!Variable) |> 
    mutate(value = 1) |> 
    spread(!!Variable, value, fill= 0) 
  
  true_dat <- dat |> 
    set_names(paste(names(dat), "_true", sep = "")) |> 
    rename(DAid = `DAid_true`)
  
  dat_prob <- final_predictions |> 
    rename_all(~stringr::str_replace_all(.,".pred_",""))
  
  prob_data <- dat_prob |> 
    set_names(paste(names(dat_prob), "_pred_glmnet", sep = ""))|> 
    rename(DAid = DAid_pred_glmnet)
  
  final_df <- 
    true_dat |> 
    left_join(prob_data, by = "DAid") |>  
    select(-DAid) |> 
    as.data.frame()
  
  auc <- multi_roc(final_df, force_diag=T)
  auc <- tibble(!!Variable := names(auc$AUC$glmnet), AUC = unlist(auc$AUC$glmnet))
  
  return(list("final_wf" = final_wf,
              "res" = res,
              "auc" = auc,
              "confusion_matrix" = cm,
              "roc_curve" = roc,
              "features" = features,
              "var_imp_plot" = var_imp_plot))
}