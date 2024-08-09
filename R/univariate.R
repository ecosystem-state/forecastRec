#' Forecast recruitment deviations as a function of covariates using GLMs and GAMs
#'
#' \code{univariate_forecast} vector of responses and dataframe of predictors and automates
#'
#' @param response A data frame of responses for the modeling (values to be forecast) containing
#' a "time" column and "dev" column
#' @param predictors A data frame of predictors used for forecasting recruitment
#' @param model_type The type of model used to link predictors to forecasted recruitment. Can
#' be "lm", "gam",
#' @param n_forecast How many years to use as a holdout / test set
#' @param n_years_ahead How many years ahead to forecast (defaults to 1)
#' 1:n_vars variables, and then results are combined and sorted to remove duplicates
#' @param max_vars The maximum number of variables to include as predictors; defaults to 3
#' @import tidymodels
#' @importFrom parsnip set_mode set_engine linear_reg gen_additive_mod fit
#' @importFrom rsample rolling_origin training testing
#' @importFrom recipes recipe
#' @importFrom purrr map2 pluck map_dbl map_df
#' @importFrom workflows add_model add_recipe workflow extract_fit_parsnip
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr left_join bind_cols
#' @importFrom broom tidy
#' @importFrom ggeffects ggpredict
#' @importFrom stats lm as.formula predict cor model.matrix na.pass
#' @importFrom mgcv gam
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom tibble tibble
#'
#' @export
#'
#' @return a list containing predictions, with elements
#'
#' * `pred`: the predictions
#' * `vars`: the variable values used to fit the models
#' * `coefs`: coefficient estimates from each year:iteration
#' @examples
#' response <- data.frame(time = 1:40, dev = rnorm(40))
#'
#' predictors <- matrix(rnorm(400), ncol = 10)
#' #colnames(predictors) = paste0("X",1:ncol(predictors))
#' predictors <- as.data.frame(predictors)
#' predictors$time <- 1:40
#' lm_example <- univariate_forecast(response,
#'      predictors,
#'      model_type = "lm",
#'      n_forecast = 10,
#'      n_years_ahead = 1,
#'      max_vars = 3)
#'
univariate_forecast = function(response,
                               predictors,
                               model_type,
                               n_forecast = 10,
                               n_years_ahead = 1,
                               max_vars = 3) {

  if (!is.data.frame(predictors) || !is.data.frame(response)) {
    stop("Error: both response and predictors must be data frames")
  }

  pred_names <- names(predictors)
  time_col <- which(names(predictors) == "time")
  combos <- create_df_predictors(names = pred_names[which(pred_names != "time")], n_vars = max_vars)
  if (class(combos) == "character") {
    combos <- data.frame(cov1 = pred_names[which(pred_names != "time")])
  }

  coef_list <- list()
  marginal_pred <- list()
  progress_bar <- txtProgressBar(min = 0, max = nrow(combos), style = 3, char = "=")

  for (i in 1:nrow(combos)) {
    setTxtProgressBar(progress_bar, value = i)
    tmp <- predictors[, c(which(pred_names %in% c(combos[i,], "time")))]
    sub <- dplyr::left_join(as.data.frame(response[, c("time", "dev")]), tmp, by = "time")
    names(sub)[3:ncol(sub)] <- paste0("cov", seq(1, length(3:ncol(sub))))
    covar_names <- names(sub)[3:ncol(sub)]
    name_df <- data.frame(orig_name = names(tmp)[which(names(tmp) != "time")], new_name = covar_names)

    if (model_type == "lm") {
      dummy_f <- as.formula(paste("dev", paste(c("0", covar_names), collapse = " + "), sep = " ~ "))
      f <- dummy_f # this is to make this compatible with mgcv in tidymodels, dumb workaround
      model_spec <- linear_reg() |> set_engine("lm") |> set_mode("regression")
    } else if (model_type == "gam") {
      covar_names_str <- paste0("s(", covar_names, ",k=4,bs='ps')")
      dummy_f <- as.formula(paste("dev", paste(c("0", covar_names), collapse = " + "), sep = " ~ ")) # this is just for tidymodels
      f <- as.formula(paste("dev", paste(c("0", covar_names_str), collapse = " + "), sep = " ~ "))
      model_spec <- gen_additive_mod() |> set_engine("mgcv")|> set_mode("regression")
    } else {
      stop("Unsupported model_type. Use 'lm' or 'gam'")
    }

    # formula needs to get re-added to drop global intercept
    # https://stackoverflow.com/questions/69004818/how-to-fit-a-model-without-an-intercept-using-r-tidymodels-workflow
    recipe <- recipe(dummy_f, data = sub)
    workflow <- workflow() |>
      add_recipe(recipe) |>
      add_model(model_spec, formula = f)

    min_yr <- max(sub$time) - n_forecast + 1
    max_yr <- max(sub$time)

    kfold_splits <- rolling_origin(data = sub, initial = min_yr-1, assess = n_years_ahead)

    # Function to fit a model, predict and calculate RMSE, coefficients, SEs, and marginal effects
    fit_model <- function(split, slice_id, workflow_input) {
      train_data <- training(split)

      if(n_years_ahead == 0) {
        test_data <- train_data
      } else {
        test_data <- testing(split)
      }

      # Fit the workflow to the training data
      fitted_workflow <- workflow_input |> fit(data = train_data)

      # Extract the fitted model from the workflow
      model <- fitted_workflow |> extract_fit_parsnip() |> pluck("fit")

      # Extract model coefficients
      coefficients <- tidy(model)

      # Predictions and SEs on training data
      train_preds <- predict(fitted_workflow, new_data = train_data, type = c("numeric"))#predict(fitted_workflow, new_data = train_data) |>
      ci <- predict(fitted_workflow, new_data = train_data, type = c("conf_int"))
      train_preds <- bind_cols(train_preds, ci)
      train_rmse <- sqrt(mean((train_data$dev - train_preds$.pred)^2))
      train_r2 <- cor(train_data$dev, train_preds$.pred, use="pairwise.complete.obs") ^ 2

      # Predictions and SEs on assessment data
      test_preds <- predict(fitted_workflow, new_data = test_data, type = c("numeric"))
      ci <- predict(fitted_workflow, new_data = test_data, type = c("conf_int"))
      test_preds <- bind_cols(test_preds, ci)
      test_rmse <- sqrt(mean((test_data$dev - test_preds$.pred)^2))
      test_r2 <- cor(test_data$dev, test_preds$.pred, use="pairwise.complete.obs") ^ 2

      # Marginal effects for all covariates
      marg_all <- data.frame()
      for (ii in 1:length(covar_names)) {
        marg <- ggpredict(model, covar_names[ii])
        marg$slice_id <- slice_id
        marg$cov <- name_df$orig_name[ii]
        marg$orig_cov <- name_df$orig_name[ii]
        marg_all <- rbind(marg_all, marg)
      }

      list(
        model = model,
        coefficients = coefficients,
        train_rmse = train_rmse,
        test_rmse = test_rmse,
        train_r2 = train_r2,
        test_r2 = test_r2,
        train_preds = tibble(index = rownames(train_data),
                             pred = train_preds$.pred,
                             pred_lower = train_preds$.pred_lower,
                             pred_upper = train_preds$.pred_upper,
                             #se = train_preds$.pred_lower,
                             slice_id = slice_id,
                             time = train_data$time,
                             type = "train"),
        test_preds = tibble(index = rownames(test_data),
                            pred = test_preds$.pred,
                            pred_lower = test_preds$.pred_lower,
                            pred_upper = test_preds$.pred_upper,
                            slice_id = slice_id,
                            time = test_data$time,
                            type = "test"),
        marginal_effects = marg_all
      )
    }

    # Apply the function to each fold and collect results
    results <- map2(kfold_splits$splits, seq_along(kfold_splits$splits), ~ fit_model(.x, .y, workflow))

    # Extract RMSE and R2 values for train and test data
    train_rmse_values <- map_dbl(results, "train_rmse")
    test_rmse_values <- map_dbl(results, "test_rmse")
    train_r2_values <- map_dbl(results, "train_r2")
    test_r2_values <- map_dbl(results, "test_r2")

    # Combine predictions with original data
    train_preds <- map_df(results, "train_preds")
    test_preds <- map_df(results, "test_preds")

    # Combine coefficients and marginal effects
    coefficients <- map_df(results, "coefficients", .id = "slice_id")
    marginal_effects <- map_df(results, "marginal_effects")

    # Merge the predictions back to the original data frame
    combined_preds <- rbind(train_preds, test_preds)
    dat_with_preds <- sub |> left_join(combined_preds, by = "time")
    indx <- which(names(combos) %in% names(dat_with_preds) == FALSE)
    if(length(indx) > 0) {
      for(ii in 1:length(indx)) dat_with_preds[[names(combos)[indx[ii]]]] <- NA
    }
    # print(train_rmse_values) # Print the RMSE values
    # print(test_rmse_values) # Print the RMSE values
    # print(dat_with_preds)# Print the data with predictions
    # print(coefficients)# Print coefficients
    # print(marginal_effects)# Print marginal effects

    # Extract RMSE values for train and test data
    summary_df <- data.frame(i = i, train_rmse_values = train_rmse_values,
                             test_rmse_values = test_rmse_values,
                             train_r2_values = train_r2_values,
                             test_r2_values = test_r2_values)
    dat_with_preds$i <- i
    coefficients$i <- i
    marginal_effects$i <- i

    if(i == 1) {
      all_preds <- dat_with_preds
      summary <- summary_df
      all_coefficients <- coefficients
      all_marginal_effects <- marginal_effects
    } else {
      all_preds <- rbind(all_preds, dat_with_preds)
      summary <- rbind(summary, summary_df)
      all_coefficients <- rbind(all_coefficients, coefficients)
      all_marginal_effects <- rbind(all_marginal_effects, marginal_effects)
    }

  }

  return(list(pred = all_preds, summary = summary, coefs = all_coefficients, marginals = all_marginal_effects))
}
