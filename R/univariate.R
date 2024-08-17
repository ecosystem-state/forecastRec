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
#' @importFrom rsample rolling_origin
#' @importFrom recipes recipe
#' @importFrom tune control_resamples fit_resamples
#' @importFrom purrr map_dfr map_dbl
#' @importFrom workflows add_model add_recipe workflow extract_fit_parsnip
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join
#' @importFrom broom tidy
#' @importFrom stats lm as.formula predict cor
#' @importFrom mgcv gam
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom tibble tibble
#' @importFrom yardstick metric_set rmse rsq
#' @importFrom ggeffects ggpredict
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

    # annoying -- but formula needs to get re-added to drop global intercept
    # https://stackoverflow.com/questions/69004818/how-to-fit-a-model-without-an-intercept-using-r-tidymodels-workflow
    recipe <- recipe(dummy_f, data = sub)
    workflow <- workflow() |>
      add_recipe(recipe) |>
      add_model(model_spec, formula = f)

    if(n_forecast == 0 && n_years_ahead == 0) {
      fitted_workflow <- fit(workflow, data = sub)

      # Extract the fitted model
      fitted_model <- extract_fit_parsnip(fitted_workflow)

      # Predictions and RMSE/R2 calculations on the full dataset
      pred <- predict(fitted_model, new_data = sub)$.pred
      resid <- sub$dev - pred
      rmse_train <- sqrt(mean(resid^2))
      rsq_train <- cor(sub$dev, pred)^2

      # Collect coefficients and predictions
      # this gets messy with everything that isn't lms
      # coefficients <- tidy(fitted_model$fit) |>
      #   select(term, estimate) |>
      #   tidyr::pivot_wider(names_from = term, values_from = estimate)
      # coefficients$i <- i

      dat_with_preds <- sub
      dat_with_preds$.pred <- pred
      indx <- which(names(combos) %in% names(dat_with_preds) == FALSE)
      if(length(indx) > 0) {
        for(ii in 1:length(indx)) {
          #coefficients[[names(combos)[indx[ii]]]] <- NA
          dat_with_preds[[names(combos)[indx[ii]]]] <- NA
        }
      }
      dat_with_preds$i <- i

      # Extract marginals
      covariate_names <- attr(terms(fitted_model$fit), "term.labels")
      marg_all <- data.frame()
      for (ii in 1:length(covariate_names)) {
        marg <- ggpredict(fitted_model$fit, covariate_names[ii])
        marg_all <- rbind(marg_all, marg)
      }
      marg_all$i <- i

      summary_df <- data.frame(id = "full_dataset",
                               .metric = c("rmse", "rsq"),
                               .estimator = "standard",
                               .estimate = c(rmse_train, rsq_train),
                               type = "train")
      summary_df$i <- i
    } else {
      min_yr <- max(sub$time) - n_forecast + 1
      max_yr <- max(sub$time)

      kfold_splits <- rolling_origin(data = sub, initial = min_yr-1, assess = n_years_ahead)
      # Function to fit a model, predict and calculate RMSE, coefficients, SEs, and marginal effects
      fit_model <- function (x) {
        fitted <- extract_fit_parsnip(x)
        resid <- fitted$fit$residuals
        pred <- fitted$fit$fitted.values
        y <- pred + resid
        rmse_train <- sqrt(mean(resid^2))
        rsq_train <- cor(y,pred)^2

        covariate_names <- attr(terms(fitted$fit), "term.labels")
        marg_all <- data.frame()
        for (ii in 1:length(covariate_names)) {
          marg <- ggpredict(fitted$fit, covariate_names[ii])
          marg_all <- rbind(marg_all, marg)
        }

        return(list(model = fitted$fit$model,
                    #coefficients = fitted$fit$coefficients,
                    marginals = marg_all,
                    pred = pred,
                    rmse_train = rmse_train,
                    rsq_train = rsq_train))
      }

      control <- control_resamples(save_pred = TRUE, extract = fit_model)

      results <- suppressMessages(fit_resamples(workflow,
                                                resamples = kfold_splits,
                                                metrics = metric_set(rmse, rsq),
                                                control = control))

      # Extract RMSE and R2 values for test data
      test_values <- tune::collect_metrics(results, summarize = FALSE)
      test_values$type <- "test"
      test_rmse_values <- dplyr::filter(test_values,
                                        .metric=="rmse")
      test_r2_values <- dplyr::filter(test_values,
                                      .metric=="rsq")
      # Extract RMSE and R2 values for training data
      rmse_train_values <- map_dbl(results$.extracts, function(df) {
        df$.extracts[[1]]$rmse_train
      })
      rsq_train_values <- map_dbl(results$.extracts, function(df) {
        df$.extracts[[1]]$rsq_train
      })
      train_values <- data.frame(id = rep(results$id,2),
                                 .metric = sort(rep(c("rmse","rsq"), length(results$id))),
                                 .estimator = "standard",
                                 .estimate = c(rmse_train_values, rsq_train_values),
                                 .config = test_values$.config[1]) %>%
        arrange(id, .metric)
      train_values$type <- "train"

      # marginal effects from training models
      marg_all <- map_dfr(results$.extracts, function(df) {
        # Extract the coefficients directly from the nested structure
        df$.extracts[[1]]$marginals
      })
      marg_all$id <- results$id
      marg_all$i <- i

      # coefficients from training models
      # coefficients <- map_dfr(results$.extracts, function(df) {
      #   # Extract the coefficients directly from the nested structure
      #   df$.extracts[[1]]$coefficients
      # })
      # coefficients$id <- results$id
      # coefficients$i <- i

      # Collect predictions
      test_pred <- tune::collect_predictions(results, summarize = FALSE) |>
        dplyr::select(-dev,-.config) |>
        dplyr::rename(time = .row)
      test_pred$type <- "test"

      # train_pred_df <- map_dfr(results$.extracts, function(df) {
      #   # Extract the predictions directly from the nested structure
      #   df$.extracts[[1]]$pred
      # })
      train_pred_df <- map_dfr(results$.extracts, function(df) {
        pred_obj <- df$.extracts[[1]]$pred

        # If pred_obj is a vector (like from mgcv), we need to reshape it
        if (is.vector(pred_obj)) {
          # Convert to a tibble and add row identifiers if needed
          pred_obj <- tibble(.pred = pred_obj)
          if (ncol(pred_obj) == 1) {
            pred_obj <- pred_obj %>%
              mutate(row = row_number()) %>%
              pivot_wider(names_from = row, values_from = .pred)
          }
        }
        # If pred_obj is already a data frame (like from lm), just return it
        if (is.data.frame(pred_obj)) {
          return(pred_obj)
        }
        # If structure is still not correct, return an empty tibble
        return(tibble())
      })
      train_pred <- train_pred_df |>
        mutate(id = paste0("Slice",row_number())) |>  # Create a row ID to identify each row
        pivot_longer(
          cols = -id,  # All columns except row_id pivoted to long format
          names_to = "time",  # Name of the new column
          values_to = ".pred"  # Name of the new column
        )
      train_pred$time <- as.numeric(train_pred$time)
      train_pred$type <- "train"

      # Combine coefficients and marginal effects
      #marginal_effects <- map_df(results, "marginal_effects")

      # Merge the predictions back to the original data frame
      combined_preds <- rbind(train_pred, test_pred)
      dat_with_preds <- sub |> left_join(combined_preds, by = "time")
      indx <- which(names(combos) %in% names(dat_with_preds) == FALSE)
      if(length(indx) > 0) {
        for(ii in 1:length(indx)) {
          #coefficients[[names(combos)[indx[ii]]]] <- NA
          dat_with_preds[[names(combos)[indx[ii]]]] <- NA
        }
      }
      dat_with_preds$i <- i

      # Extract RMSE values for train and test data
      summary_df <- rbind(test_values, train_values)
      summary_df$i <- i
    }
    #marginal_effects$i <- i

    if(i == 1) {
      all_preds <- dat_with_preds
      summary <- summary_df
      #all_coefficients <- coefficients
      all_marginal_effects <- marg_all
    } else {
      all_preds <- rbind(all_preds, dat_with_preds)
      summary <- rbind(summary, summary_df)
      #all_coefficients <- rbind(all_coefficients, coefficients)
      all_marginal_effects <- rbind(all_marginal_effects, marg_all)
    }

  }

  return(list(pred = all_preds,
              summary = summary,
              #coefs = all_coefficients,
              covariate_combos = combos,
              marginal_effects = all_marginal_effects))
}
