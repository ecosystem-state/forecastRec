#' Forecast recruitment deviations as a function of covariates using GLMs and GAMs
#'
#' \code{multivariate_forecast} Takes data frame or matrix of responses and dataframe of predictors and automates predictions
#'
#' @param response A matrix or data frame of responses for the modeling (values to be forecast) containing
#' a "time" column, "species" column, and "dev" column
#' @param predictors A data frame of predictors used for forecasting recruitment
#' @param model_type The type of model used to link predictors to forecasted recruitment. Can
#' be "lm", "gam",
#' @param n_forecast How many years to use as a holdout / test set
#' @param n_years_ahead How many years ahead to forecast (defaults to 1)
#' 1:n_vars variables, and then results are combined and sorted to remove duplicates
#' @param max_vars The maximum number of variables to include as predictors; defaults to 3
#' @param formula Optional formula for passing to gam(), lmer(), randomForest(), etc.
#' @import tidymodels
#' @import multilevelmod
#' @importFrom parsnip set_mode set_engine linear_reg gen_additive_mod fit
#' @importFrom rsample rolling_origin
#' @importFrom recipes recipe
#' @importFrom tune control_resamples fit_resamples
#' @importFrom purrr map_dfr map_dbl
#' @importFrom workflows add_model add_recipe workflow extract_fit_parsnip
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join arrange mutate row_number bind_rows
#' @importFrom broom tidy
#' @importFrom stats lm as.formula predict cor
#' @importFrom mgcv gam
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom tibble tibble
#' @importFrom yardstick metric_set rmse rsq
#' @importFrom ggeffects ggpredict
#' @importFrom stats terms
#'
#' @export
#'
#' @return a list containing predictions and the variables used in making predictions
#'
multivariate_forecast = function(response,
                               predictors,
                               model_type,
                               n_forecast = 10,
                               n_years_ahead = 1,
                               max_vars = 3,
                               formula = NULL) {

  # create a dataframe of predictors
  pred_names = names(predictors)
  time_col = which(names(predictors)=="time")

  combos = create_df_predictors(names = pred_names[which(pred_names!="time")],
                                n_vars = max_vars)
  if(class(combos)=="character") {
    # catch the case where a single character is returned
    combos = data.frame(cov1 = pred_names[which(pred_names!="time")])
  }

  coef_list <- list() # empty list for storing coefficients
  marginal_pred <- list() # empty list for storing marginal predictions

  # add progress bar
  progress_bar <- txtProgressBar(min = 0, max = nrow(combos), style = 3, char = "=")

  for(i in 1:nrow(combos)) {
    setTxtProgressBar(progress_bar, value = i)
    # keep time and
    tmp <- predictors[, c(which(pred_names %in% c(combos[i,], "time")))]
    sub = dplyr::left_join(as.data.frame(response[,c("time","dev","species")]), tmp, by="time")
    # remove spaces if they exist to help with formula parsing
    names(sub)[4:ncol(sub)] = paste0("cov", seq(1,length(4:ncol(sub))))
    covar_names = names(sub)[4:ncol(sub)]
    name_df = data.frame(orig_name = names(tmp)[which(names(tmp) != "time")], new_name = covar_names)

    # Add formulas -- no intercept because rec devs are already standardized, as are predictors
    if(model_type=="lm") {
      dummy_f <- as.formula(paste("dev",
                            paste(c("0",paste0("species:",covar_names)), collapse = " + "),
                            sep = " ~ "))
      # combined species:covariate because ":" isn't allowed in tidymodels recipe
      model_mat <- model.matrix(dummy_f, sub)
      colnames(model_mat) <- gsub(":", "_", colnames(model_mat))
      sub <- cbind(sub, model_mat)

      dummy_f <- as.formula(paste("dev",
                                  paste(colnames(model_mat), collapse = " + "),
                                  sep = " ~ "))
      mod_form <- dummy_f # this is to make this compatible with mgcv in tidymodels, dumb workaround
      model_spec <- linear_reg() |> set_engine("lm") |> set_mode("regression")
    }
    if(model_type=="gam") {
      sub$species = as.factor(sub$species)
      covar_names_str <- paste0("s(",covar_names,", species,k=4,bs='fs',m=2)")
      dummy_f <- as.formula(paste("dev",
                            paste(c("0",covar_names), collapse = " + "),
                            sep = " ~ "))
      mod_form <- as.formula(paste("dev", paste(c("0", covar_names_str), collapse = " + "), sep = " ~ "))
      model_spec <- gen_additive_mod() |> set_engine("mgcv")|> set_mode("regression")
    }
    if(model_type=="glmm") {
      sub$species = as.factor(sub$species)
      covar_names_str = paste0("(0+",covar_names," | species)")
      mod_form <- as.formula(paste("dev",
                            paste(c("0",covar_names_str), collapse = " + "),
                            sep = " ~ "))
      model_spec <- model_spec <- linear_reg() |> set_engine("lmer") |> set_mode("regression")
    }
    # annoying -- but formula needs to get re-added to drop global intercept
    # https://stackoverflow.com/questions/69004818/how-to-fit-a-model-without-an-intercept-using-r-tidymodels-workflow
    recipe <- recipe(dev ~ ., data = sub) # empty recipe

    workflow <- workflow() |>
      add_recipe(recipe) |>
      add_model(model_spec, formula = mod_form)

    if(n_forecast == 0 && n_years_ahead == 0) {
      fitted_workflow <- suppressMessages(suppressWarnings(fit(workflow, data = sub)))

      # Extract the fitted model
      fitted_model <- extract_fit_parsnip(fitted_workflow)

      # Predictions and RMSE/R2 calculations on the full dataset
      pred <- predict(fitted_model, new_data = sub)$.pred
      if(model_type == "glmm") pred <- as.numeric(predict(fitted_model$fit, newdata=sub))
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
      covariate_names <- covariate_names[which(covariate_names != "species")] # drop factors for GAMs

      marg_all <- data.frame()
      if(model_type != "glmm") {
        for (ii in 1:length(covariate_names)) {
          marg <- ggpredict(fitted_model$fit, covariate_names[ii])
          marg$cov_name <- covariate_names[ii] # add covariate name
          marg_all <- rbind(marg_all, marg)
        }
        marg_all$i <- i
      }

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

        resid <- as.numeric(residuals(fitted$fit))#fitted$fit$residuals
        pred <- as.numeric(fitted(fitted$fit))#fitted$fit$fitted.values
        y <- pred + resid
        rmse_train <- sqrt(mean(resid^2))
        rsq_train <- cor(y,pred)^2

        marg_all <- data.frame()
        if(model_type != "glmm") {
          covariate_names <- attr(terms(fitted$fit), "term.labels")
          covariate_names <- covariate_names[which(covariate_names != "species")]
          for (ii in 1:length(covariate_names)) {
            marg <- ggpredict(fitted$fit, covariate_names[ii])
            marg$cov_name <- covariate_names[ii] # add covariate name
            marg_all <- bind_rows(marg_all, marg)
          }
        }

        #mod <- NULL
        #coef <- NULL
        #if(class(fitted$fit) !="lmer") mod = fitted$fit$model

        return(list(#model = mod,
                    #coefficients <- coef,
                    #coefficients = fitted$fit$coefficients,
                    marginals = marg_all,
                    pred = pred,
                    rmse_train = rmse_train,
                    rsq_train = rsq_train))
      }

      control <- control_resamples(save_pred = TRUE, extract = fit_model)

      results <- suppressWarnings(suppressMessages(fit_resamples(workflow,
                                                resamples = kfold_splits,
                                                metrics = metric_set(rmse, rsq),
                                                control = control)))

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
                                 .config = test_values$.config[1]) |>
        arrange(id, .metric)
      train_values$type <- "train"

      # marginal effects from training models
      marg_all <- map_dfr(results$.extracts, function(df) {
        # Extract the coefficients directly from the nested structure
        df$.extracts[[1]]$marginals
      })
      if(nrow(marg_all) > 0) marg_all$i <- i # null for lmer objects

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
            pred_obj <- pred_obj |>
              mutate(row = row_number()) |>
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

      # Merge the predictions back to the original data frame
      combined_preds <- bind_rows(train_pred, test_pred)
      dat_with_preds <- sub |> left_join(combined_preds, by = "time")

      dat_with_preds$i <- i

      # Extract RMSE values for train and test data
      summary_df <- bind_rows(test_values, train_values)
      summary_df$i <- i
    }
    #marginal_effects$i <- i

    if(i == 1) {
      all_preds <- dat_with_preds
      summary <- summary_df
      #all_coefficients <- coefficients
      all_marginal_effects <- marg_all
    } else {
      all_preds <- bind_rows(all_preds, dat_with_preds)
      summary <- bind_rows(summary, summary_df)
      #all_coefficients <- rbind(all_coefficients, coefficients)
      all_marginal_effects <- bind_rows(all_marginal_effects, marg_all)
    }

  }

  # add a unique idenfier for the marginal effects across variables and folds
  if(model_type != "glmm") {
    breaks <- c(1, ifelse(diff(all_marginal_effects$x) < 0, 1, 0))
    all_marginal_effects$var_fold <- cumsum(breaks)
  }
  return(list(pred = all_preds,
              summary = summary,
              #coefs = all_coefficients,
              covariate_combos = combos,
              marginal_effects = all_marginal_effects))
}
