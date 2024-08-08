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
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr left_join
#' @importFrom broom tidy
#' @importFrom ggeffects ggpredict
#' @importFrom stats lm as.formula predict cor model.matrix na.pass
#' @importFrom mgcv gam
#' @importFrom utils txtProgressBar setTxtProgressBar
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

    sub$est <- NA
    sub$se <- NA
    sub$train_r2 <- NA
    sub$train_rmse <- NA
    sub$id <- i
    min_yr <- max(sub$time) - n_forecast + 1
    max_yr <- max(sub$time)

    # This is the the leave future out cross validation, which could replace the loop below
    # The major challenges with implementing this are
    # 1. tidymodels won't let us generate prediction metrics for the training dataset
    # 2. tidymodels won't return things like standard errors on predictions (for training or test data)
    # kfold_splits <- rsample::rolling_origin(data = sub,
    #                                   initial = min_yr - 1, #
    #                                   assess = 1)
    # # Fit the workflow to each resample
    # resample_results <- fit_resamples(
    #   workflow,
    #   resamples = kfold_splits,
    #   metrics = metric_set(rmse, rsq),
    #   control = control_resamples(save_pred = TRUE)
    # )
    # predictions
    # predictions <- resample_results %>%
    #   collect_predictions()
    # Collect the metrics
    # metrics <- resample_results %>%
    #   collect_metrics(summarize=FALSE) # if TRUE, averaged across slices

    for (yr in min_yr:max_yr) {
      sub_dat <- sub[which(sub$time < yr - n_years_ahead + 1),]
      fit <- try(fit(workflow, data = sub_dat), silent = TRUE)

      if (class(fit)[1] != "try-error") {
        # predict the mean
        pred <- try(predict(fit, new_data = sub[which(sub$time == yr),]), silent = TRUE)
        sub$est[which(sub$time == yr)] <- pred$.pred
        # predict SEs
        pred <- try(predict(fit, new_data = sub[which(sub$time == yr),], type = "conf_int"), silent = TRUE)
        sub$se[which(sub$time == yr)] <- (pred$.pred_upper - pred$.pred_lower) / 2

        # predict the mean of the training data
        pred_train <- try(predict(fit, sub_dat), silent = TRUE)
        # calculate the rmse and r2 for the training data
        if (class(pred_train)[1] != "try-error") {
          sub$train_r2[which(sub$time == yr)] <- cor(c(sub_dat$dev), pred_train$.pred, use = "pairwise.complete.obs") ^ 2
          sub$train_rmse[which(sub$time == yr)] <- sqrt(mean((c(sub_dat$dev) - pred_train$.pred)^2, na.rm = TRUE))
        }

        coefs <- broom::tidy(fit$fit$fit)
        coefs$yr <- yr
        coefs$orig_cov <- name_df$orig_name

        if (yr == min_yr) {
          all_coefs <- coefs
        } else {
          all_coefs <- rbind(all_coefs, coefs)
        }

        marg_all <- data.frame()
        for (ii in 1:length(covar_names)) {
          marg <- ggpredict(fit$fit$fit, covar_names[ii])
          marg$year <- yr
          marg$cov <- name_df$orig_name[ii]
          marg$orig_cov <- name_df$orig_name[ii]
          marg_all <- rbind(marg_all, marg)
        }
      }
    }

    marginal_pred[[i]] <- marg_all
    coef_list[[i]] <- all_coefs

    if (i == 1) {
      out_df <- sub
    } else {
      missing_covars <- names(out_df)[which(names(out_df) %in% names(sub) == FALSE)]
      if (length(missing_covars) > 0) {
        for (mm in 1:length(missing_covars)) {
          sub[[missing_covars[mm]]] <- NA
        }
      }
      out_df <- rbind(out_df, sub)
    }
  }

  combos$id <- seq(1, nrow(combos))
  out <- list("pred" = out_df,
              "vars" = combos,
              "coefs" = coef_list,
              "marginal" = marginal_pred)
  return(out)
}
