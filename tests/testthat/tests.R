test_that("1D linear model works", {
  set.seed(123)

  response <- data.frame(time = 1:40, dev = rnorm(40))

  predictors <- matrix(rnorm(400), ncol = 10)
  #colnames(predictors) = paste0("X",1:ncol(predictors))
  predictors <- as.data.frame(predictors)
  predictors$time <- 1:40
  lm_example <- univariate_forecast(response,
       predictors,
       model_type = "lm",
       n_forecast = 0,
       n_years_ahead = 0,
       max_vars = 3)
  expect_equal("pred", names(lm_example)[1])
  expect_equal("summary", names(lm_example)[2])
  expect_equal("coefs", names(lm_example)[3])

  sub <- predictors[1:40,]
  sub$dev <- response$dev[1:40]
  f <- lm(dev[1:40] ~ -1 + V1 + V2 + V3, data=sub)
  pred <- predict(f, newdata=predictors[40,1:3])

  predicted <- lm_example$pred |>
    dplyr::filter(i==1, type=="test", time==40)

  expect_equal(predicted$pred, as.numeric(pred), tolerance = 0.001)

  # Evaluate 1-step ahead comparisons
  lm_example <- univariate_forecast(response,
                                    predictors,
                                    model_type = "lm",
                                    n_forecast = 1,
                                    n_years_ahead = 1,
                                    max_vars = 3)
  sub <- predictors[1:39,]
  sub$dev <- response$dev[1:39]
  f <- lm(dev[1:39] ~ -1 + V1 + V2 + V3, data=sub)
  pred <- predict(f, newdata=predictors[40,1:3])
  expect_equal(lm_example$pred$pred[40], as.numeric(pred), tolerance = 0.001)

})

test_that("1D gam model works", {
  set.seed(123)

  response <- data.frame(time = 1:40, dev = rnorm(40))

  predictors <- matrix(rnorm(400), ncol = 10)
  #colnames(predictors) = paste0("X",1:ncol(predictors))
  predictors <- as.data.frame(predictors)
  predictors$time <- 1:40
  gam_example <- univariate_forecast(response,
                                    predictors,
                                    model_type = "gam",
                                    n_forecast = 0,
                                    n_years_ahead = 0,
                                    max_vars = 3)

  predicted <- gam_example$pred |>
    dplyr::filter(i==1, type=="test", time==40)

  expect_equal(predicted$pred, -0.08842347, tolerance = 0.001)

  set.seed(123)
  gam_example <- univariate_forecast(response,
                                     predictors,
                                     model_type = "gam",
                                     n_forecast = 4,
                                     n_years_ahead = 0,
                                     max_vars = 3)

  predicted <- gam_example$pred |>
    dplyr::filter(i==1, type=="test", time > 36, time == 37, slice_id==2)
  expect_equal(predicted$pred, -0.0003308037, tolerance = 0.001)

  predicted <- gam_example$pred |>
    dplyr::filter(i==1, type=="test", time > 36, time == 38, slice_id==3)
  expect_equal(predicted$pred, 0.2696979393, tolerance = 0.001)

  predicted <- gam_example$pred |>
    dplyr::filter(i==1, type=="test", time > 36, time == 39, slice_id==4)
  expect_equal(predicted$pred, -0.0094620398, tolerance = 0.001)

  predicted <- gam_example$pred |>
    dplyr::filter(i==1, type=="test", time > 36, time == 40, slice_id==5)
  expect_equal(predicted$pred, -0.0884234683, tolerance = 0.001)

})

#
# test_that("1D lasso model works", {
#   set.seed(123)
#
#   response <- data.frame(time = 1:40, dev = rnorm(40))
#
#   predictors <- matrix(rnorm(400), ncol = 10)
#   #colnames(predictors) = paste0("X",1:ncol(predictors))
#   predictors <- as.data.frame(predictors)
#   predictors$time <- 1:40
#   lm_example <- univariate_forecast_ml(response,
#        predictors,
#        model_type = "glmnet",
#        n_forecast = 1,
#        n_years_ahead = 1)
#   expect_equal(lm_example$pred$est[40], -0.3595902, tolerance = 0.001)
#
#   set.seed(123)
#   lm_example <- univariate_forecast_ml(response,
#                                        predictors,
#                                        model_type = "glmnet",
#                                        n_forecast = 4,
#                                        n_years_ahead = 1)
#   expect_equal(lm_example$pred$est[37], -0.2298422, tolerance = 0.001)
#   expect_equal(lm_example$pred$est[38], -0.3201412, tolerance = 0.001)
#   expect_equal(lm_example$pred$est[39], 1.3507670, tolerance = 0.001)
#   expect_equal(lm_example$pred$est[40], -0.3595902, tolerance = 0.001)
# })
#
#
# test_that("1D random forest model works", {
#   set.seed(123)
#
#   response <- data.frame(time = 1:40, dev = rnorm(40))
#
#   predictors <- matrix(rnorm(400), ncol = 10)
#   #colnames(predictors) = paste0("X",1:ncol(predictors))
#   predictors <- as.data.frame(predictors)
#   predictors$time <- 1:40
#   lm_example <- univariate_forecast_ml(response,
#                                        predictors,
#                                        model_type = "randomForest",
#                                        n_forecast = 1,
#                                        n_years_ahead = 1)
#   expect_equal(lm_example$pred$est[40], -0.02658667, tolerance = 0.001)
#
#   set.seed(123)
#   lm_example <- univariate_forecast_ml(response,
#                                        predictors,
#                                        model_type = "randomForest",
#                                        n_forecast = 4,
#                                        n_years_ahead = 1)
#   expect_equal(lm_example$pred$est[37], -0.07783784, tolerance = 0.001)
#   expect_equal(lm_example$pred$est[38], 0.1210044, tolerance = 0.001)
#   expect_equal(lm_example$pred$est[39], 0.4661449, tolerance = 0.001)
#   expect_equal(lm_example$pred$est[40], 0.1016875, tolerance = 0.001)
# })
