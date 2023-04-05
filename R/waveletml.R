#' @title Wavelet Decomposition-Based ARIMA-GARCH-ANN Hybrid Modeling
#'
#' @param Y Univariate time series
#' @param ratio Ratio of number of observations in training and testing sets
#' @param n_lag  Lag of the provided time series data
#' @param l  Level of decomposition
#' @param f  Filter of decomposition
#' @import utils, stats, wavelets, tseries, forecast, fGarch, aTSA, FinTS, LSTS, earth, caret, neuralnet, pso
#' @return
#' \itemize{
#'   \item Train_fitted: Train fitted result
#'   \item Test_predicted: Test predicted result
#'   \item Accuracy: Accuracy
#'   }
#' @export
#'
#' @examples
#' Y <- rnorm(100, 100, 10)
#' result <- warigaan(Y, ratio = 0.8, n_lag = 4)
#' @references
#' \itemize{
#'   \item Paul, R. K., & Garai, S. (2021). Performance comparison of wavelets-based machine learning technique for forecasting agricultural commodity prices. Soft Computing, 25(20), 12857-12873.
#'   \item Paul, R. K., & Garai, S. (2022). Wavelets based artificial neural network technique for forecasting agricultural prices. Journal of the Indian Society for Probability and Statistics, 23(1), 47-61.
#'   \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
#'   }

warigaan <- function(Y, ratio = 0.9, n_lag = 4, l = 6, f = 'haar'){
  optimize_weights<-NULL
  objective<-NULL
  all_metrics<-NULL
  # some default values have been set
  # Y is a vector
  # ratio is train:test
  # n_lag is number of lags in the data
  # Wavelet parameters #### (Level = l, Filter = f)
  # Level <- c(2:floor(log2(length(Y))))
  # Filter<-c('haar', 'bl14','la8', 'c6')
  ######################################
  # Data preparation ####

  # embedding for finding log return
  diff.1 <- embed(Y, 2)
  Yt <- diff.1[,1]
  Yt_1 <- diff.1[,2]
  y <- log(Yt/Yt_1)
  # Compute the average of non-zero contents in the data
  nonzero_data <- y[y != 0 & !is.na(y)]
  average_nonzero <- mean(nonzero_data)

  # Replace NaNs with the average of non-zero contents in the data
  y[is.nan(y)] <- average_nonzero

  # Check the result
  y
  # embedding for finding lag series of actual series
  #n_lag <- 4
  embed_size_y <- n_lag+1 # same lag (n_lag_y-1) for every sub series
  diff.2 <- embed(Yt,embed_size_y)
  dim(diff.2)
  Y_actual <- diff.2[,1]
  Y_actual_1 <- diff.2[,2]
  # train-test split
  #ratio <- 0.8
  n <- length(Y_actual) # this is already (original length-embed_y)
  Y_train <- Y_actual[1:(n*ratio)]
  Y_test <- Y_actual[(n*ratio+1):n]
  Y_train_1 <- Y_actual_1[1:(n*ratio)]
  Y_test_1 <- Y_actual_1[(n*ratio+1):n]
  # embedding for finding lag series of log return series
  diff.3 <- embed(y, embed_size_y)
  y_actual <- diff.3[,1]
  y_train <- y_actual[1:(n*ratio)]
  y_test <- y_actual[(n*ratio+1):n]

  # Data pre-processing ####
  #f=Filter[1]
  #l=Level[1]
  # wavelet
  wavelet_y <- modwt(y, filter=f,n.levels=l)
  wavelet_df <- w<-cbind(as.data.frame(wavelet_y@W),V=as.data.frame(wavelet_y@V)[,l])
  # lag matrices for wavelet decomposed series
  wavelet_lag_matrices <- list()
  for (col_name in colnames(wavelet_df)) {
    embedded <- embed(wavelet_df[[col_name]], embed_size_y)
    wavelet_lag_matrices_name <- paste0("embed_", col_name)
    wavelet_lag_matrices[[wavelet_lag_matrices_name]] <- embedded
  }
  wavelet_lag_train <- list()
  for (names in names(wavelet_lag_matrices)) {
    wavelet_train <- wavelet_lag_matrices[[names]][1:(n*ratio),]
    wavelet_lag_train_name <- names
    wavelet_lag_train[[wavelet_lag_train_name]] <- wavelet_train
  }

  wavelet_lag_test <- list()
  for (names in names(wavelet_lag_matrices)) {
    wavelet_test <- wavelet_lag_matrices[[names]][(n*ratio+1):n,]
    wavelet_lag_test_name <- names
    wavelet_lag_test[[wavelet_lag_test_name]] <- wavelet_test
  }

  # model fitting ####

  # ARIMA
  arima <- auto.arima(wavelet_lag_train$embed_V[,1])
  order <- arimaorder(arima)
  model_arima<-arima(wavelet_train[,ncol(wavelet_train)],order=c(order[1], order[2], order[3]))
  pred_arima <-arima$fitted
  forecast_arima <- data.frame(predict(arima,n.ahead=((n-(n*ratio)))))
  forecast_arima <-forecast_arima$pred

  #GARCH Model ####
  ARCH_pvalue <- as.numeric(FinTS::ArchTest(model_arima$residuals)$p.value)
  #ARCH_pvalue<-1
  if(ARCH_pvalue<=0.05){
    garch.fit <- fGarch::garchFit(~garch(1, 1), data = y_train, trace = FALSE)
    pred_V <- garch.fit@fitted
    forecast_V <- predict(garch.fit, n.ahead=(n-(n*ratio)))
    forecast_V <- forecast_V$meanForecast

    Resid_V <- garch.fit@residuals
    for_resid<-as.ts(y_test-forecast_V)
  }else {
    pred_V <- pred_arima
    forecast_V <- forecast_arima

    Resid_V <- as.ts(model_arima$residuals)
    for_resid<-as.vector(y_test-as.vector(forecast_arima))
  }

  names(for_resid)<-"Resid"
  names(Resid_V)<-"Resid"
  bl.test <- Box.Ljung.Test(Resid_V) #not white noise
  if(max(bl.test$data$y)<=0.05){
    Resid_ind<-1
  } else {
    Resid_ind<-2
  }

  # mars_data preparation
  mars_train <- list()
  for (names in names(wavelet_lag_train)) {
    wavelet_train_data <- cbind(wavelet_lag_train[[names]][,1])
    wavelet_train_data_name <- names
    mars_train[[wavelet_train_data_name]] <- wavelet_train_data
  }
  mars_test <- list()
  for (names in names(wavelet_lag_test)) {
    wavelet_test_data <- cbind(wavelet_lag_test[[names]][,1])
    wavelet_test_data_name <- names
    mars_test[[wavelet_test_data_name]] <- wavelet_test_data
  }

  # if arima/garch residual is not white noise it will be used as mars predictor
  # if(Resid_ind==2){
  #   predictors_train <- cbind(do.call(cbind, head(mars_train, -1)), as.vector(Resid_V))
  #colnames(predictors_train)<-c(names(mars_train[-length(mars_train)]), "Resid")
  #   predictors_test <- cbind(do.call(cbind, head(mars_test, -1)), as.vector(for_resid))
  #colnames(predictors_test)<-c(names(mars_train[-length(mars_train)]), "Resid")
  # } else {
  predictors_train<-data.frame(do.call(cbind, head(mars_train, -1)))
  colnames(predictors_train)<-c(names(mars_train[-length(mars_train)]))
  predictors_test<-data.frame(do.call(cbind, head(mars_test, -1)))
  colnames(predictors_test)<-c(names(mars_train[-length(mars_train)]))
  # }
  # mars
  data_MARS<-data.frame(cbind(y=y_train,predictors_train))

  MARS <- earth(y~., data = data_MARS)

  selcted_var<-MARS$namesx

  # train and test data for ann and svr
  data_train <- wavelet_lag_train[selcted_var]
  data_test <- wavelet_lag_test[selcted_var]
  for(i in seq_along(data_train)) {
    # determine the name of the current element
    list_name <- names(data_train)[i]

    # determine the number of columns for the current matrix
    num_cols <- ncol(data_train[[i]])

    # create a vector of column names for the current matrix
    col_names <- paste(list_name, "col", 1:num_cols, sep = "_")

    # assign the column names to the current matrix
    colnames(data_train[[i]]) <- col_names
  }
  for(i in seq_along(data_test)) {
    # determine the name of the current element
    list_name <- names(data_test)[i]

    # determine the number of columns for the current matrix
    num_cols <- ncol(data_test[[i]])

    # create a vector of column names for the current matrix
    col_names <- paste(list_name, "col", 1:num_cols, sep = "_")

    # assign the column names to the current matrix
    colnames(data_test[[i]]) <- col_names
  }

  # ann
  # create empty data frames to store predicted values for each element's response
  train_predicted_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(train_predicted_df) <- c("element", "train_predicted")
  test_predicted_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(test_predicted_df) <- c("element", "test_predicted")
  for (i in seq_along(data_train)) {
    # extract the predictors and response from the current matrix
    predictors <- data_train[[i]][, -1]
    response <- data_train[[i]][, 1]
    allvars <- colnames(predictors)
    predictorvars <- paste(allvars, collapse="+")
    form <- as.formula(paste0(names(data_train)[i], "_col_1~", predictorvars))

    # train a neural network with one hidden layer
    nn <- neuralnet(form, data_train[[i]], hidden = c(4))

    # predict the response for the current matrix in the training set
    train_predicted <- predict(nn, predictors)

    # create a data frame with the element name and predicted values for the training set
    train_predicted_element_df <- data.frame(element = rep(names(data_train[i]), nrow(train_predicted)),
                                             train_predicted = train_predicted)

    # add the predicted values for the current element in the training set to the overall data frame
    train_predicted_df <- rbind(train_predicted_df, train_predicted_element_df)
    # save the model for this element
    assign(paste0("nn_", names(data_train[i])), nn)
  }
  # create empty data frames to store predicted values for each element's response
  test_predicted_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(test_predicted_df) <- c("element", "test_predicted")
  test_predicted_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(test_predicted_df) <- c("element", "test_predicted")
  for (i in seq_along(data_test)) {
    # extract the predictors and response from the current matrix
    predictors <- data_test[[i]][, -1]
    response <- data_test[[i]][, 1]

    # predict the response for the current matrix in the testing set
    nn <- get(paste0("nn_", names(data_test[i])))
    test_predicted <- predict(nn, predictors)

    # create a data frame with the element name and predicted values for the testing set
    test_predicted_element_df <- data.frame(element = rep(names(data_test[i]), nrow(test_predicted)),
                                            test_predicted = test_predicted)

    # add the predicted values for the current element in the testing set to the overall data frame
    test_predicted_df <- rbind(test_predicted_df, test_predicted_element_df)
  }

  # print the predicted values for each element in the training and test sets
  # Extract unique category names
  #train
  colnames (train_predicted_df) <- c('Category', 'Column2')
  category_names <- unique(train_predicted_df$Category)
  # Create a list of data frames split by category
  train_predicted_df_list <- lapply(category_names, function(cat) {
    train_predicted_df[train_predicted_df$Category == cat, "Column2", drop = FALSE]
  })

  # Rename each data frame with its corresponding category name
  names(train_predicted_df_list) <- category_names

  # Create a named list of data frames
  train_predicted_df_named_list <- as.list(setNames(train_predicted_df_list, category_names))

  # Assign each data frame in the named list to a separate variable
  for (i in seq_along(train_predicted_df_named_list)) {
    assign(paste0("train_predicted_df", i), data.frame(train_predicted_df_named_list[[i]]))
  }
  extra_df <- data.frame(Column2=pred_V)
  train_predicted_df_list <- c(train_predicted_df_list, list(extra_df))
  names(train_predicted_df_list)[length(train_predicted_df_list)] <- "pred_V"
  train_predicted_df_named_list <- as.list(setNames(train_predicted_df_list, c(category_names,'pred_V')))

  # test
  colnames (test_predicted_df) <- c('Category', 'Column2')
  category_names <- unique(test_predicted_df$Category)
  # Create a list of data frames split by category
  test_predicted_df_list <- lapply(category_names, function(cat) {
    test_predicted_df[test_predicted_df$Category == cat, "Column2", drop = FALSE]
  })

  # Rename each data frame with its corresponding category name
  names(test_predicted_df_list) <- category_names

  # Create a named list of data frames
  test_predicted_df_named_list <- as.list(setNames(test_predicted_df_list, category_names))
  # Assign each data frame in the named list to a separate variable
  for (i in seq_along(test_predicted_df_named_list)) {
    assign(paste0("test_predicted_df", i), data.frame(test_predicted_df_named_list[[i]]))
  }
  extra_df <- data.frame(Column2=forecast_V)
  test_predicted_df_list <- c(test_predicted_df_list, list(extra_df))
  names(test_predicted_df_list)[length(test_predicted_df_list)] <- "forecast_V"
  test_predicted_df_named_list <- as.list(setNames(test_predicted_df_list, c(category_names,'forecast_V')))

  # pso_train
  y_actual <- y_train
  df_list <- train_predicted_df_named_list
  # Function to optimize weights using PSO
  optimize_weights <- function(df_list, y_actual) {
    # Flatten the list into a matrix
    df_mat <- as.matrix(do.call(cbind, df_list))

    # Define the objective function for PSO
    objective <- function(weights) {
      # Convert weights to numeric vector
      weights <- as.numeric(weights)
      # Calculate weighted sum of predicted values
      y_pred <- df_mat %*% weights

      # Calculate mean squared error
      mse <- mean((y_pred - y_actual)^2)

      # Return the mean squared error
      return(mse)
    }
    # Initialize PSO solver

    solver <- psoptim(
      lower = rep(0, ncol(df_mat)),
      upper = rep(1, ncol(df_mat)),
      par = runif(ncol(df_mat), 0, 1),
      fn = objective
    )

    # Extract the optimized weights
    weights <- solver$par

    # Calculate optimized prediction
    y_pred <- df_mat %*% weights

    # Return the optimized prediction
    return(y_pred)
  }

  # Call the function to obtain optimized prediction
  y_optimized_train <- optimize_weights(df_list, y_actual)

  # pso_test
  y_actual <- y_test
  df_list <- test_predicted_df_named_list

  y_optimized_test <- optimize_weights(df_list, y_actual)

  final_y_optimized_train <- exp(y_optimized_train)*Y_train_1
  final_y_optimized_test <- exp(y_optimized_test)*Y_test_1
  # all metrics
  all_metrics <- function(actual, predicted) {
    # Calculate the residuals
    residuals <- actual - predicted
    abs_residuals <- abs(actual - predicted)
    scaled_abs_residuals <- abs_residuals/actual
    lag_frame <- data.frame(embed(actual,2))
    diff <- lag_frame[,1]-lag_frame[,2]
    abs_diff <- abs(diff)
    # Calculate simple metrics
    mse <- mean(residuals^2)
    rmse <- sqrt(mse)
    rrmse <- 100*rmse/mean(actual)
    mae <- mean(abs_residuals)
    mape <- 100*mean(scaled_abs_residuals)
    mase <- mae/mean(abs_diff)
    # calculate complex matrics
    nse <- 1- (mse/(mean(actual^2)-(mean(actual))^2))
    wi <- 1- (mse/mean((abs(actual-mean(actual))+abs(predicted-mean(actual)))^2))
    lme <- 1- mae/mean(abs(actual-mean(actual)))
    # creating the data frame
    AllMetrics <- data.frame(cbind(c('RMSE', 'RRMSE',
                                     'MAE', 'MAPE',
                                     'MASE','NSE',
                                     'WI', 'LME'),
                                   c(round(rmse,3), round(rrmse,3),
                                     round(mae,3), round(mape,3),
                                     round(mase,3), round(nse,3),
                                     round(wi,3), round(lme,3))))
    colnames(AllMetrics) <- c('Metrics','Values')
    dimnames(AllMetrics)
    dim(AllMetrics)
    # returning the table containing all the metrics
    return(AllMetrics)
  }

  # validation_warigaan
  metrics_warigaan_train <- data.frame(all_metrics( Y_train, final_y_optimized_train))
  metrics_warigaan_test <- data.frame(all_metrics( Y_test, final_y_optimized_test))

  metrics <- cbind(as.numeric(metrics_warigaan_train$Values),
                   as.numeric(metrics_warigaan_test$Values))
  colnames(metrics) <- c('WARIGAAN_Train',
                         'WARIGAAN_Test')
  row.names(metrics)<-c('RMSE', 'RRMSE',
                        'MAE', 'MAPE',
                        'MASE','NSE',
                        'WI', 'LME')
  predict_compare <- data.frame(cbind(train_actual = Y_train,
                                      predicted = final_y_optimized_train))
  colnames(predict_compare) <- c('train_actual', 'train_predicted')
  forecast_compare <- data.frame(cbind(test_actual = Y_test,
                                       forecast = final_y_optimized_test))
  colnames(forecast_compare) <- c('test_actual', 'test_predicted')

  return(list(Train_fitted=predict_compare,
              Test_predicted=forecast_compare,
              Accuracy=metrics))
}

#' @title Wavelet Decomposition-Based ARIMA-GARCH-SVR Hybrid Modeling
#'
#' @param Y Univariate time series
#' @param ratio Ratio of number of observations in training and testing sets
#' @param n_lag  Lag of the provided time series data
#' @param l  Level of decomposition
#' @param f  Filter of decomposition
#' @import stats wavelets, tseries, forecast, fGarch, aTSA, FinTS, LSTS, earth, caret, neuralnet, pso
#' @return
#' \itemize{
#'   \item Train_fitted: Train fitted result
#'   \item Test_predicted: Test predicted result
#'   \item Accuracy: Accuracy
#'   }
#' @export
#'
#' @examples
#' Y <- rnorm(100, 100, 10)
#' result <- warigas(Y, ratio = 0.8, n_lag = 4)
#' @references
#' \itemize{
#'   \item Paul, R. K., & Garai, S. (2021). Performance comparison of wavelets-based machine learning technique for forecasting agricultural commodity prices. Soft Computing, 25(20), 12857-12873.
#'   \item Paul, R. K., & Garai, S. (2022). Wavelets based artificial neural network technique for forecasting agricultural prices. Journal of the Indian Society for Probability and Statistics, 23(1), 47-61.
#'   \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
#'   }
#'
warigas <- function(Y, ratio = 0.9, n_lag = 4, l = 6, f = 'haar'){
  optimize_weights<-NULL
  objective<-NULL
  all_metrics<-NULL
  # some default values have been set
  # Y is a vector
  # ratio is train:test
  # n_lag is number of lags in the data
  # Wavelet parameters #### (Level = l, Filter = f)
  # Level <- c(2:floor(log2(length(Y))))
  # Filter<-c('haar', 'bl14','la8', 'c6')
  ######################################
  # Data preparation ####

  # embedding for finding log return
  diff.1 <- embed(Y, 2)
  Yt <- diff.1[,1]
  Yt_1 <- diff.1[,2]
  y <- log(Yt/Yt_1)
  # Compute the average of non-zero contents in the data
  nonzero_data <- y[y != 0 & !is.na(y)]
  average_nonzero <- mean(nonzero_data)

  # Replace NaNs with the average of non-zero contents in the data
  y[is.nan(y)] <- average_nonzero

  # Check the result
  y
  # embedding for finding lag series of actual series
  # n_lag <- 4
  embed_size_y <- n_lag+1 # same lag (n_lag_y-1) for every sub series
  diff.2 <- embed(Yt,embed_size_y)
  dim(diff.2)
  Y_actual <- diff.2[,1]
  Y_actual_1 <- diff.2[,2]
  # train-test split
  # ratio <- 0.8
  n <- length(Y_actual) # this is already (original length-embed_y)
  Y_train <- Y_actual[1:(n*ratio)]
  Y_test <- Y_actual[(n*ratio+1):n]
  Y_train_1 <- Y_actual_1[1:(n*ratio)]
  Y_test_1 <- Y_actual_1[(n*ratio+1):n]
  # embedding for finding lag series of log return series
  diff.3 <- embed(y, embed_size_y)
  y_actual <- diff.3[,1]
  y_train <- y_actual[1:(n*ratio)]
  y_test <- y_actual[(n*ratio+1):n]

  # Data pre-processing ####

  # wavelet
  wavelet_y <- modwt(y, filter=f,n.levels=l)
  wavelet_df <- w<-cbind(as.data.frame(wavelet_y@W),V=as.data.frame(wavelet_y@V)[,l])
  # lag matrices for ceemdan decomposed series
  wavelet_lag_matrices <- list()
  for (col_name in colnames(wavelet_df)) {
    embedded <- embed(wavelet_df[[col_name]], embed_size_y)
    wavelet_lag_matrices_name <- paste0("embed_", col_name)
    wavelet_lag_matrices[[wavelet_lag_matrices_name]] <- embedded
  }
  wavelet_lag_train <- list()
  for (names in names(wavelet_lag_matrices)) {
    wavelet_train <- wavelet_lag_matrices[[names]][1:(n*ratio),]
    wavelet_lag_train_name <- names
    wavelet_lag_train[[wavelet_lag_train_name]] <- wavelet_train
  }

  wavelet_lag_test <- list()
  for (names in names(wavelet_lag_matrices)) {
    wavelet_test <- wavelet_lag_matrices[[names]][(n*ratio+1):n,]
    wavelet_lag_test_name <- names
    wavelet_lag_test[[wavelet_lag_test_name]] <- wavelet_test
  }

  # model fitting ####

  # ARIMA
  arima <- auto.arima(wavelet_lag_train$embed_V[,1])
  order <- arimaorder(arima)
  model_arima<-arima(wavelet_train[,ncol(wavelet_train)],order=c(order[1], order[2], order[3]))
  pred_arima <-arima$fitted
  forecast_arima <- data.frame(predict(arima,n.ahead=((n-(n*ratio)))))
  forecast_arima <-forecast_arima$pred

  #GARCH Model ####
  ARCH_pvalue <- as.numeric(FinTS::ArchTest(model_arima$residuals)$p.value)
  #ARCH_pvalue<-1
  if(ARCH_pvalue<=0.05){
    garch.fit <- fGarch::garchFit(~garch(1, 1), data = y_train, trace = FALSE)
    pred_V <- garch.fit@fitted
    forecast_V <- predict(garch.fit, n.ahead=(n-(n*ratio)))
    forecast_V <- forecast_V$meanForecast

    Resid_V <- garch.fit@residuals
    for_resid<-as.ts(y_test-forecast_V)
  }else {
    pred_V <- pred_arima
    forecast_V <- forecast_arima

    Resid_V <- as.ts(model_arima$residuals)
    for_resid<-as.vector(y_test-as.vector(forecast_arima))
  }
  names(for_resid)<-"Resid"
  names(Resid_V)<-"Resid"
  bl.test <- Box.Ljung.Test(Resid_V) #not white noise
  if(max(bl.test$data$y)<=0.05){
    Resid_ind<-1
  } else {
    Resid_ind<-2
  }

  # mars_data preparation
  mars_train <- list()
  for (names in names(wavelet_lag_train)) {
    wavelet_train_data <- cbind(wavelet_lag_train[[names]][,1])
    wavelet_train_data_name <- names
    mars_train[[wavelet_train_data_name]] <- wavelet_train_data
  }
  mars_test <- list()
  for (names in names(wavelet_lag_test)) {
    wavelet_test_data <- cbind(wavelet_lag_test[[names]][,1])
    wavelet_test_data_name <- names
    mars_test[[wavelet_test_data_name]] <- wavelet_test_data
  }
  # if arima/garch residual is not white noise it will be used as mars predictor
  # if(Resid_ind==2){
  #   predictors_train <- cbind(mars_train[-length(mars_train)],Resid=Resid_V)
  #   colnames(predictors_train)<-c(colnames(mars_train[-length(mars_train)]), "Resid")
  #   predictors_test <- cbind(mars_test[-length(mars_test)],Resid=for_resid)
  #   colnames(predictors_test)<-c(colnames(mars_test[-length(mars_test)]), "Resid")
  # } else {
  predictors_train<-data.frame(do.call(cbind, head(mars_train, -1)))
  colnames(predictors_train)<-c(names(mars_train[-length(mars_train)]))
  predictors_test<-data.frame(do.call(cbind, head(mars_test, -1)))
  colnames(predictors_test)<-c(names(mars_train[-length(mars_train)]))
  # }
  # mars
  data_MARS<-data.frame(cbind(y=y_train,predictors_train))
  names(data_MARS)
  MARS <- earth(y~., data = data_MARS)
  selcted_var<-MARS$namesx
  # train and test data for ann and svr
  data_train <- wavelet_lag_train[selcted_var]
  data_test <- wavelet_lag_test[selcted_var]
  for(i in seq_along(data_train)) {
    # determine the name of the current element
    list_name <- names(data_train)[i]

    # determine the number of columns for the current matrix
    num_cols <- ncol(data_train[[i]])

    # create a vector of column names for the current matrix
    col_names <- paste(list_name, "col", 1:num_cols, sep = "_")

    # assign the column names to the current matrix
    colnames(data_train[[i]]) <- col_names
  }
  for(i in seq_along(data_test)) {
    # determine the name of the current element
    list_name <- names(data_test)[i]

    # determine the number of columns for the current matrix
    num_cols <- ncol(data_test[[i]])

    # create a vector of column names for the current matrix
    col_names <- paste(list_name, "col", 1:num_cols, sep = "_")

    # assign the column names to the current matrix
    colnames(data_test[[i]]) <- col_names
  }

  # svr
  # create empty data frames to store predicted values for each element's response
  train_predicted_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(train_predicted_df) <- c("element", "train_predicted")
  test_predicted_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(test_predicted_df) <- c("element", "test_predicted")
  for (i in seq_along(data_train)) {
    # extract the predictors and response from the current matrix
    predictors <- data_train[[i]][, -1]
    response <- data_train[[i]][, 1]
    allvars <- colnames(predictors)
    predictorvars <- paste(allvars, collapse="+")
    form <- as.formula(paste0(names(data_train)[i], "_col_1~", predictorvars))

    # train a svm model with radial basis function as kernel
    svmmodel <- svm(form, data_train[[i]], type="eps-regression", kernel = 'radial')

    # predict the response for the current matrix in the training set
    train_predicted <- data.frame(predict(svmmodel, predictors))

    # create a data frame with the element name and predicted values for the training set
    train_predicted_element_df <- data.frame(element = rep(names(data_train[i]), nrow(train_predicted)),
                                             train_predicted = train_predicted)

    # add the predicted values for the current element in the training set to the overall data frame
    train_predicted_df <- rbind(train_predicted_df, train_predicted_element_df)
    # save the model for this element
    assign(paste0("svmmodel_", names(data_train[i])), svmmodel)
  }

  # create empty data frames to store predicted values for each element's response
  test_predicted_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(test_predicted_df) <- c("element", "test_predicted")
  test_predicted_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(test_predicted_df) <- c("element", "test_predicted")
  for (i in seq_along(data_test)) {
    # extract the predictors and response from the current matrix
    predictors <- data_test[[i]][, -1]
    response <- data_test[[i]][, 1]

    # predict the response for the current matrix in the testing set
    svmmodel <- get(paste0("svmmodel_", names(data_test[i])))
    test_predicted <- data.frame(predict(svmmodel, predictors))

    # create a data frame with the element name and predicted values for the testing set
    test_predicted_element_df <- data.frame(element = rep(names(data_test[i]), nrow(test_predicted)),
                                            test_predicted = test_predicted)

    # add the predicted values for the current element in the testing set to the overall data frame
    test_predicted_df <- rbind(test_predicted_df, test_predicted_element_df)
  }

  # print the predicted values for each element in the training and test sets
  # Extract unique category names
  #train
  colnames (train_predicted_df) <- c('Category', 'Column2')
  category_names <- unique(train_predicted_df$Category)
  # Create a list of data frames split by category
  train_predicted_df_list <- lapply(category_names, function(cat) {
    train_predicted_df[train_predicted_df$Category == cat, "Column2", drop = FALSE]
  })

  # Rename each data frame with its corresponding category name
  names(train_predicted_df_list) <- category_names

  # Create a named list of data frames
  train_predicted_df_named_list <- as.list(setNames(train_predicted_df_list, category_names))

  # Assign each data frame in the named list to a separate variable
  for (i in seq_along(train_predicted_df_named_list)) {
    assign(paste0("train_predicted_df", i), data.frame(train_predicted_df_named_list[[i]]))
  }
  extra_df <- data.frame(Column2=pred_V)
  train_predicted_df_list <- c(train_predicted_df_list, list(extra_df))
  names(train_predicted_df_list)[length(train_predicted_df_list)] <- "pred_V"
  train_predicted_df_named_list <- as.list(setNames(train_predicted_df_list, c(category_names,'pred_V')))

  # test
  colnames (test_predicted_df) <- c('Category', 'Column2')
  category_names <- unique(test_predicted_df$Category)
  # Create a list of data frames split by category
  test_predicted_df_list <- lapply(category_names, function(cat) {
    test_predicted_df[test_predicted_df$Category == cat, "Column2", drop = FALSE]
  })

  # Rename each data frame with its corresponding category name
  names(test_predicted_df_list) <- category_names

  # Create a named list of data frames
  test_predicted_df_named_list <- as.list(setNames(test_predicted_df_list, category_names))
  # Assign each data frame in the named list to a separate variable
  for (i in seq_along(test_predicted_df_named_list)) {
    assign(paste0("test_predicted_df", i), data.frame(test_predicted_df_named_list[[i]]))
  }
  extra_df <- data.frame(Column2=forecast_V)
  test_predicted_df_list <- c(test_predicted_df_list, list(extra_df))
  names(test_predicted_df_list)[length(test_predicted_df_list)] <- "forecast_V"
  test_predicted_df_named_list <- as.list(setNames(test_predicted_df_list, c(category_names,'forecast_V')))

  # pso_train
  y_actual <- y_train
  df_list <- train_predicted_df_named_list
  # Function to optimize weights using PSO
  optimize_weights <- function(df_list, y_actual) {
    # Flatten the list into a matrix
    df_mat <- as.matrix(do.call(cbind, df_list))

    # Define the objective function for PSO
    objective <- function(weights) {
      # Convert weights to numeric vector
      weights <- as.numeric(weights)
      # Calculate weighted sum of predicted values
      y_pred <- df_mat %*% weights

      # Calculate mean squared error
      mse <- mean((y_pred - y_actual)^2)

      # Return the mean squared error
      return(mse)
    }
    # Initialize PSO solver

    solver <- psoptim(
      lower = rep(0, ncol(df_mat)),
      upper = rep(1, ncol(df_mat)),
      par = runif(ncol(df_mat), 0, 1),
      fn = objective
    )

    # Extract the optimized weights
    weights <- solver$par

    # Calculate optimized prediction
    y_pred <- df_mat %*% weights

    # Return the optimized prediction
    return(y_pred)
  }

  # Call the function to obtain optimized prediction
  y_optimized_train <- optimize_weights(df_list, y_actual)

  # pso_test
  y_actual <- y_test
  df_list <- test_predicted_df_named_list

  y_optimized_test <- optimize_weights(df_list, y_actual)

  final_y_optimized_train <- exp(y_optimized_train)*Y_train_1
  final_y_optimized_test <- exp(y_optimized_test)*Y_test_1
  # all metrics
  all_metrics <- function(actual, predicted) {
    # Calculate the residuals
    residuals <- actual - predicted
    abs_residuals <- abs(actual - predicted)
    scaled_abs_residuals <- abs_residuals/actual
    lag_frame <- data.frame(embed(actual,2))
    diff <- lag_frame[,1]-lag_frame[,2]
    abs_diff <- abs(diff)
    # Calculate simple metrics
    mse <- mean(residuals^2)
    rmse <- sqrt(mse)
    rrmse <- 100*rmse/mean(actual)
    mae <- mean(abs_residuals)
    mape <- 100*mean(scaled_abs_residuals)
    mase <- mae/mean(abs_diff)
    # calculate complex matrics
    nse <- 1- (mse/(mean(actual^2)-(mean(actual))^2))
    wi <- 1- (mse/mean((abs(actual-mean(actual))+abs(predicted-mean(actual)))^2))
    lme <- 1- mae/mean(abs(actual-mean(actual)))
    # creating the data frame
    AllMetrics <- data.frame(cbind(c('RMSE', 'RRMSE',
                                     'MAE', 'MAPE',
                                     'MASE','NSE',
                                     'WI', 'LME'),
                                   c(round(rmse,3), round(rrmse,3),
                                     round(mae,3), round(mape,3),
                                     round(mase,3), round(nse,3),
                                     round(wi,3), round(lme,3))))
    colnames(AllMetrics) <- c('Metrics','Values')
    dimnames(AllMetrics)
    dim(AllMetrics)
    # returning the table containing all the metrics
    return(AllMetrics)
  }

  # validation_warigas
  metrics_warigas_train <- data.frame(all_metrics( Y_train, final_y_optimized_train))
  metrics_warigas_test <- data.frame(all_metrics( Y_test, final_y_optimized_test))

  metrics <- cbind(as.numeric(metrics_warigas_train$Values),
                   as.numeric(metrics_warigas_test$Values))
  colnames(metrics) <- c('WARIGAS_Train',
                         'WARIGAS_Test')
  row.names(metrics)<-c('RMSE', 'RRMSE',
                        'MAE', 'MAPE',
                        'MASE','NSE',
                        'WI', 'LME')
  predict_compare <- data.frame(cbind(train_actual = Y_train,
                                      predicted = final_y_optimized_train))
  colnames(predict_compare) <- c('train_actual', 'train_predicted')
  forecast_compare <- data.frame(cbind(test_actual = Y_test,
                                       forecast = final_y_optimized_test))
  colnames(forecast_compare) <- c('test_actual', 'test_predicted')

  return(list(Train_fitted=predict_compare,
              Test_predicted=forecast_compare,
              Accuracy=metrics))
}
