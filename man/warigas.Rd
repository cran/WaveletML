\name{warigas}
\alias{warigas}
\title{Wavelet Decomposition-Based ARIMA-GARCH-SVR Hybrid Modeling}
\usage{
warigas(Y, ratio = 0.9, n_lag = 4, l = 6, f = 'haar')
}
\arguments{
\item{Y}{Univariate time series}

\item{ratio}{Ratio of number of observations in training and testing sets}

\item{n_lag}{Lag of the provided time series data}

\item{l}{Level of decomposition}

\item{f}{Filter of decomposition}
}
\value{
\itemize{
  \item Train_fitted: Train fitted result
  \item Test_predicted: Test predicted result
  \item Accuracy: Accuracy
  }
}
\description{
Wavelet Decomposition-Based ARIMA-GARCH-SVR Hybrid Modeling
}
\examples{
Y <- rnorm(100, 100, 10)
result <- warigas(Y, ratio = 0.8, n_lag = 4)
}
\references{
\itemize{
  \item Paul, R. K., & Garai, S. (2021). Performance comparison of wavelets-based machine learning technique for forecasting agricultural commodity prices. Soft Computing, 25(20), 12857-12873.
  \item Paul, R. K., & Garai, S. (2022). Wavelets based artificial neural network technique for forecasting agricultural prices. Journal of the Indian Society for Probability and Statistics, 23(1), 47-61.
  \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
  }
}
