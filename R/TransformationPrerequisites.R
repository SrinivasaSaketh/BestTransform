#' Performs Box-cox Transformation on a vector
#'
#' Takes in a vector and transforms it using a Box-cox Power Trasformation
#' @param one.col Any numeric vector that needs to be transformed
#' @param returnLambda A boolean value which returns the optimal lambda value for a vector when se to TRUE
#' @return Returns the optimal lambda value for Box-cox transformation when "returnLambda" is set to TRUE and a Box-cox transformed vector otherwise.
#' @export
#' @usage boxcoxtransForecast(one.col, returnLambda = FALSE)
#' @example
#' a <- c(12,34,234,23,678, 768, 34, 34 ,78)
#' boxcoxtransForecast(a)
#'

boxcoxtransForecast <- function(one.col, returnLambda = FALSE)
{
  # require(forecast)
  if (returnLambda == TRUE)
    return(BoxCox.lambda(one.col))
  return(BoxCox(one.col, BoxCox.lambda(one.col)))
}

#' Performs Yeo-Johnson Transformation on a vector
#'
#' Takes in a vector and transforms it using a Yeo-Johnson Power Trasformation
#' @param one.col Any numeric vector that needs to be transformed
#' @param returnLambda A boolean value which returns the optimal lambda value for a vector when se to TRUE
#' @return Returns the optimal lambda value for Yeo-Johnson transformation when "returnLambda" is set to TRUE and a Yeo-Johnson transformed vector otherwise.
#' @export
#' @usage yeoJohnsonVGAM(one.col, returnLambda = FALSE)
#' @example
#' a <- c(12,34,234,23,678, 768, 34, 34 ,78)
#' yeoJohnsonVGAM(a)
#'

yeoJohnsonVGAM <- function(one.col, returnLambda = FALSE)
{
  # require(forecast)
  # require(VGAM)
  if (returnLambda == TRUE)
    return(BoxCox.lambda(one.col))
  return(yeo.johnson(one.col, BoxCox.lambda(one.col)))
}

#' Performs Lambert Transformation on a vector
#'
#' Takes in a vector and transforms it using a Lmabert-W Transformation
#' @param one.col Any numeric vector that needs to be transformed
#' @param return_tau_mat A boolean value which returns the tau_mat vector when set to TRUE
#' @return Returns the tau_mat vector for Lambert-W transformation when "return_tau_mat" is set to TRUE and a Lambert-W transformed vector otherwise.
#' @export
#' @usage GaussianizeLambertW(one.col, return_tau_mat = FALSE)
#' @example
#' a <- c(12,34,234,23,678, 768, 34, 34 ,78)
#' GaussianizeLambertW(a, return_tau_mat = FALSE)
#'
GaussianizeLambertW <- function(one.col, return_tau_mat = FALSE)
{
  # require(LambertW)
  if ((sort(table(one.col), decreasing = T)[1]/length(one.col) * 100) < 50)
  {
    data.col <- data.frame(Gaussianize(one.col))
    tau_mat <- Gaussianize(one.col, return.tau.mat = TRUE)$tau.mat
  } else
  {
    data.col <- one.col
    tau_mat <- NULL
  }
  if (return_tau_mat == TRUE)
    return(tau_mat)
  return(data.col)
}

#' Performsturkey power transformation on a vector
#'
#' Takes in a vector and conducts Tukey's Ladder of Powers on a vector of values to produce a more-normally distributed vector of values
#' @param x A vector of values
#' @param start The starting value of lambda to try
#' @param end The ending value of lambda to try
#' @param int The interval between lambda values to try
#' @param plotit If TRUE, produces plots of Shapiro-Wilks W or Anderson-Darling A vs. lambda, a histogram of transformed values, and a quantile-quantile plot of transformed values
#' @param verbose If TRUE, prints extra output for Shapiro-Wilks W or Anderson-Darling A vs. lambda
#' @param quiet If TRUE, doesn't print any output to the screen
#' @param statistic If 1, uses Shapiro-Wilks test. If 2, uses Anderson-Darling test
#' @param returnLambda If TRUE, returns only the lambda value, not the vector of transformed values
#'
#' @return The transformed vector of values.
#' @export
#' @usage transformTurkey(x, start, end, int, plotit, verbose, quiet, statistic, returnLambda)
transformTurkey <- function (x, start = -10, end = 10, int = 0.025, plotit = TRUE,
                             verbose = FALSE, quiet = FALSE, statistic = 1, returnLambda = FALSE)
{
  n = (end - start)/int
  lambda = as.numeric(rep(0, n))
  W = as.numeric(rep(0, n))
  Shapiro.p.value = as.numeric(rep(0, n))
  if (statistic == 2) {
    A = as.numeric(rep(1000, n))
    Anderson.p.value = as.numeric(rep(0, n))
  }
  for (i in (1:n)) {
    lambda[i] = signif(start + (i - 1) * int, digits = 4)
    if (lambda[i] > 0) {
      TRANS = x^lambda[i]
    }
    if (lambda[i] == 0) {
      TRANS = log(x)
    }
    if (lambda[i] < 0) {
      TRANS = -1 * x^lambda[i]
    }
    W[i] = NA
    if (statistic == 2) {
      A[i] = NA
    }
    ifelse(length(TRANS) > 5000, TRANS <- sample(TRANS, 4999), TRANS)
    if (any(is.infinite(TRANS)) == FALSE & any(is.nan(TRANS)) ==
        FALSE) {
      W[i] = signif(shapiro.test(TRANS)$statistic, digits = 4)
      Shapiro.p.value[i] = signif(shapiro.test(TRANS)$p.value,
                                  digits = 4)
      if (statistic == 2) {
        A[i] = signif(ad.test(TRANS)$statistic, digits = 4)
        Anderson.p.value[i] = signif(ad.test(TRANS)$p.value,
                                     digits = 4)
      }
    }
  }
  if (statistic == 2) {
    df = data.frame(lambda, W, Shapiro.p.value, A, Anderson.p.value)
  }
  if (statistic == 1) {
    df = data.frame(lambda, W, Shapiro.p.value)
  }
  if (verbose == TRUE) {
    print(df)
  }
  if (plotit == TRUE) {
    if (statistic == 1) {
      plot(lambda, W, col = "black")
    }
    if (statistic == 2) {
      plot(lambda, A, col = "blue")
    }
  }
  if (statistic == 1) {
    df2 = df[with(df, order(-W)), ]
  }
  if (statistic == 2) {
    df2 = df[with(df, order(A)), ]
  }
  if (quiet == FALSE) {
    cat("\n")
    print(df2[1, ])
    cat("\n")
    cat("if (lambda >  0){TRANS = x ^ lambda}", "\n")
    cat("if (lambda == 0){TRANS = log(x)}", "\n")
    cat("if (lambda <  0){TRANS = -1 * x ^ lambda}", "\n")
    cat("\n")
  }
  lambda = df2[1, "lambda"]
  if (lambda > 0) {
    TRANS = x^lambda
  }
  if (lambda == 0) {
    TRANS = log(x)
  }
  if (lambda < 0) {
    TRANS = -1 * x^lambda
  }
  if (plotit == TRUE) {
    plotNormalHistogram(TRANS, xlab = "Transformed variable",
                        linecol = "red", col = "lightgray")
  }
  if (plotit == TRUE) {
    qqnorm(TRANS)
    qqline(TRANS, col = "red")
  }
  if (returnLambda == FALSE) {
    return(TRANS)
  }
  if (returnLambda == TRUE) {
    names(lambda) = "lambda"
    return(lambda)
  }
}

#' Perform Linear model on the provided dataset
#'
#' Takes in a dataset gives the perfoormance measures by taking a 75:25 train-test split. Performance measures are considered based the distribution of the DV column
#' @param data A list of datasets
#' @param num Index of the dataset to be accessed
#' @return The performance measures of the dataset
#' @export
#' @usage model_my_data(data, num)
model_my_data<-function(data, num)
{
  # require(caret)
  # library(mlbench)
  data <- as.data.frame(data[[num]])
  set.seed(280)
  index <- createDataPartition(data$dvcol, p = .75, list = FALSE)
  df_tr <- data[index, ]
  df_te <- data[index, ]

  dist = data_distribution(data, "dvcol")
  dvtype<-as.character(dist[dist$names=="dvcol",3])


  if(dvtype=="Continous")
  {
    lm_fit <- train(dvcol ~ .,
                    data = df_tr,
                    method = "lm")
    dvcol_pred <- predict(lm_fit, df_te)
    q<-postResample(pred = dvcol_pred, obs = df_te$dvcol)
    # q<-data.frame(value=postResample(pred = dvcol_pred, obs = df_te$dvcol),names=names(postResample(pred = dvcol_pred, obs = df_te$dvcol)))

  }
  else
  {
    if(length(unique(data$dvcol)) > 2)
    {
      df_tr$dvcol<-as.factor(df_tr$dvcol)
      df_te$dvcol<-as.factor(df_te$dvcol)
      colnames(df_te)[ncol(df_te)] <- "obs"
      lm_fit <- train(dvcol~ .,
                      data = df_tr,
                      method = "LogitBoost")
      df_te$pred <- predict(lm_fit, df_te)
      q<- multiClassSummary(df_te, lev = levels(df_te$obs))
    } else
    {
      require(forcats)
      df_tr$dvcol<-fct_rev(as.factor(df_tr$dvcol))
      df_te$dvcol<-fct_rev(as.factor(df_te$dvcol))
      colnames(df_te)[ncol(df_te)] <- "obs"
      lm_fit <- train(dvcol~ .,
                      data = df_tr,
                      method = "LogitBoost")
      df_te$pred <- predict(lm_fit, df_te)
      q<- prSummary(df_te, lev = levels(df_te$obs))
    }


  }
  return (q)
}
