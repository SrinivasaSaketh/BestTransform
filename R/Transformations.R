#' Transforms the data vector using the best possible transformation technique
#'
#' Takes in a data vector and performs the best possible transformation to each of the columns in the data frame
#' @param data_col Any continuous data vector that has to transformed into normal form
#' @param best_trans_metric Considers the best transformation based on one of the three values ("Shapiro P Value", "Pearson P Value", "Min skewness")
#' @description VariableTransform initially pads the data inorder to eliminates all negative and zero values and then categorizes the data into normal, positive skewed and negative skewed based on the skewness score of each column of the given data frame
#'
#' \strong{Transformation for Positive skewed data:}
#'
#' Different transformation techiniques considered to transform positive skewed data are as follows:
#' \itemize{
#' \item Log Transformation
#' \item Square Root Transformation
#' \item Cube Root Transformation
#' \item Turkey Power Transformation
#' \item Box-Cox Power Transformation
#' \item Yeo-Johnson Transformation
#' \item Order-Norm Transformation
#' \item Lambert-W Transformation
#' }
#' Of all these techniques, Best technique is choosed based on the "ChooseBestTrans" argument provided.Default method will be Shapiro P-Value
#'
#' \strong{Transformation for Negative skewed data:}
#'
#' Different transformation techiniques considered to transform negative skewed data are as follows:
#' \itemize{
#' \item Square Transformation
#' \item Cube Transformation
#' \item Turkey Power Transformation
#' \item Box-Cox Power Transformation
#' \item Yeo-Johnson Transformation
#' \item Order-Norm Transformation
#' \item Lambert-W Transformation

#' }
#' Of all these techniques, Best technique is choosed based on the "ChooseBestTrans" argument provided.Default method will be Shapiro P-Value
#'
#' \strong{Scaling the dataset:}
#'
#' After merging all the datasets created (normally distributed data, positive skewd data, neggative skewed data), entire dataset is normalized and stored in an other dataframe.
#'
#' @return
#' Returns a list of 7 objects:
#' \describe{
#' \item{transformed_df}{Tranformed Dataset for all Continuous variables}
#' \item{scaled_df}{Scaled Dataset for all Continuous variables}
#' \item{original_dist}{Distribution of the dataset provided}
#' \item{neg_skew_trans}{Intermediate scores (for all the techniques) for Positive skewed data}
#' \item{pos_skew_trans}{Intermediate scores (for all the techniques) for Positive skewed data}
#' \item{trans_fit_model}{A list of model fit file, OrderNorm objects and tau_mat objects}
#' \item{complete_data}{Combination of transformed continuous data and categorical data}
#' }
#'
#' @export
#' @usage VariableTransform(data_col, best_trans_metric = "Shapiro P Value")
VariableTransform<-function(data_col, best_trans_metric = "Shapiro P Value")
{
  # require(rcompanion)
  # require(plyr)
  # require(data.table)
  # require(bestNormalize)

  if(!is.numeric(data_col))
    stop("data_col is not numeric")

  #Padding for Negetive & 0 Values. Also Storing the minvalue for a Reverse pad
  min_val <- min(data_col)
  data_col <- data_col+1-min_val

  no_transform<-data_col
  arcsin_transformed <- asinh(data_col)

  turkey_transformed<-transformTurkey(data_col,plotit=FALSE,quiet = T)
  turkey_lamda<-transformTurkey(data_col,plotit=FALSE,returnLambda = T,quiet = T)

  boxcox_transformed<-boxcoxtransForecast(data_col)
  boxcox_lamda<-boxcoxtransForecast(data_col,returnLambda = TRUE)

  yeo_johnson_transformed<-yeoJohnsonVGAM(data_col)
  yeo_johnson_lamda<-yeoJohnsonVGAM(data_col,returnLambda = TRUE)

  order_norm_transformed <- bestNormalize::orderNorm(data_col, warn = FALSE)$x.t
  order_norm_object <- bestNormalize::orderNorm(data_col, warn = FALSE)

  lambertw_transformed<-data.frame(GaussianizeLambertW(data_col))[,1]
  lambertw_tau_mat <- GaussianizeLambertW(data_col,return_tau_mat = TRUE)

  skewness <- e1071::skewness(data_col)
  if (skewness >= 0.3)
  {
    log_transformed<-log(data_col)
    sqrt_transformed<-sqrt(data_col)
    cube_root_transformed<-sign(data_col)*abs(data_col)^(1/3)
    all_transforms <- data.frame(no_transform, arcsin_transformed, log_transformed, sqrt_transformed, cube_root_transformed, turkey_transformed, boxcox_transformed, yeo_johnson_transformed, order_norm_transformed, lambertw_transformed)
  } else if (skewness <= -0.3)
  {
    sqr_transformed<-data_col^2
    cube_transformed<-sign(data_col)*abs(data_col)^(3)
    all_transforms <- data.frame(no_transform, arcsin_transformed, sqr_transformed, cube_transformed, turkey_transformed, boxcox_transformed, yeo_johnson_transformed, order_norm_transformed, lambertw_transformed)
  } else
  {
    all_transforms <- data.frame(no_transform, arcsin_transformed, turkey_transformed, boxcox_transformed, yeo_johnson_transformed, order_norm_transformed, lambertw_transformed)
  }

  if (best_trans_metric == "Pearson P Value")
  {
    pearson_scores <- data.frame(t(data.table(all_transforms)[,lapply(.SD, function(x) {nortest::pearson.test(x)$statistic/nortest::pearson.test(x)$df}), .SDcols=colnames(all_transforms)]))
    colnames(pearson_scores) <- c("PearsonScores")
    pearson_scores$method <- row.names(pearson_scores)
    best_transformation <- pearson_scores[pearson_scores$PearsonScores == min(pearson_scores$PearsonScores),]$method[1]
  } else if (best_trans_metric == "Min skewness")
  {
    skewness_scores <- data.frame(t(data.table(all_transforms)[,lapply(.SD, e1071::skewness), .SDcols=colnames(all_transforms)]))
    colnames(skewness_scores) <- c("SkewnessScores")
    skewness_scores$method <- row.names(skewness_scores)
    best_transformation <- skewness_scores[skewness_scores$SkewnessScores == min(skewness_scores$SkewnessScores),]$method[1]
  } else
  {
    shapiro_scores <- data.frame(t(data.table(all_transforms)[,lapply(.SD, isnormal, TRUE), .SDcols=colnames(all_transforms)]))
    colnames(shapiro_scores) <- c("ShapiroScores")
    shapiro_scores$method <- row.names(shapiro_scores)
    best_transformation <- shapiro_scores[shapiro_scores$ShapiroScores == min(shapiro_scores$ShapiroScores),]$method[1]
  }

  ##Reverse padding
  transformed_data <-  get(best_transformation)
  transformed_data <- transformed_data-1+min_val

  return (list(best_transformation = best_transformation, transformed_data = transformed_data,  possible_fits = list(turkey_lamda = turkey_lamda, boxcox_lamda = boxcox_lamda, yeo_johnson_lamda = yeo_johnson_lamda, order_norm_object = order_norm_object, lambertw_tau_mat = lambertw_tau_mat, min_val_padding = min_val)))

}

#' Transforms the dataset using the best possible transformation technique
#'
#' Takes in a data frame and performs the best possible transformation to each of the columns in the data frame
#' @param data Any data frame that has atlest one column with continuous data and that has to transformed into normal form
#' @param dv Dependent variable in the dataset.
#' @description
#' Considers the best transformed dataset based on three metrics
#' \itemize{
#' \item Shapiro P Value
#' \item Pearson P Value
#' \item Min skewness
#' }
#' Best possible transformation is applied to all the continous columns in the dataset for each of the above mentioned metrics and three datasets are collected respectively. Best dataset is chosen based on the performance metrics of the dataset on a base model. Original dataset is also considered for the performance measure to actually check if transforming the data is necessary.
#' @return Returns a list of 2 objects:
#' \describe{
#' \item{trans_data}{Tranformed Dataset for all Continuous variables}
#' \item{trans_fit_model}{Model fit file to fit the test data}
#' }
#' @export
#' @usage BestTransform(data, dv)
BestTransform <- function(data, dv)
{
  dist <- data_distribution(data, dv)
  cont_cols <- dist[(dist$distribution == "Continous" & dist$is_dv == FALSE),]$names
  dvcol <- data[,dv]
  ppv <- data.table(data)[,lapply(.SD,function(x){VariableTransform(x, "Pearson P Value")$transformed_data}),.SDcols=cont_cols]
  ppv1 <- cbind(ppv, dvcol)

  spv <- data.table(data)[,lapply(.SD,function(x){VariableTransform(x, "Shapiro P Value")$transformed_data}),.SDcols=cont_cols]
  spv1 <- cbind(spv, dvcol)

  skew_v <- data.table(data)[,lapply(.SD,function(x){VariableTransform(x, "Min skewness")$transformed_data}),.SDcols=cont_cols]
  skew_v1 <- cbind(skew_v, dvcol)

  # original <- data[,cont_cols]
  # original1 <- cbind(original, dvcol)
  original1 <- data

  total_matrix_all <- list("Pearson P Value" = ppv1, "Shapiro P Value" = spv1, "Min skewness" = skew_v1, "Original" = original1)
  # library(pbmcapply)
  output<-do.call(rbind,pbmclapply(seq(1:4), model_my_data,
                                   data = total_matrix_all,
                                   mc.cores = 1))
  if (dist[dist$is_dv == T, ]$distribution == "Continous")
  {
    output <- data.frame(output[, c("Rsquared", "MAE", "RMSE")])
    var <- c("Rsquared", "MAE", "RMSE")
  } else
  {
    output <- data.frame(output[, c("Mean_F1", "Mean_Precision", "Mean_Recall")])
    var <- c("Mean_F1", "Mean_Precision", "Mean_Recall")
  }

  # var<-nearZeroVar(output, saveMetrics = TRUE)
  # var<-rownames(var[(var$zeroVar=="FALSE"),])
  # out<-data.frame(output)[,var]
  output$Method <- names(total_matrix_all)
  best_trans_metric <- names(total_matrix_all[which(output[,1] == max(output[,1]))[1]])
  # best_trans_metric <- "Pearson P Value"  ## Delete this line after testing
  if (best_trans_metric == "Original")
  {
    trans_data <- total_matrix_all[[best_trans_metric]]
    trans_fit_model <- NULL
  } else
  {
    trans_data <- total_matrix_all[[best_trans_metric]]
    possible_fits <- list()
    chosen_transforms <- list()
    for (i in 1:length(cont_cols))
    {
      # print(i)
      res <- VariableTransform(data[,cont_cols[i]], best_trans_metric)
      chosen_transforms[[cont_cols[i]]] <- res$best_transformation
      possible_fits[[cont_cols[i]]] <- res$possible_fits
    }

    trans_fit_model <- list(chosen_transforms = chosen_transforms, possible_fits = possible_fits)
  }
  output <- output[,c("Method", var)]
  return(list(trans_data = trans_data, trans_fit_model = trans_fit_model, model_perf_metrics = output))

}

