#' Performs a Inverse Turkey power tranformation with predefined Lambda
#'
#' Takes in a vector and a lambda value to perform Inverse Turkey power transfomation
#' @param x A vector which needs to be inverse-tranformed
#' @param lambda_fit A fixed lambda value for turkey power transformation
#' @return Returns exponential tranformattion if lambda equals zero and a inverse power transormation otherwise
#' @export
#' @usage InvTurkey(x, lambda_fit)
#'
#' @example
#' a <- c(12,34,234,23,678, 768, 34, 34 ,78)
#' InvTurkey(a, -5)
#'
InvTurkey <- function(x, lambda_fit)
{
  lambda <- lambda_fit$turkey_lamda
  if (lambda > 0) {
    TRANS = x^(1/lambda)
  }
  if (lambda == 0) {
    TRANS = exp(x)
  }
  if (lambda < 0) {
    TRANS = -1 * x^(1/lambda)
  }
  return(TRANS)
}

#' Fit function that performs inverse transformations to the dataset
#'
#' Takes in transformed dataset and inverse-transform it based on the model file provided
#' @param data Transformed dataset
#' @param trans_fit_model A model file that captures the details of transformations done to given data. This is returned as a list element by VariableTransform function
#' @return An inverse transformed dataset
#' @export
#' @usage InvTransform(data, trans_fit_model)
InvTransform <- function(data, trans_fit_model)
{
  # require(data.table)
  # require(dplyr)

  if(is.null(trans_fit_model))
    return(data)

  chosen_trans <- data.frame(t(data.frame(trans_fit_model$chosen_transforms)))
  colnames(chosen_trans) <- c("Chosen_Trans")
  chosen_trans$column_names <- row.names(chosen_trans)

  data_org <- data
  ## Transformations
  no_trans_list<-chosen_trans[chosen_trans$Chosen_Trans == "no_transform",]$column_names
  arcsin_list<-chosen_trans[chosen_trans$Chosen_Trans == "arcsin_transformed",]$column_names
  log_list<-chosen_trans[chosen_trans$Chosen_Trans == "log_transformed",]$column_names
  sqrt_list<-chosen_trans[chosen_trans$Chosen_Trans == "sqrt_transformed",]$column_names
  cu_rt_list<-chosen_trans[chosen_trans$Chosen_Trans == "cube_root_transformed",]$column_names
  turkey_list<-chosen_trans[chosen_trans$Chosen_Trans == "turkey_transformed",]$column_names
  boxcox_list<-chosen_trans[chosen_trans$Chosen_Trans == "boxcox_transformed",]$column_names
  sqr_list<-chosen_trans[chosen_trans$Chosen_Trans == "sqr_transformed",]$column_names
  cube_list<-chosen_trans[chosen_trans$Chosen_Trans == "cube_transformed",]$column_names
  yeo_johnson_list <- chosen_trans[chosen_trans$Chosen_Trans == "yeo_johnson_transformed",]$column_names
  lambertw_list <- chosen_trans[chosen_trans$Chosen_Trans == "lambertw_transformed",]$column_names
  order_norm_list <- chosen_trans[chosen_trans$Chosen_Trans == "order_norm_transformed",]$column_names

  ## Reverse Transormation
  no_trans<-data.table(data)[,lapply(.SD,function(x){x}),.SDcols=no_trans_list]
  arcsin_trans<-data.table(data)[,lapply(.SD,function(x){sinh(x)}),.SDcols=arcsin_list]
  log_trans<-data.table(data)[,lapply(.SD,function(x){exp(x)}),.SDcols=log_list]
  sqrt_trans<-data.table(data)[,lapply(.SD,function(x){x^2}),.SDcols=sqrt_list]
  cube_root_trans<-data.table(data)[,lapply(.SD,function(x){sign(x)*abs(x)^(3)}),.SDcols=cu_rt_list]
  sqr_trans<-data.table(data)[,lapply(.SD,function(x){(sqrt(x))}),.SDcols=sqr_list]
  cube_trans<-data.table(data)[,lapply(.SD,function(x){sign(x)*abs(x)^(1/3)}),.SDcols=cube_list]

  data <- data.frame(data)
  turkey_trans<- data.frame(mapply(InvTurkey, data[,turkey_list],trans_fit_model$possible_fits[turkey_list]))
  colnames(turkey_trans) <- turkey_list

  boxcox_trans<- data.frame(mapply(function(x,y){forecast::InvBoxCox(x,y$boxcox_lamda)}, data[,boxcox_list],trans_fit_model$possible_fits[boxcox_list]))
  colnames(boxcox_trans) <- boxcox_list



  ifelse(length(yeo_johnson_list) > 0, yeo_johnson_trans<- data.frame(mapply(function(x, y) {yeo.johnson(x, y$yeo_johnson_lamda, inverse=TRUE)}, data[,yeo_johnson_list], trans_fit_model$possible_fits[yeo_johnson_list])), yeo_johnson_trans <- data.frame())
  # yeo_johnson_trans<- data.frame(mapply(yeo.johnson, data[,yeo_johnson_list],lambda=data.table(trans_fit)[chosen_method=="yeo_johnson_transformed",]$yeo_johnson_lamda, inverse=TRUE))
  colnames(yeo_johnson_trans) <- yeo_johnson_list

  # lambertw_trans<- data.frame(mapply(Gaussianize, data[,lambertw_list], tau_mat = tau_mat_fit[lambertw_list]))
  for (i in 1:length(lambertw_list))
  {
    # print(i)
    if (length(lambertw_list) == 0)
    {
      lambertw_trans <-  data.frame()
      break
    }
    tau_mat_fit <- trans_fit_model$possible_fits[lambertw_list]
    if(!is.null(tau_mat_fit[[lambertw_list[i]]]$lambertw_tau_mat))
    {
      if (i == 1)
      {
        lambertw_trans <- Gaussianize(data[,lambertw_list[i]],tau.mat = tau_mat_fit[[lambertw_list[i]]]$lambertw_tau_mat, inverse = TRUE)
      }else
      {
        lambertw_trans <- cbind(lambertw_trans,Gaussianize(data[,lambertw_list[i]],tau.mat = tau_mat_fit[[lambertw_list[i]]]$lambertw_tau_mat, inverse = TRUE))
      }
    }else
    {
      if (i == 1)
      {
        lambertw_trans <- data[,lambertw_list[i]]
      }else
      {
        lambertw_trans <- cbind(lambertw_trans,data[,lambertw_list[i]])
      }
    }

  }
  lambertw_trans <- data.frame(lambertw_trans)
  colnames(lambertw_trans) <- lambertw_list

  order_norm_trans<- data.frame(mapply(function(x, y){predict(newdata = x, object = y$order_norm_object, inverse = TRUE)}, data[,order_norm_list], trans_fit_model$possible_fits[order_norm_list]))
  colnames(order_norm_trans) <- order_norm_list

  q <- c(no_trans, arcsin_trans, lambertw_trans, order_norm_trans, turkey_trans, boxcox_trans, log_trans, sqrt_trans, cube_root_trans, sqr_trans, cube_trans, yeo_johnson_trans)
  transformed_df <- data.frame(q[order(sapply(q,length),decreasing=T)])

  # transformed_df <- bind_cols(no_trans, arcsin_trans, lambertw_trans, order_norm_trans, turkey_trans, boxcox_trans, log_trans, sqrt_trans, cube_root_trans, sqr_trans, cube_trans, yeo_johnson_trans)
  # complete_data<-data.frame(cbind(transformed_df,data.frame(data_org)[,(!colnames(data_org) %in% colnames(transformed_df))]))

  ## Reverse padding..
  # padding_list <- as.character(data.table(trans_fit)$names)
  transformed_df <- data.frame(transformed_df)
  padding_list <- names(trans_fit_model$possible_fits)
  transformed_df <- data.frame(mapply(function(x,fit_min){x-1+fit_min$min_val_padding}, transformed_df[,padding_list], trans_fit_model$possible_fits))

  # return(complete_data)
  return(transformed_df)

}
