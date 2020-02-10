#' Performs a Turkey power tranformation with predefined Lambda
#'
#' Takes in a vector and a lambda value to perform Turkey power transfomation
#' @param x A vector which neds to be tranformed
#' @param lambda_fit A fixed lambda value for turkey power transformation
#' @return Returns log tranformattion if lambda equals zero and a power transormation otherwise
#' @export
#' @usage TurkeyFit(x, lambda_fit)
#'
#' @example
#' a <- c(12,34,234,23,678, 768, 34, 34 ,78)
#' TurkeyFit(a, -5)
#'
TurkeyFit <- function(x, lambda_fit)
{
  lambda <- lambda_fit$turkey_lamda
  if (lambda > 0) {
    TRANS = x^lambda
  }
  if (lambda == 0) {
    TRANS = log(x)
  }
  if (lambda < 0) {
    TRANS = -1 * x^lambda
  }
  return(TRANS)
}

#' Fit function for test data transformations
#'
#' Takes in test data and fits the data based on the model file provided
#' @param data Test data to be tranformed
#' @param trans_fit_model A model file that captures the details of transformations done to train data. This is returned as a list element by VariableTransform function
#' @return A transformed dataset
#' @export
#' @usage TransformFit(data, trans_fit_model)
TransformFit <- function(data, trans_fit_model)
{
  # require(data.table)
  # require(dplyr)

  if(is.null(trans_fit_model))
    return(data)
  data <- data.frame(data)
  chosen_trans <- data.frame(t(data.frame(trans_fit_model$chosen_transforms)))
  colnames(chosen_trans) <- c("Chosen_Trans")
  chosen_trans$column_names <- row.names(chosen_trans)
  # trans_fit <- trans_fit_model$trans_fit
  # tau_mat_fit <- trans_fit_model$lambertw_tau_mat
  # order_norm_obj <- trans_fit_model$order_norm_obj

  data_org <- data
  ## padding..
  padding_list <- as.character(chosen_trans$column_names)
  data <- data.frame(mapply(function(x,fit_min){x+1-fit_min$min_val_padding}, data[,padding_list], trans_fit_model$possible_fits))

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

  no_trans<-data.table(data)[,lapply(.SD,function(x){x}),.SDcols=no_trans_list]
  arcsin_trans<-data.table(data)[,lapply(.SD,function(x){asinh(x)}),.SDcols=arcsin_list]
  log_trans<-data.table(data)[,lapply(.SD,function(x){log(x)}),.SDcols=log_list]
  sqrt_trans<-data.table(data)[,lapply(.SD,function(x){sqrt(x)}),.SDcols=sqrt_list]
  cube_root_trans<-data.table(data)[,lapply(.SD,function(x){sign(x)*abs(x)^(1/3)}),.SDcols=cu_rt_list]
  sqr_trans<-data.table(data)[,lapply(.SD,function(x){(x^2)}),.SDcols=sqr_list]
  cube_trans<-data.table(data)[,lapply(.SD,function(x){sign(x)*abs(x)^(3)}),.SDcols=cube_list]

  turkey_trans<- data.frame(mapply(TurkeyFit, data[,turkey_list],trans_fit_model$possible_fits[turkey_list]))
  colnames(turkey_trans) <- turkey_list

  boxcox_trans<- data.frame(mapply(function(x,y){forecast::BoxCox(x,y$boxcox_lamda)}, data[,boxcox_list],trans_fit_model$possible_fits[boxcox_list]))
  colnames(boxcox_trans) <- boxcox_list

  yeo_johnson_trans<- data.frame(mapply(function(x, y) {yeo.johnson(x, y$yeo_johnson_lamda)}, data[,yeo_johnson_list], trans_fit_model$possible_fits[yeo_johnson_list]))
  colnames(yeo_johnson_trans) <- yeo_johnson_list

  order_norm_trans<- data.frame(mapply(function(x, y){predict(newdata = x, object = y$order_norm_object)}, data[,order_norm_list], trans_fit_model$possible_fits[order_norm_list]))
  colnames(order_norm_trans) <- order_norm_list

  # lambertw_trans<- data.frame(mapply(Gaussianize, data[,lambertw_list], tau_mat = tau_mat_fit[lambertw_list]))
  for (i in 1:length(lambertw_list))
  {
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
        lambertw_trans <- Gaussianize(data[,lambertw_list[i]],tau.mat = tau_mat_fit[[lambertw_list[i]]]$lambertw_tau_mat)
      }else
      {
        lambertw_trans <- cbind(lambertw_trans,Gaussianize(data[,lambertw_list[i]],tau.mat = tau_mat_fit[[lambertw_list[i]]]$lambertw_tau_mat))
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

  q <- c(no_trans, arcsin_trans, lambertw_trans, order_norm_trans, turkey_trans, boxcox_trans, log_trans, sqrt_trans, cube_root_trans, sqr_trans, cube_trans, yeo_johnson_trans)
  transformed_df <- data.frame(q[order(sapply(q,length),decreasing=T)])

  # transformed_df <- bind_cols(no_trans, arcsin_trans, lambertw_trans, order_norm_trans, turkey_trans, boxcox_trans, log_trans, sqrt_trans, cube_root_trans, sqr_trans, cube_trans, yeo_johnson_trans)
  complete_data<-data.frame(cbind(transformed_df,data.frame(data_org)[,(!colnames(data_org) %in% colnames(transformed_df))]))

  return(complete_data)
  # return(transformed_df)
}
