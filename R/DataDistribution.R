#' Check for the normality of the data
#'
#' Takes in a data vector and checks for the normality of the data using Shapiro Test
#' @param data_vector Any data vector that has to be tested for normality
#' @param returnPValue A boolean value which when TRUE returns Shapiro P Value
#' @return
#' \itemize{
#' \item If returnPValue = TRUE, returns the Shapiro P-Value
#' \item If returnPValue = FALSE, returns TRUE if the data vector is normally dirstributed and FALSE otherwise
#' }
#'
#' @export
#' @usage isnormal(data_vector, returnPValue=FALSE)
#' @example
#' a <- c(12,34,234,23,678, 768, 34, 34 ,78)
#' isnormal(a, FALSE)
isnormal <- function(data_vector, returnPValue = FALSE)
{
  data <- data_vector
  ifelse(length(data_vector) > 5000, data_vector <- sample(data, 4999), data_vector)
  # print(length(unique(data_vector)))
  while (length(unique(data_vector)) == 1)
  {
    # print(TRUE)
    data_vector <- sample(data, 4999)
    # print(length(unique(data_vector)))
  }
  # print(length(unique(data_vector)))
  st <- shapiro.test(data_vector)
  if (returnPValue == TRUE){return(st$p.value)}
  return (ifelse(st$p.value > 0.05, TRUE, FALSE))
}



# isnormal <- function(data_vector, returnPValue = FALSE)
# {
#   ifelse(length(data_vector) > 5000, data_vector <- sample(data_vector, 4999), data_vector)
#   st <- shapiro.test(data_vector)
#   if (returnPValue == TRUE){return(st$p.value)}
#   return (ifelse(st$p.value > 0.05, TRUE, FALSE))
# }



#' Find the distribution of a dataframe
#'
#' Takes in a data frame and identifies the data distribution
#' @param data A data frame for which the distribution needs to be identified
#' @param dv A string value specifying the column name of DV
#' @return A data frame with different distribution metrics such as variable type, is_normal, shapiro p value, skewness score, kurtosis score, pearson score, etc.
#' @export
#' @usage data_distribution(data, "dv")

data_distribution<-function(data, dv)
{
  # require(e1071)
  # require(data.table)
  # require(dplyr)
  # require(nortest)
  data<-data.frame(data)
  a<-sapply(data, typeof)
  data_types<-data.frame(a)
  data_types$names<-rownames(data_types)
  for(i in 1:ncol(data))
  {
    a<-data[,i]
    j_t<-data.frame(colnames(data)[i])
    j_t$unique<-length(unique(data[,i]))
    if(i==1)
    {
      j<-j_t
    }
    else{j<-rbind(j,j_t)}
  }

  colnames(j)[1]<-"names"
  colnames(j)[2]<-"distinct_no"
  colnames(data_types)[1]<-"data_type"
  j$names<-as.character(j$names)

  if(nrow(j)==nrow(data_types))
  {
    main_out<-left_join(data_types,j ,by ="names")
  }

  main_out$numeric <- t(data.frame(lapply(data, is.numeric)))
  main_out$distribution <- ifelse((main_out$distinct_no <= 10) | (main_out$data_type == "character") | (main_out$numeric == FALSE), "Categorical", "Continous")
  main_out$is_dv <- ifelse(main_out$names == dv, TRUE, FALSE)

  main_out$is_normal <- NA
  main_out$shapiro_pval <- NA
  main_out$skewness_res <- NA
  main_out$kurtosis_res <- NA
  main_out$pearson_score <- NA

  continuous_cols <- main_out[main_out$distribution == "Continous",]$names

  main_out[main_out$distribution == "Continous",]$is_normal <- t(data.table(data)[,lapply(.SD, isnormal), .SDcols=continuous_cols])

  # for (i in continuous_cols)
  # {
  #   print(i)
  #   isnormal(data.frame(data)[, i])
  # }

  main_out[main_out$distribution == "Continous",]$shapiro_pval <- t(data.table(data)[,lapply(.SD, isnormal, TRUE), .SDcols=continuous_cols])

  main_out[main_out$distribution == "Continous",]$skewness_res <- t(data.table(data)[,lapply(.SD, e1071::skewness), .SDcols=continuous_cols])

  main_out[main_out$distribution == "Continous",]$kurtosis_res <- t(data.table(data)[,lapply(.SD, e1071::kurtosis), .SDcols=continuous_cols])

  main_out[main_out$distribution == "Continous",]$pearson_score <- t(data.table(data)[,lapply(.SD, function(x) {nortest::pearson.test(x)$statistic/nortest::pearson.test(x)$df}), .SDcols=continuous_cols])

  return(main_out[, c("names", "is_dv", "distribution", "is_normal", "shapiro_pval", "kurtosis_res", "pearson_score")])
}

# data_distribution<-function(data, dv)
# {
#   # require(e1071)
#   # require(data.table)
#   # require(dplyr)
#   # require(nortest)
#   data<-data.frame(data)
#   a<-sapply(data, typeof)
#   data_types<-data.frame(a)
#   data_types$names<-rownames(data_types)
#   for(i in 1:ncol(data))
#   {
#     a<-data[,i]
#     j_t<-data.frame(colnames(data)[i])
#     j_t$unique<-length(unique(data[,i]))
#     if(i==1)
#     {
#       j<-j_t
#     }
#     else{j<-rbind(j,j_t)}
#   }
#
#   colnames(j)[1]<-"names"
#   colnames(j)[2]<-"distinct_no"
#   colnames(data_types)[1]<-"data_type"
#   j$names<-as.character(j$names)
#
#   if(nrow(j)==nrow(data_types))
#   {
#     main_out<-left_join(data_types,j ,by ="names")
#   }
#
#   main_out$numeric <- t(data.frame(lapply(data, is.numeric)))
#   main_out$distribution <- ifelse((main_out$distinct_no <= 10) | (main_out$data_type == "character") | (main_out$numeric == FALSE), "Categorical", "Continous")
#   main_out$is_dv <- ifelse(main_out$names == dv, TRUE, FALSE)
#
#   main_out$is_normal <- NA
#   main_out$shapiro_pval <- NA
#   main_out$skewness_res <- NA
#   main_out$kurtosis_res <- NA
#   main_out$pearson_score <- NA
#
#   continuous_cols <- main_out[main_out$distribution == "Continous",]$names
#
#   main_out[main_out$distribution == "Continous",]$is_normal <- t(data.table(data)[,lapply(.SD, isnormal), .SDcols=continuous_cols])
#
#   main_out[main_out$distribution == "Continous",]$shapiro_pval <- t(data.table(data)[,lapply(.SD, isnormal, TRUE), .SDcols=continuous_cols])
#
#   main_out[main_out$distribution == "Continous",]$skewness_res <- t(data.table(data)[,lapply(.SD, e1071::skewness), .SDcols=continuous_cols])
#
#   main_out[main_out$distribution == "Continous",]$kurtosis_res <- t(data.table(data)[,lapply(.SD, e1071::kurtosis), .SDcols=continuous_cols])
#
#   main_out[main_out$distribution == "Continous",]$pearson_score <- t(data.table(data)[,lapply(.SD, function(x) {nortest::pearson.test(x)$statistic/nortest::pearson.test(x)$df}), .SDcols=continuous_cols])
#
#   return(main_out[, c("names", "is_dv", "distribution", "is_normal", "shapiro_pval", "skewness_res", "kurtosis_res", "pearson_score")])
# }

