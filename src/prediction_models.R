library(ranger)
library(ROCR)
library(glmnet)
library(e1071)
library(tidyverse)

prediction_models <- function(training = train.d_orgin, test = NULL,
                              using = partial_mouth, 
                              method = c("naive", "svm", "rf", "penalized"),
                              metrics = c("Prevalence", "AUC", "Sensitivity", 
                                          "Specificity", "F1-score", 
                                          "Diagnostic Odds Ratio"),
                              outcome = "MS", verbose = TRUE){
  # Description
  # This function will calculate the predictive power of some competing 
  # partial-mouth diagnostic protocols, including naive method, svm, 
  # random forest (rf), penalized regression.
  
  # Arguments
  # training: training dataset, obs with NA measurement in `using` will be excluded
  # test: test dataset, obs with NA measurement in `using` will be excluded
  # using: vector of string, features used for the model
  # method: string, or string vector. one of "naive", "svm", "rf", "penalized"
  # metrics: string or string vector. one of "Prevalence", "AUC", "Sensitivity", 
  #"Specificity", "F1-score", "Diagnostic Odds Ratio"
  # outcome: string, the name of the outcome variable. the outcome variable can 
  #be integer or factor or logic
  # verbose: logic, TRUE/FALSE, to print the details
  
  # Details
  # The parameters used in the ML-based protocols:
  # SVM: the default parameters, which works well in the NHANES dataset
  # random forest: the default parameters, which works well.
  # penalized regression: alpha = 0.5, 0.5 L1 penalty + 0.5 L2 penalty, which
  # works well. The lambda values was determined through cross validation.
  
  # Value
  # A list of 
  # metrics_mat: predictive power matrix
  # pred_mat: predicted values and truth
  # n_train, n_test: number of obs in the training and test datasets
  
  training = drop_na(training[, c(using, outcome)])
  n_train = nrow(training)
  n_cases = table(training[[outcome]])
  n_CV = 10
  if(verbose){
    cat("the training data set have", n_train, "complete cases.\n")
    cat("the outcome variable:", print(n_cases))
  }
  n_test = NULL
  pred = NULL
  if(!is.null(test)){
    test = drop_na(test[, c(using, outcome)])
    n_test = nrow(test)
    if(verbose){
      cat("the test data set have", n_test, "complete cases.\n")
    }
  }
  
  training[[outcome]] = as.numeric(training[[outcome]])
  if(!is.null(test)) 
    test[[outcome]] = as.numeric(test[[outcome]])
  metrics_mat <- matrix(nrow = length(method),
                        ncol = length(metrics),
                        dimnames = list(method , metrics))
  if(!is.null(test)){
    split_no <- NULL
    pred_mat <- matrix(nrow = n_test,
                       ncol = length(method)+1,
                       dimnames = list(NULL, c(method, "truth")))
    pred_mat[,"truth"] <- test[[outcome]]
  } else {
    split_no <- sample(cut(1:n_train, breaks = n_CV, labels = FALSE)) # assign 1:n_CV to each row
    pred_mat <- matrix(nrow = n_train,
                       ncol = length(method)+1,
                       dimnames = list(NULL, c(method, "truth")))
    pred_mat[,"truth"] <- training[[outcome]]
  }
  
  if("naive" %in% method){
    if(!is.null(test)){
      n_teeth <- length(unique(str_extract(using,"\\d\\d")))
      temp_d <- test[, using]
      temp_d$"seqn" <- 1:nrow(temp_d)
      d <- case_def(temp_d, n_teeth_total = n_teeth) %>%
        mutate(MS = (mod_12|sev_12))
      pred <- as.numeric(d[[outcome]])
      pred_mat[, "naive"] <- pred
      metrics_mat["naive", metrics] <- 
        calculate_metrics2(truth = test[[outcome]], 
                           pred = pred, threshold = NA_real_)[metrics]
    } else {
      n_teeth <- length(unique(str_extract(using,"\\d\\d")))
      temp_d <- training[, using]
      temp_d$"seqn" <- 1:nrow(temp_d)
      d <- case_def(temp_d, n_teeth_total = n_teeth) %>%
        mutate(MS = (mod_12|sev_12))
      pred <- as.numeric(d[[outcome]])
      pred_mat[, "naive"] <- pred
      metrics_mat["naive", metrics] <- 
        calculate_metrics2(truth = training[[outcome]], 
                           pred = pred, threshold = NA_real_)[metrics]
    }
  }
  if("svm" %in% method){
    if(!is.null(test)){
      mdl <- svm(x = training[, using], y = training[[outcome]])
      pred <- predict(object = mdl, newdata = test[, using], type="response")
      pred_mat[, "svm"] <- pred
      metrics_mat["svm", metrics] <- 
        calculate_metrics2(truth = test[[outcome]], 
                           pred = pred, threshold = "auto")[metrics]
    } else {
      for(ii in 1:n_CV){
        mdl <- svm(x = training[split_no != ii, using], 
                   y = training[[outcome]][split_no != ii])
        pred <- predict(object = mdl, type="response",
                        newdata = training[split_no == ii, using])
        pred_mat[split_no == ii, "svm"] <- pred
      }
      pred <- pred_mat[, "svm"]
      metrics_mat["svm", metrics] <- 
        calculate_metrics2(truth = training[[outcome]], 
                           pred = pred, threshold = "auto")[metrics]
    }
  } 
  if("rf" %in% method){
    if(!is.null(test)){
      mdl <- ranger(x = training[, using], y = training[[outcome]],
                    classification = TRUE, probability = TRUE, 
                    always.split.variables = NULL)
      pred <- predict(object = mdl, data = test, type="response")
      # for random forest, when using 0-1 numeric values as binary outcome,
      # the first column in prediction indicates the outcome type
      # in the first row in the training, so it can be 0 or 1.
      pred <- pred$predictions[,ifelse(training[[outcome]][1] == 1, 1, 2)]
      pred_mat[, "rf"] <- pred
      metrics_mat["rf", metrics] <- 
        calculate_metrics2(truth = test[[outcome]], 
                           pred = pred, threshold = "auto")[metrics]
    } else {
      for(ii in 1:n_CV){
        mdl <- ranger(x = training[split_no != ii, using], 
                      y = training[[outcome]][split_no != ii],
                      classification = TRUE, probability = TRUE, 
                      always.split.variables = NULL)
        pred <- predict(object = mdl, type="response",
                        data = training[split_no == ii, using])
        pred <- pred$predictions[,ifelse(training[[outcome]][split_no != ii][1] == 1, 1, 2)]
        pred_mat[split_no == ii, "rf"] <- pred
      }
      pred <- pred_mat[, "rf"]
      metrics_mat["rf", metrics] <- 
        calculate_metrics2(truth = training[[outcome]], 
                           pred = pred, threshold = "auto")[metrics]
    }
  }
  if("penalized" %in% method){
    if(!is.null(test)){
      mdl <- glmnet(x = as.matrix(training[, using]), y = training[[outcome]],
                    family = "binomial", alpha = 0.5)
      mdl_cv <- cv.glmnet(x = as.matrix(training[, using]), 
                          y = training[[outcome]],
                          type.measure="class", family = "binomial", alpha = 0.5)
      pred <- predict(mdl, newx = as.matrix(test[, using]), 
                      s = mdl_cv$lambda.min, type = "response")
      pred_mat[, "penalized"] <- pred
      metrics_mat["penalized", metrics] <- 
        calculate_metrics2(truth = test[[outcome]], 
                           pred = pred, threshold = "auto")[metrics]
    } else {
      for(ii in 1:n_CV){
        mdl <- glmnet(x = as.matrix(training[split_no != ii, using]), 
                      y = training[[outcome]][split_no != ii],
                      family = "binomial", alpha = 0.5)
        mdl_cv <- cv.glmnet(x = as.matrix(training[split_no != ii, using]), 
                            y = training[[outcome]][split_no != ii],
                            type.measure="class", family = "binomial", alpha = 0.5)
        pred <- predict(mdl, newx = as.matrix(training[split_no == ii, using]), 
                        s = mdl_cv$lambda.min, type = "response")
        pred_mat[split_no == ii, "penalized"] <- pred
      }
      pred <- pred_mat[, "penalized"]
      metrics_mat["penalized", metrics] <- 
        calculate_metrics2(truth = training[[outcome]], 
                           pred = pred, threshold = "auto")[metrics]
    }
  }
  return(list(metrics_mat = metrics_mat, pred_mat = pred_mat,
              n_train=n_train, n_test=n_test, split_no = split_no))
}