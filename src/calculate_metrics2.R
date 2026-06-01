# this is a recovered copy from the memory. Had better to retrieve one from HPC
calculate_metrics2 <- function(truth, pred = NULL, threshold = c(0.5, "auto", NA_real_), 
           metrics = c("Prevalence", "AUC", "Sensitivity", 
                       "Specificity", "F1-score", "Diagnostic Odds Ratio"),
           specificity_thre = 0.95)
  {
    # Description
    # This function will calculate the metrics for predictive accuracy
    # given the truth (0/1) and predictions (in probabilities).
    
    # Arguments
    # truth: vector. True values in {0, 1}
    # pred: vector. Predictions in probabilities in [0, 1]
    # threshold: a single numeric value or choose from c(0.5, "auto", NA_real_). 
    # the cut-off value to calculate the prevalence and other accuracy matrics. 
    # AUC will be only calculated when threshold == "auto". 
    # Selection of threshold:
    # "auto", using a threshold which can achieve a specificity of 0.95;
    # NA_real_, the `pred` are all 0s or 1s but not probabilities in between.
    # metrics: accuracy metrics to be calculated.

    res <- rep(NA, length = length(metrics) + 1)
    names(res) <- c(metrics, "Threshold")
    
    if(is.null(pred)){
      if("Prevalence" %in% metrics)
        res["Prevalence"] <- mean(truth)
      return(res)
    }
    
    if(is.na(threshold)){
      temp <- table(predictions = pred, true = truth)
      TN = temp[1,1]
      FP = temp[2,1]
      TP = temp[2,2]
      FN = temp[1,2]
    } else if( (!is.na(threshold)) & is.numeric(threshold)){
      pred_temp = (pred >= threshold) *1
      temp <- table(predictions = pred_temp, true = truth)
      TN = temp[1,1]
      FP = temp[2,1]
      TP = temp[2,2]
      FN = temp[1,2]
    } else if(threshold == "auto"){
      for(thre in 1:99){
        pred_temp <- 1*(pred >= thre/100)
        if(length(unique(pred_temp)) < 2)
          next
        temp <- table(predictions = pred_temp, true = truth)
        TN = temp[1,1]
        FP = temp[2,1]
        TP = temp[2,2]
        FN = temp[1,2]
        Specificity = TN/(TN+FP)
        if(Specificity >= specificity_thre-0.0005) break
      }
      for(thre in seq(thre-1, thre+1, 0.1)){
        pred_temp <- 1*(pred >= thre/100)
        if(length(unique(pred_temp)) < 2)
          next
        temp <- table(predictions = pred_temp, true = truth)
        TN = temp[1,1]
        FP = temp[2,1]
        TP = temp[2,2]
        FN = temp[1,2]
        Specificity = TN/(TN+FP)
        if(Specificity >= specificity_thre-0.0005) break
      }
    }
    
    if("Prevalence" %in% metrics)
      res["Prevalence"] <- (TP+FP)/(TP+TN+FP+FN)
    if("Sensitivity" %in% metrics)
      res["Sensitivity"] <- TP/(TP+FN)
    if("Specificity" %in% metrics)
      res["Specificity"] <- TN/(TN+FP)
    if("F1-score" %in% metrics)
      res["F1-score"] <- 2*TP/(2*TP+FP+FN)
    if("Diagnostic Odds Ratio" %in% metrics)
      res["Diagnostic Odds Ratio"] <- (TP/FN)/(FP/TN)
    if("AUC" %in% metrics)
      if( (!is.na(threshold) && threshold == "auto") |
          (!is.na(threshold) && is.numeric(threshold)  )) {
        perf <- prediction(predictions = pred, labels = truth) %>% 
          performance(., "auc")
        res["AUC"] <- perf@y.values[[1]]
      } else res["AUC"] <- NA 
    if( !is.na(threshold) && !is.numeric(threshold) && threshold == "auto" ){
      res["Threshold"] <- thre/100
    } else res["Threshold"] <- threshold
    
    
    return(res)
  }
