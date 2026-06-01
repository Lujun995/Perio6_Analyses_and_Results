library(tidyverse)
library(ordinalNet)

predict.ordinalNet <- 
  function (mdl, newx = NULL, whichLambda = NULL, 
            empirical_max = NULL, empirical_min = NULL,
            criteria = c("aic", "bic"), 
            type = c("response", "class", "link", "beta_X"), ...) {
    # Function to predict based on an ordinalNet object. Modified from ordinalNet:::predict.ordinalNet
    # Inputs: 
    # mdl: model fit from ordinalNet.
    # newx: new data. a numeric matrix
    # whichLambda" which lambda value to be used in the ordinalNet fit
    # empirical_max, empirical_min: these values restrict the output in a range
    criteria <- match.arg(criteria)
    type <- match.arg(type)
    if (is.null(newx) && !mdl$args$keepTrainingData) 
      stop(paste0("Model was fit with keepTrainingData=FALSE, so training data was not saved. ", 
                  "A newx argument is required."))
    nLev <- mdl$nLev
    xNames <- mdl$xNames
    parallelTerms <- mdl$args$parallelTerms
    nonparallelTerms <- mdl$args$nonparallelTerms
    family <- mdl$args$family
    link <- mdl$args$link
    reverse <- mdl$args$reverse
    linkfun <- if (is.null(mdl$args$customLink)) 
      ordinalNet:::makeLinkfun(family, link) else 
        mdl$args$customLink
    if (is.null(newx)) 
      newx <- mdl$args$x
    
    coefs <- ordinalNet:::coef.ordinalNet(mdl, matrix = FALSE, 
                                          whichLambda = whichLambda, 
                                          criteria = criteria)
    intercept_index <- grepl(x = names(coefs), pattern = "\\(Intercept\\)")
    coefs_inter <- coefs[intercept_index]
    coefs_var <- coefs[!intercept_index]
    coefs_var_n0 <- coefs_var[coefs_var != 0]
    # it will automatically cause error of "subscript out of bound" 
    # if not all coefs_var_n0 are in newx
    beta_X <- newx[,names(coefs_var_n0), drop = FALSE] %*% coefs_var_n0 
    if(!is.null(empirical_min))
      beta_X <- pmax(beta_X, empirical_min)
    if(!is.null(empirical_max))
      beta_X <- pmin(beta_X, empirical_max)
    if(type == "beta_X")
      return(beta_X)
    etaMat <- sapply(coefs_inter, FUN = `+`, beta_X, simplify = TRUE)
    
    # real_beta_X <- as.matrix(nhanes_complete_num[,names(coefs_var_n0)]) %*% coefs_var_n0
    # empirical_min = quantile(real_beta_X, 0.001) # for real-world data, the 0.1% percentile is -20.45009
    # # min(real_beta_X) # -25.0045, still way larger than beta_X
    # empirical_max = quantile(real_beta_X, 0.999) # for real-world data, the 99.9% percentile is 0.1517
    # # the intercept for 6 is 20.2823, while the intercept for 1 is 2.417980
    
    deltaNames <- ordinalNet:::getDeltaNames(family, reverse, nLev) 
    # colnames(betaMat) <- paste0(link, "(", deltaNames, ")") # from ordinalNet:::coef.ordinalNet
    # colnames(etaMat) <- colnames(betaMat)
    colnames(etaMat) <- paste0(link, "(", deltaNames, ")")
    
    ## Original algorithm:
    # betaMat <- ordinalNet:::coef.ordinalNet(mdl, matrix = TRUE,
    #                                         whichLambda = whichLambda,
    #                                         criteria = criteria)
    # etaMat2 <- cbind(1, newx) %*% betaMat
    # all.equal(etaMat, etaMat2) # TRUE
    if (type == "link") 
      return(etaMat) # it will break the function here
    probMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, 
    ])))
    probMat <- cbind(probMat, 1 - rowSums(probMat))
    if (reverse) 
      probMat <- probMat[, nLev:1]
    colnames(probMat) <- paste0("P[Y=", 1:nLev, "]")
    if (type == "response") 
      return(probMat)
    class <- c(apply(probMat, 1, which.max))
    class
  }

check.validity <- function(object, gen_output, training_data = NULL, 
                           upper_quantile = 0.999, lower_quantile = 0.001){
  # calculate how many generated observations may exceed the empirical 
  # maxima and mimima. This function can detect improper generated observations
  # due to the self-excitation during Gibbs sampling
  exceed_upper_quantile <- exceed_lower_quantile <-
    matrix(FALSE, nrow = nrow(gen_output), ncol = ncol(gen_output),
           dimnames = list(NULL, colnames(gen_output)))
  if(object$call$definite_model){
    vars <- names(object[["data=1"]][["step=1"]])
    for (y_name in vars) {
      mdl <- object[["data=1"]][["step=1"]][[y_name]]
      coefs <- ordinalNet:::coef.ordinalNet(mdl, matrix = FALSE, 
                                            whichLambda = NULL, 
                                            criteria = "aic")
      intercept_index <- grepl(x = names(coefs), pattern = "\\(Intercept\\)")
      # coefs_inter <- coefs[intercept_index]
      coefs_var <- coefs[!intercept_index]
      coefs_var_n0 <- coefs_var[coefs_var != 0]
      if(is.null(training_data)){
        empirical_min <- object$empirical_min[y_name]
        empirical_max <- object$empirical_max[y_name]
      } else if(length(coefs_var_n0) == 0){
        empirical_min <- -Inf
        empirical_max <-  Inf
      } else {
        real_beta_X <- training_data[, names(coefs_var_n0)] %*% coefs_var_n0 
        empirical_min <- quantile(real_beta_X, lower_quantile) 
        empirical_max <- quantile(real_beta_X, upper_quantile) 
      }
      
      beta_X <- gen_output[, names(coefs_var_n0)] %*% coefs_var_n0
      exceed_lower_quantile[, y_name] <- beta_X < empirical_min
      exceed_upper_quantile[, y_name] <- beta_X > empirical_max
    }
  }
  
  return(list(exceed_upper_quantile = exceed_upper_quantile,
              exceed_lower_quantile = exceed_lower_quantile))
}

gibbs.sampler <- function(object, newdata, input_teeth, mis = NULL, ohx_max = 6,
                          dataset.no = 1:100, n.loop = 20, verbose = TRUE){
  # Function to predict the unobserved measurements using Gibbs sampling.
  # these measurements can be further used to predict the periodontitis status.
  #
  # newdata: df, with names like "ohx02pcd" ... assuming the input values 
  # are ordinal but in a numeric form
  # input_teeth: measurements observed. character (variable names)
  # object: a list from fit.gen, including a list of models (dataset -> loop -> vars -> models),
  # call, initial_values, and empirical_min and empirical_max
  #  (the 0.001 and 0.999 quantile of posterior predictive P(\hat{X_j}|X_{-j})
  # of each var using empirical data in data_mat. They are not the quantiles in X_j.
  # The empirical_min and empirical_max can be useful to detect self-excitation in the Gibbs sampling.
  # When self-excitation, X_{-j} will be very large/small, therefore, P(\hat{X_j}|X_{-j}) can be extreme.)
  # dataset.no: datasets No.s in the object to be used for Gibbs sampling. When object is a definite_model, dataset.no does not matter
  # n.loop: loops in the Gibbs sampling.
  # mis: missing data matrix (at least have all at input_teeth)
  # 
  # output: df_full, the complete dataset with unobserved values filled in
  
  n.dataset = max(dataset.no)
  if(!object$call$definite_model){
    # if the model is definite,such as ordinalNet in complete cases
    # we can generate datasets for abitrary times
    # otherwise, the predict function can only predict less than object$n.loop and object$n.dataset times
    if(n.loop > object$call$n.loop)
      stop("n.loop > object$call$n.loop")
    if(n.dataset > object$call$n.dataset)
      stop("n.dataset > object$call$n.dataset")
  } 
  df_full <- vector(mode = "list", length = length(dataset.no))
  names(df_full) <- paste0("imp", dataset.no)
  for(dtst in dataset.no){
    if(verbose) cat("Dataset No.", dtst, "\n")
    if(object$call$definite_model) {
      current_dtst <- object[["data=1"]]
    } else {
      current_dtst <- object[[str_glue("data={dtst}")]]
    }
    # initialize
    # it works for both missingness in input_teeth 
    # and unobserved teeth (all values are missing due to masking in predict.gen)
    # this will be controlled by the random seed in the environment
    newdata_imp <- newdata
    for(varr in colnames(newdata_imp)){
      var_x <- newdata_imp[, varr]
      na.ind <- is.na(var_x)
      if(sum(na.ind) > length(var_x)*0.7){# if most of the data (>70%) are NAs.
        # Sample from initial values in NHANES
        var_x[na.ind] <- sample(x = 0:ohx_max, size = sum(na.ind), replace = TRUE,
                                prob = object$initial_values[varr,])
      } else {
        # sample from the data
        var_x[na.ind] <- sample(x = var_x[!na.ind], size = sum(na.ind),
                                replace = TRUE)
      }
      newdata_imp[, varr] <- var_x
    }
    ####################################################
    # set.seed(seed= 123456)
    # newdata_imp <- newdata[rep(5, 10),]
    # for(varr in names(newdata_imp)){
    #   var_x <- newdata_imp[[varr]]
    #   na.ind <- is.na(var_x)
    #   if(sum(na.ind) == length(var_x)){
    #     # all NAs. Sample from initial values in NHANES
    #     var_x <- sample(x = 0:6, size = sum(na.ind), replace = TRUE)
    #                     # prob = object$initial_values[varr,])
    #   } else {
    #     var_x[na.ind] <- sample(x = var_x[!na.ind], size = sum(na.ind),
    #                             replace = TRUE)
    #   }
    #   newdata_imp[[varr]] <- var_x
    # }
    ##########################################################
    for(loop in 1:n.loop){
      if(verbose) cat("Loop No.", loop, "\t")
      #########################################
      # if(verbose) cat("Loop No.", loop, "\n")
      #########################################
      if(object$call$definite_model){
        current_loop <- current_dtst[["step=1"]]
      } else {
        current_loop <- current_dtst[[str_glue("step={loop}")]]
      }
      vars <- names(current_loop)
      # predict repeatedly using the estimated parameters (models)
      for(y_name in vars){
        if((y_name %in% input_teeth) & is.null(mis))
          next # no need to impute on input_teeth
        newx <- newdata_imp[, current_loop[[y_name]]$xNames]
        nodes <- predict.ordinalNet(mdl = current_loop[[y_name]], 
                                    newx = newx, type = "response") 
        # #################################################
        # if(y_name == "ohx06las"){
        #   print(newdata_imp[[y_name]])
        #   print(nodes)
        # }
        # #################################################
        nodes <- apply(nodes, MARGIN = 2, FUN = \(x) {
          x[x < 1e-10] <- 0
          return(x)
        })
        # nodes: the prob. of different levels of the ordinal variable y
        imputes <- apply(nodes, MARGIN = 1, FUN = function(s){
          # the mdl will have intercepts for unused levels
          # so the nodes always have prob of 0:ohx_max
          sample(0:ohx_max, size = 1, prob = s)
        }) # imputes: are integers ranging 0-ohx_max
        # ####################################
        # if(y_name == "ohx06las"){
        #   print(imputes)
        # }
        # ####################################
        # update the dataset
        if(!(y_name %in% input_teeth)) # unobserved teeth, update all
          newdata_imp[, y_name] <- imputes else # input teeth, update unobserved values
            newdata_imp[mis[, y_name], y_name] <- imputes[mis[, y_name]] # update the dataset
      }
    }
    df_full[[paste0("imp", dtst)]] <- newdata_imp
  }
  return(df_full)
}