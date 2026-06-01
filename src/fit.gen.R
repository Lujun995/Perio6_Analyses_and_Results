library(tidyverse)
library(ordinalNet)

#remove.lindep <- function()
# There is a remove.lindep function used in `mice` but I feel it does not make
# much sense here. 
# First, we will using penalization which can partially solve multicolinearity.
# Second, we are interested in the predictive power of the models. Multicolinearity
# is of less concern for prediction problem.
# Third, the remove.lindep method used in `mice` is not reasonable enough as it depends
# on a greedy algorithm by sequentially removing the variables which have the 
# highest loading in the eigen vectors of the smallest eigen values. This algorithm
# may not work.
# Consider a group of correlated predictors, LASSO may randomly select a variable 
# from them. It may make sense to prioritize measurements at known/observed teeth.
# For example, eliminate the variables having a colinearity with the measurement
# at the known teeth from the model. 
# However, this "attention" to the observed teeth will destroy model "symmetry",
# which will make the influential network hard to interpret and the model less
# flexible to accept an altenative tooth combination or missing tooth.

fit.gen <- function(data_mat, method = "ordinalNet",
                    definite_model = TRUE, ohx_max = max(data_mat),
                    n.loop = ifelse(definite_model, 1, 20), 
                    n.dataset = ifelse(definite_model, NA_real_, 100), 
                    n.cores = ifelse(definite_model, NA_real_, 50), 
                    random.seed = 123456, verbose = TRUE){
  # Function to fit a generative classifier. 
  # data_mat: input numeric data matrix for training. Can only be complete cases at the moment.
  #Colnames are periodontal measurements, rownames are ids.
  # method: method to construct the predictive model. Other methods are not currently supported. 
  # ohx_max: a universal maximum value used to the fit the model across different variables. 
  #This parameter will ensure every variable will consistently have 0:ohx_max levels 
  # random.seed: random seed used. Default, 123456
  # verbose: If TRUE, the funtion will print out detailed progresses
  # other parameters are currently not supported
  # return value: a list of models (dataset -> loop -> vars -> models),
  # call, initial_values (distribution of periodontal measurements in each var in data_mat),
  # empirical_min and empirical_max (the 0.001 and 0.999 quantile of posterior predictive P(\hat{X_j}|X_{-j})
  # of each var using empirical data in data_mat. They are not the quantiles in X_j.
  # The empirical_min and empirical_max can be useful to detect self-excitation in the Gibbs sampling.
  # When self-excitation, X_{-j} will be very large/small, therefore, P(\hat{X_j}|X_{-j}) can be extreme.)
  library(tidyverse)
  set.seed(seed = random.seed)
  
  stopifnot(method == "ordinalNet")
  stopifnot(definite_model == TRUE)
  stopifnot(inherits(data_mat, "matrix"))
  stopifnot(is.numeric(data_mat))
  
  #ohx_max = max(data_mat)
  initial_values <- sapply(0:ohx_max, data_mat = data_mat, FUN = \(x, data_mat){
    return(colSums(data_mat == x, na.rm = TRUE))
  })
  colnames(initial_values) <- as.character(0:ohx_max)
  
  if(definite_model){
    library(ordinalNet)
    
    vars <- colnames(data_mat)
    vars_list <- vector(mode = "list", length = length(vars))
    names(vars_list) <- vars
    empirical_min <- setNames(rep(-Inf, times = length(vars)), nm = vars)
    empirical_max <- setNames(rep( Inf, times = length(vars)), nm = vars)
    for(y_name in vars){
      if(verbose) cat(y_name, "\n")
      y = factor(data_mat[, y_name], levels = 0:ohx_max, ordered = TRUE)
      fml <- as.formula(
        str_glue("{y_name}~{str_c(base::setdiff(x = vars, y = y_name), collapse = '+')}")
      )
      x <- data_mat[, labels(terms(fml))]
      # x <- remove.lindep(x)
      
      mdl <- ordinalNet(x = x, y = y, family = "cumulative", link = "logit", 
                        alpha=1,
                        # lambdaVals = c(0.005, 0.008, 0.01, 0.015),
                        lambdaVals = c(0.006, 0.008, 0.01, 0.015, 0.02, 0.025),
                        # lambdaVals = 0.01,
                        parallelTerms = TRUE, nonparallelTerms = FALSE,
                        keepTrainingData = FALSE)
      # print(object.size(mdl), units = "auto") # 28 Kb, looks good to save on disk
      vars_list[[y_name]] <- mdl
      
      coefs <- ordinalNet:::coef.ordinalNet(mdl, matrix = FALSE, 
                                            whichLambda = NULL, 
                                            criteria = "aic")
      intercept_index <- grepl(x = names(coefs), pattern = "\\(Intercept\\)")
      # coefs_inter <- coefs[intercept_index]
      coefs_var <- coefs[!intercept_index]
      coefs_var_n0 <- coefs_var[coefs_var != 0]
      if(length(coefs_var_n0) != 0){
        real_beta_X <- x[, names(coefs_var_n0), drop = FALSE] %*% coefs_var_n0 
        empirical_min[y_name] <- quantile(real_beta_X, 0.001) 
        empirical_max[y_name] <- quantile(real_beta_X, 0.999)
      }
    }
    
    res <- list("data=1" = list("step=1" = vars_list), 
                "call" = list(method = method,
                              definite_model = TRUE,
                              n.loop = n.loop, 
                              n.dataset = n.dataset, 
                              n.cores = n.cores, 
                              random.seed = random.seed, 
                              verbose = verbose),
                "initial_values" = initial_values,
                "empirical_min" = empirical_min,
                "empirical_max" = empirical_max)
    
  } #else {}
  
  return(res)
}