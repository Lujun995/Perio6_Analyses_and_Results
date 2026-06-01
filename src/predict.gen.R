library(tidyverse)
library(ordinalNet)
library(parallel)

pred.gen <- function(cal_mdl, ppd_mdl, new_data_cal, new_data_ppd, input_teeth,
                     check_validity = TRUE, valid_threshold = 0.05,
                     dataset_no = 1:50, n_loop = 20, verbose = FALSE,
                     type = c("data", "results"), mis = NULL){
  # wrapper function running in a single thread
  # using cal_mdl, ppd_mdl from fit.gen
  # new_data_cal, new_data_ppd: a data matrix with colnames like "ohx02lad" for 
  # new_data_cal and "ohx02pcd" for new_data_ppd and row names of ids. 
  # the input values are numbers/integers ranging 0-6. 
  # input_teeth: observed measurement, whose values will not be imputed
  # dataset_no: which dataset to be used in cal_mdl and ppd_mdl
  # type: return values. data, imputed data; results, MS and severe probabilities
  
  n_gen = length(dataset_no)
  
  # generating
  new_data_cal_gen <- 
    gibbs.sampler(object = cal_mdl, newdata = new_data_cal,
                  input_teeth = input_teeth, mis = mis,
                  dataset.no = dataset_no, n.loop = n_loop, verbose = verbose)
  new_data_cal_gen <- do.call(rbind, new_data_cal_gen)
  new_data_ppd_gen <- 
    gibbs.sampler(object = ppd_mdl, newdata = new_data_ppd,
                  input_teeth = input_teeth, mis = mis,
                  dataset.no = dataset_no, n.loop = n_loop, verbose = verbose)
  new_data_ppd_gen <- do.call(rbind, new_data_ppd_gen)
  
  if(TRUE){
    # check for self-excitation
    new_data_cal_gen_valid <- check.validity(object = cal_mdl,
                                             gen_output = new_data_cal_gen)
    new_data_ppd_gen_valid <- check.validity(object = ppd_mdl,
                                             gen_output = new_data_ppd_gen)
    
    # discard invalid samples
    valid_index_cal <- 
      rowMeans(new_data_cal_gen_valid$exceed_lower_quantile + 
                 new_data_cal_gen_valid$exceed_upper_quantile) < valid_threshold
    valid_index_ppd <- 
      rowMeans(new_data_ppd_gen_valid$exceed_lower_quantile + 
                 new_data_ppd_gen_valid$exceed_upper_quantile) < valid_threshold
    
    # combine the dataset
    stopifnot(identical(rownames(new_data_ppd_gen),rownames(new_data_cal_gen)))
    new_data_gen <- cbind(new_data_cal_gen, 
                          new_data_ppd_gen) %>%
      as_tibble(rownames = NA) %>%
      rownames_to_column(var = "seqn") %>%
      mutate(valid_index = valid_index_cal & valid_index_ppd)
    if(type == "data"){
      return(new_data_gen)
    }
    res <- new_data_gen %>%
      mutate(temp_id = row_number()) %>%
      case_def(id = "temp_id") %>%
      mutate(MS = (mod_12|sev_12)) %>%
      group_by(seqn) %>%
      summarise(n_val = sum(valid_index), pred_MS = mean(MS[valid_index]), 
                pred_sev_12 = mean(sev_12[valid_index])) %>%
      mutate(pred_MS = if_else(n_val <= n_gen/5, 1, pred_MS),
             pred_sev_12 = if_else(n_val <= n_gen/10, 1, pred_sev_12))
  }
  return(res)
}

predict.gen <- 
  function(cal_mdl, ppd_mdl, new_data, input_teeth, cl = NULL,
           n_cores = 5, n_rep = 50, split_size = min(100, nrow(new_data)),
           check_validity = TRUE, valid_threshold = 0.05, 
           type = c("data", "results"), ...){
  # Function to predict periodontitis status based on object and newdata
  #
  # cal_mdl, ppd_mdl: object created from fit.gen
  # new_data: a data matrix!!! with colnames like "ohx02pcd" and row names of ids. 
  # the input values are numbers/integers (ordinal but in a numeric form)
  # input_teeth: teeth/measurements in the newdata to be used for prediction
  # n_cores: number of threads to use
  # split_size: dataset will be split into smaller blocks of split_size and submitted to different cores
  # check_validity, valid_threshold: to check the validity of the Gibbs sampling 
  # output, avoiding self excitation in the sampling
  # n_rep: total number of output datasets in Gibbs sampling generated. 
  # This determines the digits of predicted probability (like, how many digits in P-values do you want
  # in permutation test).
  # ...: other parameters supported in gibbs.sampler
  #
  # output: an dataframe similar with newdata but have predicted periodontitis status 
    
    type <- match.arg(type)
    if(n_cores == 1 | nrow(new_data) == 1){
      res <- predict.gen.single(cal_mdl, ppd_mdl, new_data, input_teeth,
                                n_rep, split_size = min(100, n_rep * nrow(new_data)), # override the split_size
                                check_validity, valid_threshold, type)
      return(res)
    }
      
    library(parallel)
    # inside the function
    cl_flag = !is.null(cl)
    if(!cl_flag){
      cl <- makeCluster(n_cores)
    }
    
    
    all_ohx <- colnames(new_data)
    test_data <- new_data
    test_data[, !(colnames(test_data) %in% input_teeth)] <- NA # masking
    mis_all <- is.na(test_data[, input_teeth])
    if(sum(mis_all) == 0)
      mis_all <- NULL # complete observation and no missing
    # replicate the test_data by n.rep.per.core times. The rownames can be used to 
    # identify the results from the same individual
    # test_data <- do.call(rbind, replicate(n.rep.per.core, test_data, simplify=FALSE))
    # predict CAL and PPD seperately:
    test_data_cal <- test_data[, str_subset(all_ohx, pattern = "ohx\\d\\dla[dspaml]")]
    test_data_ppd <- test_data[, str_subset(all_ohx, pattern = "ohx\\d\\dpc[dspaml]")]
    
    nr <- nrow(new_data)
    n_split <- ceiling(nr/split_size)
    
    clusterExport(cl, envir = environment(),
                  c("test_data_cal", "test_data_ppd", "cal_mdl", "ppd_mdl",
                    "input_teeth", "mis_all", "pred.gen", "case_def", 
                    "predict.ordinalNet", "check.validity", "gibbs.sampler"))
    clusterEvalQ(cl, {
      library(tidyverse)
      NULL
    })
    temp_func <- \(x){ # x is a control number for dataset_no and split_no. x = 1:(n_rep*n_split)
      dataset_no = (x-1)%/%n_split + 1 # 1:n_rep
      split_no = (x-1)%%n_split + 1 # 1:n_split
      # View(cbind(dataset_no, split_no))
      current_split <- ((split_no-1)*split_size+1) : min(split_no*split_size, nr)
      mis <- if(is.null(mis_all)) NULL else mis_all[current_split,]
      test_data_gen <- 
        pred.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, 
                 new_data_cal = test_data_cal[current_split,], 
                 new_data_ppd = test_data_ppd[current_split,],
                 input_teeth = input_teeth, mis = mis,
                 check_validity = TRUE, valid_threshold = 0.05,
                 dataset_no = dataset_no, n_loop = 20, verbose = FALSE, type = "data")
      # if dataset_no > the dataset_no in the model, gibbs.sampler will return error
      return(test_data_gen)
    }
    
    res <- parSapplyLB(cl = cl, X = 1:(n_rep*n_split), FUN = temp_func, simplify = FALSE)
    res <- bind_rows(res)
    if(cl_flag){
      clusterEvalQ(cl, {
        rm(list = ls())
        NULL
      })
    } else {
      stopCluster(cl)
    }
    
    if(type == "data")
      return(res)
    res <- res %>%
      mutate(temp_id = row_number()) %>%
      case_def(id = "temp_id") %>%
      mutate(MS = (mod_12|sev_12)) %>%
      group_by(seqn) %>%
      summarise(n_val = sum(valid_index), pred_MS = mean(MS[valid_index]), 
                pred_sev_12 = mean(sev_12[valid_index])) %>%
      mutate(pred_MS = if_else(n_val <= n_rep/5, 1, pred_MS),
             pred_sev_12 = if_else(n_val <= n_rep/10, 1, pred_sev_12))
    return(res)
}

predict.gen.single <- 
  function(cal_mdl, ppd_mdl, new_data, input_teeth,
           n_rep = 50, split_size = min(100, n_rep * nrow(new_data)),
           check_validity = TRUE, valid_threshold = 0.05, 
           type = c("data", "results"), ...){
    # this function need to be updated in that it assumes the model is always a deterministic model
    type <- match.arg(type)
    library(tidyverse)
    
    if(!cal_mdl$call$definite_model | !ppd_mdl$call$definite_model)
      stop("predict.gen.single only allows deterministic model for now!")
    
    all_ohx <- colnames(new_data)
    test_data <- new_data
    test_data[, !(colnames(test_data) %in% input_teeth)] <- NA # masking
    # replicate the test_data. The rownames can be used to identify the results from the same individual
    #this assumes the model was deterministic
    test_data <- do.call(rbind, replicate(n_rep, test_data, simplify=FALSE)) 
    mis_all <- is.na(test_data[, input_teeth])
    if(sum(mis_all) == 0)
      mis_all <- NULL
    # predict CAL and PPD seperately:
    test_data_cal <- test_data[, str_subset(all_ohx, pattern = "ohx\\d\\dla[dspaml]"), drop = FALSE]
    test_data_ppd <- test_data[, str_subset(all_ohx, pattern = "ohx\\d\\dpc[dspaml]"), drop = FALSE]
    
    # nr <- nrow(new_data)
    nr <- nrow(new_data) * n_rep
    n_split <- ceiling(nr/split_size)
    
    temp_func <- \(x){ # x is a control number for dataset_no and split_no. x = 1:(n_rep*n_split)
      # dataset_no = (x-1)%/%n_split + 1 # 1:n_rep
      # split_no = (x-1)%%n_split + 1 # 1:n_split
      # View(cbind(dataset_no, split_no))
      dataset_no = 1
      split_no = x
      current_split <- ((split_no-1)*split_size+1) : min(split_no*split_size, nr)
      mis <- if(is.null(mis_all)) NULL else mis_all[current_split,]
      test_data_gen <- 
        pred.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, 
                 new_data_cal = test_data_cal[current_split,], 
                 new_data_ppd = test_data_ppd[current_split,],
                 input_teeth = input_teeth, mis = mis,
                 check_validity = TRUE, valid_threshold = 0.05,
                 dataset_no = dataset_no, n_loop = 20, verbose = FALSE, type = "data")
      # if dataset_no > the dataset_no in the model, gibbs.sampler will return error
      return(test_data_gen)
    }
    
    res <- sapply(X = 1:n_split, FUN = temp_func, simplify = FALSE)
    res <- bind_rows(res)

    if(type == "data")
      return(res)
    res <- res %>%
      mutate(temp_id = row_number()) %>%
      case_def(id = "temp_id") %>%
      mutate(MS = (mod_12|sev_12)) %>%
      group_by(seqn) %>%
      summarise(n_val = sum(valid_index), pred_MS = mean(MS[valid_index]), 
                pred_sev_12 = mean(sev_12[valid_index])) %>%
      mutate(pred_MS = if_else(n_val <= n_rep/5, 1, pred_MS),
             pred_sev_12 = if_else(n_val <= n_rep/10, 1, pred_sev_12))
    return(res)
  }


