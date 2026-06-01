##############################
#Find a best n-tooth combination for enhanced partial-mouth diagnosis based on
# a machine-learning classifier using evolutionary algorithm.
##############################
find_row_index <- function(record, mat) { # check if record has appeared in the mat
  idx <- which(apply(mat, 1, function(row) all(row == record)))
  if (length(idx) > 0) return(idx[1]) else return(0)
}

workdir = "C:/Users/china/Desktop/UMN Dental health"
Date = "0805" # for version control
# n.cores = 100
# n.cores = 2
# nboot = 1 # number of random splits to evaluate
#train_prop = 0.7 # 70-30 split. Not used. The model parameters in a generative classifier have been pre-specified. So, no random split required.
n_teeth = 6 # a 6-tooth-combination

library(tidyverse)
# library(parallel)
# library(ordinalNet)
# options(Ncpus = n.cores)
# library(ranger)
# library(mlr)
set.seed(seed= 1233455)
# setwd(workdir)
# source(str_glue("{workdir}/src/case_def.R"))
# source(str_glue("{workdir}/src/predict.gen.dependencies.R")) 
# source(str_glue("{workdir}/src/predict.gen.R"))
source(str_glue("{workdir}/src/calculate_metrics2.R"))

sum.na <- function(x) sum(x,na.rm=T)
NA_num <- function(x) sum(is.na(x))
case_definition <- function(d){
  # this is based on 2012 ADA/CDC case definition
  # A very dirty way...
  
  #only the pc and la columns
  num.teeth <- names(d) %>% str_subset("^ohx\\d\\d") %>% 
    str_extract("\\d\\d") %>% unique() %>% length
  pc <- d %>% select(contains("pc"))
  if(ncol(pc) != num.teeth*4) 
    stop("Number of pc measurement != 4*number of teeth")
  la <- d %>% select(contains("la"))
  if(ncol(la) != num.teeth*4) 
    stop("Number of la measurement != 4*number of teeth")
  
  # No
  ps <- rep("no", length = nrow(d)) 
  
  # Mild
  al.g3 <- la >= 3
  out <- 0
  for (i in 1:num.teeth){
    temp <- apply(al.g3[,(4*i-3):(4*i)], 1, sum.na) > 0
    out <- temp + out
  }
  mild.al <- out >= 2
  
  pc.g4 <- pc >= 4
  out <- 0
  for (i in 1:num.teeth){
    temp <- apply(pc.g4[,(4*i-3):(4*i)], 1, sum.na) > 0
    out <- temp + out
  }
  mild.pc4 <- out >= 2
  
  pc.g5 <- pc >= 5
  mild.pc5 <- apply(pc.g5, 1, sum.na) > 0
  
  # d$ps[((mild.pc4 * mild.al)+mild.pc5)>0] <- "mild"
  ps[(((mild.pc4 +mild.pc5)>0)*mild.al)>0] <- "mild"
  
  # Moderate
  al.g4 <- la >= 4
  out <- 0
  for (i in 1:num.teeth){
    temp <- apply(al.g4[,(4*i-3):(4*i)], 1, sum.na) > 0
    out <- temp + out
  }
  moderate.al <- out >= 2
  
  pc.g5 <- pc >= 5
  out <- 0
  for (i in 1:num.teeth){
    temp <- apply(pc.g5[,(4*i-3):(4*i)], 1, sum.na) > 0
    out <- temp + out
  }
  moderate.pc <- out >= 2
  
  ps[(moderate.pc+moderate.al) > 0] <- "mod"
  
  # Severer
  al.g6 <- la >= 6
  out <- 0
  for (i in 1:num.teeth){
    temp <- apply(al.g6[,(4*i-3):(4*i)], 1, sum.na) > 0
    out <- temp + out
  }
  severer.al <- out >= 2
  
  pc.g5 <- pc >= 5
  severer.pc5 <- apply(pc.g5, 1, sum.na) > 0
  
  ps[(severer.al*severer.pc5)>0] <- "sev"
  S <- ps=="sev"
  MS <- ps=="sev" | ps=="mod"
  d <- d %>% mutate(ps = ps, S = S, MS = MS)
  return(d)
}

prefix = "ohx"
suffix = c(str_c("pc", c("d", "m", "s", "p", "l", "a"), sep = ""), 
           str_c("la", c("d", "m", "s", "p", "l", "a"), sep = "")) # pcd, ..., lad, ...
all_tooth <- sprintf("%02d", c(seq(2,15), seq(18,31))) # all teeth in NHANES
all_ohx <- str_c(prefix, 
                 rep(all_tooth, each = length(suffix)), 
                 rep(suffix, times = length(all_tooth)))
training_data_raw <- readRDS(str_glue("{workdir}/res/nhanes_complete.rds"))
training_data <- training_data_raw %>%
  mutate(across(.cols = starts_with("ohx"),
                .fns = \(x) case_match(x, "0"~0, "1"~1, "2"~2, "3"~3,
                                       "4"~4, "5"~5, "6/6+"~6))) # %>%
# column_to_rownames(var = "seqn") %>%
# select(starts_with("ohx")) %>%
# as.matrix()
#rand_index <- sapply(X = 1:nboot, 
#                     FUN = \(x) return(sample.int(
#                       n = nrow(training_data),
#                       size = round(train_prop * nrow(training_data)),
#                       replace = FALSE, prob = NULL, useHash = FALSE)))

mtrcs = c("Prevalence", "AUC", "Sensitivity","Specificity", "F1-score",
          "Diagnostic Odds Ratio")
# cal_mdl <- readRDS(str_glue("{workdir}/res/nhanes_cal_genmdl_ord.rds"))
# ppd_mdl <- readRDS(str_glue("{workdir}/res/nhanes_ppd_genmdl_ord.rds"))
# well, these model were trained on the whole NHANES dataset, so there can be some data leakage issue here.
# But that should be no problem because we wants to find the best tooth combinations in the whole NHANES dataset
# using these models. The random split was just to introduce some randomness to avoid EA falling into some 
# regional maxima.

# specify the ML method to use
# get_metrics_vec_rf <- function(partial_mouth){
#   # a wrapper function, takes current node as input, and outputs the metrics_vec
#   training_data_temp <- training_data %>% 
#     select(seqn, MS, all_of(partial_mouth)) %>%
#     mutate(seqn = row_number()) %>%
#     as.data.frame()
#   my_task <- makeClassifTask(id = "seqn", data = training_data_temp, 
#                              target = "MS")
#   lrn <- makeLearner("classif.ranger", predict.type = "prob", 
#                      keep.inbag = TRUE)
#   mdl <- train(lrn, my_task)
#   oob_predictions <- getOOBPreds(mdl, my_task)
#   metrics_vec <- 
#     calculate_metrics2(truth = (oob_predictions$data$truth == TRUE) *1, 
#                        pred = oob_predictions$data$prob.TRUE, 
#                        metrics = mtrcs, 
#                        threshold = "auto")
#   return(metrics_vec)
# }

get_metrics_vec_naive <- function(partial_mouth){
  # a wrapper function, takes current node as input, and outputs the metrics_vec
  training_data_temp <- training_data %>% 
    select(matches("ohx\\d\\d(pc|la)[dspa]")) %>%
    select(any_of(partial_mouth)) %>%
    case_definition()
  metrics_vec <- 
    calculate_metrics2(truth = (training_data_raw$MS == TRUE) *1, 
                       pred = (training_data_temp$MS == TRUE) *1, 
                       metrics = mtrcs, 
                       threshold = 0.5)
  metrics_vec["AUC"] <- metrics_vec["Sensitivity"]
  return(metrics_vec)
}

# parameters for evolutionary algorithm
total_n_nodes = 100 # population size, n of tooth combinations to be evaluated
top_n_nodes = 20 # elite offspring, n of highest AUC tooth combinations to
mutate_rate = 0.2 # proportion of combinations that have a mutation
n_mutation = 1 # n_mutation tooth different
recombination_rate = 0.5 # the ratio of (different) teeth selected from one parent combination

max_steps = 50
top_mean_auc_vec = 0 # this is to track the changes in mean AUC of top nodes, so as to stop the iterations
thresholds <- c(2, 5) # convergence criterion, two decreases in the last five observations
stationary_steps = 3 # after the evolutionary algorithm becomes stationary, continue search around the optimal combination for some iterations

solutions <- matrix(data = "", ncol = n_teeth, nrow = total_n_nodes,
                    dimnames = list(NULL, paste("V", 1:n_teeth, sep = "")))
solution_aucs <- rep(0, total_n_nodes)
all_solutions <- tibble(
  as_tibble(solutions)[0, ],
  auc = numeric()) # record all the evaluations during the search.
# This can be use to track Because we are searching in a non-deterministic environment, this will increase the sample size and get better ranks and evaluations for the final results
total_searches_ea <- 0
stationary_count = 0

for(kkkk in 1:total_n_nodes){
  current_node = sort(sample(all_tooth, n_teeth, replace = FALSE))
  cat("Current node:", current_node, "\n")
  partial_mouth <- str_c(prefix, 
                         rep(current_node, each = length(suffix)), 
                         rep(suffix, times = length(current_node)))
  metrics_vec <- get_metrics_vec_naive(partial_mouth)
  current_auc <- metrics_vec["AUC"]
  # current_auc <- 
  #   sapply(X = 1:nboot, simplify = TRUE, 
  #          FUN = function(bbbb){
  #            # rand_ind <- sample(1:nrow(training_data), replace = TRUE)
  #            # train_temp <- training_data[rand_ind,]
  #            res <- 
  #              predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = training_data, 
  #                          input_teeth = partial_mouth, n_cores = n.cores, type ="results")
  #            metrics_vec <- 
  #              calculate_metrics2(truth = training_data_raw[["MS"]]*1, 
  #                                 pred = res$pred_MS, metrics = mtrcs, 
  #                                 threshold = "auto")
  #            return(metrics_vec["AUC"])
  #          })
  # current_auc <- mean(current_auc)
  cat("Current auc:", round(current_auc,4), "\n")
  solutions[kkkk,]<- current_node
  solution_aucs[kkkk] <- current_auc
}
solutions_dup <- duplicated.matrix(solutions, MARGIN = 1)
solutions <- solutions[!solutions_dup,]
solution_aucs <- solution_aucs[!solutions_dup]
solutions_order <- order(solution_aucs, decreasing = TRUE)
solutions <- solutions[solutions_order,]
solution_aucs <- solution_aucs[solutions_order]
all_solutions <- bind_rows(
  all_solutions,
  tibble(as_tibble(solutions), auc = solution_aucs)
)
cat("Current average AUCs for all solutions:", 
    round(mean(solution_aucs), 4), "\n")
cat("Current average AUCs for top solutions:", 
    round(mean(solution_aucs[1:top_n_nodes]), 4), "\n")
top_mean_auc_vec <- c(mean(solution_aucs[1:top_n_nodes]), top_mean_auc_vec)
total_searches_ea = total_searches_ea+total_n_nodes

for(ssss in 1:max_steps){
  parents <- solutions[1:top_n_nodes, ]
  fold <- 2 # increase the size, some recombined solutions after mutation can be still duplicated
  parent1 <- asplit(parents[sample.int(n = top_n_nodes, 
                                       size = fold*total_n_nodes, 
                                       replace = TRUE),], MARGIN = 1)
  parent2 <- asplit(parents[sample.int(n = top_n_nodes, 
                                       size = fold*total_n_nodes, 
                                       replace = TRUE),], MARGIN = 1)
  common_tooth <- mapply(FUN = base::intersect, SIMPLIFY = FALSE,
                         parent1, parent2)
  n_common_tooth <- lengths(common_tooth)
  n_choose1 <- floor((n_teeth - n_common_tooth)*recombination_rate) 
  n_choose2 <- n_teeth - n_common_tooth - n_choose1
  #it may choose 1 tooth less from parent1, 
  #but it does not matter because 1&2 are exchangeable
  distinct_tooth1 <- mapply(FUN = base::setdiff, SIMPLIFY = FALSE,
                            parent1, parent2)
  tooth_from1 <- mapply(FUN = sample, SIMPLIFY = FALSE,
                        distinct_tooth1, n_choose1)
  distinct_tooth2 <- mapply(FUN = base::setdiff, SIMPLIFY = FALSE,
                            parent2, parent1)
  tooth_from2 <- mapply(FUN = sample, SIMPLIFY = FALSE, 
                        distinct_tooth2, n_choose2)
  # recombination
  nodes_to_search <- mapply(FUN = c, SIMPLIFY = FALSE,
                            common_tooth, tooth_from1, tooth_from2)
  # to matrix
  nodes_to_search <- matrix(unlist(nodes_to_search), byrow=TRUE, 
                            nrow = length(nodes_to_search))
  # mutation
  for(ind in sample.int(fold*total_n_nodes, mutate_rate*fold*total_n_nodes))
    nodes_to_search[ind, sample.int(n_teeth, size = n_mutation)] <- 
    sample(setdiff(all_tooth, nodes_to_search[ind, ]), size = n_mutation)
  # order, add parents back, and remove duplicates
  # if we are searching in a deterministic environment, we dont need the evaluate
  # the nodes we have visited. Can create a table to prevent this.
  nodes_to_search <- t(apply(X = nodes_to_search, MARGIN = 1, 
                             simplify = TRUE, FUN = sort, decreasing = FALSE))
  nodes_to_search <- rbind(parents, nodes_to_search)
  nodes_to_search <- unique.matrix(nodes_to_search)
  if(nrow(nodes_to_search) < total_n_nodes){
    cat(str_glue("Rows in nodes_to_search = {nrow(nodes_to_search)} is smaller than total_n_nodes = {total_n_nodes}. May require to increase fold (currently {fold})"))
  } else{
    nodes_to_search <- nodes_to_search[1:total_n_nodes, ]
  }
  
  aucs <- rep(0, nrow(nodes_to_search))
  for(rrrr in 1:nrow(nodes_to_search)){
    node <- nodes_to_search[rrrr, ]
    cat("Current node:", node, "\n")
    if(TRUE){
      # search in a deterministic environment
      existed_node <- find_row_index(record = node, 
                                     mat = as.matrix(select(all_solutions, starts_with("V"))))
    }
    if(exists("existed_node") && existed_node > 0){ # longer && so that the right part will not be assessed
        current_auc <- all_solutions[["auc"]][existed_node]
        cat("Current auc (skipped):", round(current_auc,4), "\n")
    } else {
      partial_mouth <- str_c(prefix, 
                             rep(node, each = length(suffix)), 
                             rep(suffix, times = length(node)))
      metrics_vec <- get_metrics_vec_naive(partial_mouth)
      current_auc <- metrics_vec["AUC"]
      cat("Current auc:", round(current_auc,4), "\n")
    }
    aucs[rrrr] <- current_auc
  }
  total_searches_ea = total_searches_ea+nrow(nodes_to_search)
  
  solutions <- nodes_to_search
  solution_aucs <- aucs
  solutions_order <- order(solution_aucs, decreasing = TRUE)
  solutions <- solutions[solutions_order,]
  solution_aucs <- solution_aucs[solutions_order]
  if(ssss == 1){
    print("Solutions:")
    print(solutions)
    print("Solution AUCs:")
    print(solution_aucs)
  }
  all_solutions <- bind_rows(
    all_solutions,
    tibble(as_tibble(solutions), auc = solution_aucs)
  )
  cat("Current average AUCs for all solutions:", 
      round(mean(solution_aucs), 4), "\n")
  cat("Current average AUCs for top solutions:", 
      round(mean(solution_aucs[1:top_n_nodes]), 4), "\n")
  top_mean_auc_vec <- c(mean(solution_aucs[1:top_n_nodes]), top_mean_auc_vec)
  if(stationary_count == 0){
    if(length(top_mean_auc_vec) >= thresholds[2])
      if(sum(top_mean_auc_vec[1:(thresholds[2]-1)] - 
             top_mean_auc_vec[2:thresholds[2]] <= 0) >= thresholds[1]){
        cat("Evolutionary algorithm converges!\n")
        stationary_count <- 1
      }
  } else {# stationary_count >= 1
    stationary_count = stationary_count + 1
    if(stationary_count >= stationary_steps){
      break # break the ssss loop
    }
  }
}
# results
save.image(file = str_glue("{workdir}/res/EA_naive_{Date}.RData"))
saveRDS(object = all_solutions, file = str_glue("{workdir}/res/EA_naive_record_{Date}.rds"))
saveRDS(object = top_n_nodes, file = str_glue("{workdir}/res/EA_naive_ntop_{Date}.rds"))

cat("Total searched tooth combinations for EA:", total_searches_ea, "\n")
solutions <- all_solutions %>% 
  group_by(across(starts_with("V"))) %>%
  summarise(repeated_n=n(), ave_auc= mean(auc)) %>%
  ungroup() %>%
  arrange(desc(ave_auc)) %>%
  rowwise() %>%
  mutate(regimen = str_c(c_across(starts_with("V")), collapse = "")) %>%
  ungroup() 
print(solutions, n=top_n_nodes)
saveRDS(object = solutions, file = str_glue("{workdir}/res/EA_naive_solutions_{Date}.rds"))

# no recalculation required
aucvec <- solutions$ave_auc
names(aucvec) <- solutions$regimen
aucvec_top <- sort(aucvec, decreasing = TRUE)[1:top_n_nodes]
solutions_interested_auc <- 
  solutions_interested_sens <- aucvec_top
saveRDS(object = solutions_interested_auc, 
        file = str_glue("{workdir}/res/EA_naive_aucvec_{Date}.rds"))
saveRDS(object = solutions_interested_sens, 
        file = str_glue("{workdir}/res/EA_naive_senvec_{Date}.rds"))

# # recalculate the aucs of the top solutions
# n_recal = 3
# aucvec <- solutions$ave_auc
# names(aucvec) <- solutions$regimen
# aucvec_top <- sort(aucvec, decreasing = TRUE)[1:50]
# solutions_interested <- 
#   c(names(aucvec_top))
# solutions_interested_sens <- solutions_interested_auc <- 
#   do.call(rbind, 
#           replicate(n = n_recal, simplify = FALSE,
#                     setNames(object = rep(0, length(solutions_interested)), 
#                              nm = solutions_interested)))
# cat("Recalculate all the solutions:\n")
# for(recal_ind in 1:n_recal){
#   for(sol in solutions_interested){
#     current_node = str_sub(sol, start = (1:n_teeth)*2-1, end = (1:n_teeth)*2)
#     cat("Current node:", current_node, "\n")
#     partial_mouth <- str_c(prefix, 
#                            rep(current_node, each = length(suffix)), 
#                            rep(suffix, times = length(current_node)))
#     metrics_vec <- get_metrics_vec_naive(partial_mouth)
#     solutions_interested_auc[recal_ind, sol] <- metrics_vec["AUC"]
#     solutions_interested_sens[recal_ind, sol] <- metrics_vec["Sensitivity"]
#   }
# }
# solutions_interested_auc <- colMeans(solutions_interested_auc)
# solutions_interested_sens <- colMeans(solutions_interested_sens)
# saveRDS(object = solutions_interested_auc, file = str_glue("{workdir}/res/EA_naive_aucvec_{Date}.rds"))
# saveRDS(object = solutions_interested_sens, file = str_glue("{workdir}/res/EA_naive_senvec_{Date}.rds"))




