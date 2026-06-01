workdir = "/home/guanwh/zhan7474/dental_health"
Date = "0621" # for version control
n.cores = 100

library(tidyverse)
library(parallel)
library(ordinalNet)
set.seed(seed= 1233455)
source(str_glue("{workdir}/src/case_def.R"))
source(str_glue("{workdir}/src/predict.gen.dependencies.R")) 
source(str_glue("{workdir}/src/predict.gen.R"))
source(str_glue("{workdir}/src/calculate_metrics2.R"))

prefix = "ohx"
suffix = c(str_c("pc", c("d", "m", "s", "p", "l", "a"), sep = ""), 
           str_c("la", c("d", "m", "s", "p", "l", "a"), sep = "")) # pcd, ..., lad, ...
all_tooth <- sprintf("%02d", c(seq(2,15), seq(18,31))) # all teeth in NHANES
all_ohx <- str_c(prefix, 
                 rep(all_tooth, each = length(suffix)), 
                 rep(suffix, times = length(all_tooth)))
decision_threshold = 0.2

training_data_raw <- readRDS(str_glue("{workdir}/res/nhanes_complete.rds"))
training_data <- training_data_raw %>%
  mutate(across(.cols = starts_with("ohx"),
                .fns = \(x) case_match(x, "0"~0, "1"~1, "2"~2, "3"~3,
                                       "4"~4, "5"~5, "6/6+"~6))) %>%
  column_to_rownames(var = "seqn") %>%
  select(starts_with("ohx")) %>%
  as.matrix()

valid_data_raw <- readRDS(str_glue("{workdir}/res/aric_complete.rds"))
valid_data <- valid_data_raw %>%
  mutate(across(.cols = starts_with("ohx"),
                .fns = \(x) case_match(x, "0"~0, "1"~1, "2"~2, "3"~3,
                                       "4"~4, "5"~5, "6/6+"~6))) %>%
  column_to_rownames(var = "seqn") %>%
  select(starts_with("ohx")) %>%
  as.matrix()

valid_data_raw <- origins_complete
valid_data <- valid_data_raw %>%
  mutate(across(.cols = starts_with("ohx"),
                .fns = \(x) case_match(x, "0"~0, "1"~1, "2"~2, "3"~3,
                                       "4"~4, "5"~5, "6/6+"~6))) %>%
  column_to_rownames(var = "seqn") %>%
  select(starts_with("ohx")) %>%
  as.matrix()



#rand_index <- sapply(X = 1:nboot, 
#                     FUN = \(x) return(sample.int(
#                       n = nrow(training_data),
#                       size = round(train_prop * nrow(training_data)),
#                       replace = FALSE, prob = NULL, useHash = FALSE)))

mtrcs = c("Prevalence", "AUC", "Sensitivity","Specificity", "F1-score", 
          "Diagnostic Odds Ratio")
cal_mdl <- readRDS(str_glue("{workdir}/res/nhanes_cal_genmdl_ord.rds"))
ppd_mdl <- readRDS(str_glue("{workdir}/res/nhanes_ppd_genmdl_ord.rds"))

res_tot <-  list()
for(i in 1:3){
  current_node = c("03", "15", "19", "25", "30", "31")
  partial_mouth <- str_c(prefix, 
                         rep(current_node, each = length(suffix)), 
                         rep(suffix, times = length(current_node)))
  res1 <- 
    predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = training_data, 
                input_teeth = partial_mouth, n_cores = n.cores, type ="results")
  res_tot[[length(res_tot)+1]] <- res1
  metrics_vec1 <- 
    calculate_metrics2(truth = training_data_raw[["MS"]]*1, 
                       pred = res1$pred_MS, metrics = mtrcs, 
                       threshold = 0.2)
  metrics_vec1
  res2 <- 
    predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = valid_data, 
                input_teeth = partial_mouth, n_cores = n.cores, type ="results")
  res_tot[[length(res_tot)+1]] <- res2
  metrics_vec2 <- 
    calculate_metrics2(truth = valid_data_raw[["MS"]]*1, 
                       pred = res2$pred_MS, metrics = mtrcs, 
                       threshold = 0.2)
  metrics_vec2
  
  current_node = c("03", "14", "18", "25", "30", "31")
  partial_mouth <- str_c(prefix, 
                         rep(current_node, each = length(suffix)), 
                         rep(suffix, times = length(current_node)))
  res3 <- 
    predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = training_data, 
                input_teeth = partial_mouth, n_cores = n.cores, type ="results")
  res_tot[[length(res_tot)+1]] <- res3
  metrics_vec3 <- 
    calculate_metrics2(truth = training_data_raw[["MS"]]*1, 
                       pred = res3$pred_MS, metrics = mtrcs, 
                       threshold = 0.2)
  metrics_vec3
  res4 <- 
    predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = valid_data, 
                input_teeth = partial_mouth, n_cores = n.cores, type ="results")
  res_tot[[length(res_tot)+1]] <- res4
  metrics_vec4 <- 
    calculate_metrics2(truth = valid_data_raw[["MS"]]*1, 
                       pred = res4$pred_MS, metrics = mtrcs, 
                       threshold = 0.2)
  metrics_vec4
  
  current_node = c("03", "14", "18", "25", "30")
  partial_mouth <- str_c(prefix, 
                         rep(current_node, each = length(suffix)), 
                         rep(suffix, times = length(current_node)))
  res5 <- 
    predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = training_data, 
                input_teeth = partial_mouth, n_cores = n.cores, type ="results")
  res_tot[[length(res_tot)+1]] <- res5
  metrics_vec5 <- 
    calculate_metrics2(truth = training_data_raw[["MS"]]*1, 
                       pred = res5$pred_MS, metrics = mtrcs, 
                       threshold = 0.2)
  metrics_vec5
  res6 <- 
    predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = valid_data, 
                input_teeth = partial_mouth, n_cores = n.cores, type ="results")
  res_tot[[length(res_tot)+1]] <- res6
  metrics_vec6 <- 
    calculate_metrics2(truth = valid_data_raw[["MS"]]*1, 
                       pred = res6$pred_MS, metrics = mtrcs, 
                       threshold = 0.2)
  metrics_vec6
}

for(i in 1:length(res_tot)){
  if(i%%7 == 1)
    saveRDS(object = res_tot[[i]], 
            file = str_glue("{workdir}/res/nhanes_6t_set2_rep{i%/%7}.rds"))
  if(i%%7 == 2)
    saveRDS(object = res_tot[[i]], 
            file = str_glue("{workdir}/res/aric_6t_set2_rep{i%/%7}.rds"))
  if(i%%7 == 3)
    saveRDS(object = res_tot[[i]], 
            file = str_glue("{workdir}/res/nhanes_6t_set1_rep{i%/%7}.rds"))
  if(i%%7 == 4)
    saveRDS(object = res_tot[[i]], 
            file = str_glue("{workdir}/res/aric_6t_set1_rep{i%/%7}.rds"))
  if(i%%7 == 5)
    next #   res4_e
  if(i%%7 == 6)
    saveRDS(object = res_tot[[i]], 
            file = str_glue("{workdir}/res/nhanes_5t_set1_rep{i%/%7}.rds"))
  if(i%%7 == 0)
    saveRDS(object = res_tot[[i]], 
            file = str_glue("{workdir}/res/aric_5t_set1_rep{i%/%7-1}.rds"))
}

library(tidyverse)
library(parallel)
library(ordinalNet)
workdir = "/home/guanwh/zhan7474/dental_health"
load(str_glue("{workdir}/res/valid.RData"))


valid_data_raw_extend <- readRDS(str_glue("{workdir}/res/aric_extended.rds"))
aucres_6 <- readRDS(str_glue("../res/EA_aucvec_0621.rds"))
aucres_5 <- readRDS(str_glue("../res/EA_aucvec_0711.rds"))
aucres_4 <- readRDS(str_glue("../res/EA_aucvec_0712.rds"))
tooth_combintations <- c(names(sort(aucres_6, T))[1:20], 
                         # top 6 teeth combinations are almost always better than 5 and 4
                         names(sort(aucres_5, T))[1:20],
                         names(sort(aucres_4, T))[1:20])
seqn_no <- valid_data_extend_list <- vector("list", length(tooth_combintations))
names(seqn_no) <- names(valid_data_extend_list) <- tooth_combintations
for(tc in tooth_combintations){
  n_teeth = nchar(tc)/2
  current_node = str_sub(tc, start = (1:n_teeth)*2-1, end = (1:n_teeth)*2)
  partial_mouth <- str_c(prefix, 
                         rep(current_node, each = length(suffix)), 
                         rep(suffix, times = length(current_node)))
  valid_data_extend <- valid_data_raw_extend %>%
    filter(!(seqn %in% unlist(seqn_no)))
  for(tid in current_node){
    tid_var <- str_c(prefix, 
                     rep(tid, each = length(suffix)), 
                     rep(suffix, times = length(tid)))
    valid_data_extend <- valid_data_extend %>%
      filter(if_any(all_of(tid_var), ~ !is.na(.)))
  }
  if(nrow(valid_data_extend) == 0)
    next
  valid_data_extend_list[[tc]] <- 
    valid_data_extend %>%
    mutate(across(.cols = starts_with("ohx"),
                  .fns = \(x) case_match(x, "0"~0, "1"~1, "2"~2, "3"~3,
                                         "4"~4, "5"~5, "6/6+"~6))) %>%
    column_to_rownames(var = "seqn") %>%
    select(starts_with("ohx")) %>%
    as.matrix() # n = 1741
  cat(tc,":\n")
  print(nrow(valid_data_extend_list[[tc]]))
  seqn_no[[tc]] <- rownames(valid_data_extend_list[[tc]])
}
res_e <- NULL
cl <- makeCluster(n.cores)
for(jjjj in 1:length(valid_data_extend_list)){
  if(is.null(valid_data_extend_list[[jjjj]]))
    next
  tc <- names(valid_data_extend_list)[jjjj]
  n_teeth = nchar(tc)/2
  current_node = str_sub(tc, start = (1:n_teeth)*2-1, end = (1:n_teeth)*2)
  partial_mouth <- str_c(prefix, 
                         rep(current_node, each = length(suffix)), 
                         rep(suffix, times = length(current_node)))
  valid_data_extend <- valid_data_extend_list[[jjjj]]
  nr <- nrow(valid_data_extend)
  res_e <- bind_rows(res_e,
                     predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, 
                                 new_data = valid_data_extend, 
                                 input_teeth = partial_mouth, 
                                 n_cores = if(nr <= 5) 1 else n.cores,
                                 cl = cl,
                                 type ="results"))
}
stopCluster(cl)
truth = valid_data_raw_extend[["MS"]][match(res_e$seqn, 
                                            valid_data_raw_extend[["seqn"]])]
metrics_vec_e <- 
  calculate_metrics2(truth = truth*1, 
                     pred = res_e$pred_MS, metrics = mtrcs, 
                     threshold = 0.2)
metrics_vec_e
saveRDS(object = res_e, 
          file = str_glue("{workdir}/res/aric_e.rds"))
saveRDS(object = seqn_no, 
        file = str_glue("{workdir}/res/aric_e_seqn_no.rds"))
# Prevalence                   AUC           Sensitivity
# 0.5695811                    NA             0.9364428
# Specificity              F1-score Diagnostic Odds Ratio
# 0.8760897             0.9187830           104.1736505
# truth = valid_data_raw_extend$MS[match(unlist(seqn_no[1]), 
#                                        valid_data_raw_extend[["seqn"]])]
# pred = res_e$pred_MS[match(unlist(seqn_no[1]), 
#                            res_e$seqn)]
# metrics_vec_e <- 
#   calculate_metrics2(truth = truth*1, 
#                      pred = pred, metrics = mtrcs, 
#                      threshold = 0.2)
# metrics_vec_e
# # Prevalence                   AUC           Sensitivity
# # 0.5695811                    NA             0.9364428
# # Specificity              F1-score Diagnostic Odds Ratio
# # 0.8760897             0.9187830           104.1736505


save.image(str_glue("{workdir}/res/valid.RData"))

# current_node = c("03", "15", "19", "25", "30", "31")
# NHANES_complete
#           Prevalence                   AUC           Sensitivity
#            0.1851422                    NA             0.9102564
#          Specificity              F1-score Diagnostic Odds Ratio
#            0.9751949             0.9001585           398.7591837
# ARIC complete
#           Prevalence                   AUC           Sensitivity
#            0.4795396                    NA             0.9125683
#          Specificity              F1-score Diagnostic Odds Ratio
#            0.9014423             0.9014845            95.4649390
# current_node = c("03", "14", "18", "25", "30", "31")
# NHANES_complete
#           Prevalence                   AUC           Sensitivity
#            0.1892049                    NA             0.9102564
#          Specificity              F1-score Diagnostic Odds Ratio
#            0.9702339             0.8902821           330.6088435
# ARIC complete
#           Prevalence                   AUC           Sensitivity
#            0.4693095                    NA             0.9071038
#          Specificity              F1-score Diagnostic Odds Ratio
#            0.9158654             0.9058663           106.2957983
# ARIC extend
#           Prevalence                   AUC           Sensitivity
#            0.5301551                    NA             0.9328859
#          Specificity              F1-score Diagnostic Odds Ratio
#            0.8949233             0.9179967           118.3842697
# current_node = c("03", "14", "18", "25", "30")
# NHANES_complete
#           Prevalence                   AUC           Sensitivity
#            0.1897853                    NA             0.9006410
#          Specificity              F1-score Diagnostic Odds Ratio
#            0.9673990             0.8794992           268.9796634
# ARIC complete
#           Prevalence                   AUC           Sensitivity
#            0.4757033                    NA             0.9234973
#          Specificity              F1-score Diagnostic Odds Ratio
#            0.9182692             0.9159892           135.6260504




