library(tidyverse)
library(parallel)
library(ordinalNet)
workdir = "/home/guanwh/zhan7474/dental_health"
Date = "0621" # for version control
load(str_glue("{workdir}/res/EA_{Date}.RData"))
solutions <- readRDS(file = str_glue("{workdir}/res/EA_solutions_{Date}.rds"))
source(str_glue("{workdir}/src/case_def.R"))
source(str_glue("{workdir}/src/predict.gen.dependencies.R")) 
source(str_glue("{workdir}/src/predict.gen.R"))
source(str_glue("{workdir}/src/calculate_metrics2.R"))

cl <- makeCluster(n.cores)

# recalculate the aucs of the top solutions
n_recal = 5
aucvec <- solutions$ave_auc
names(aucvec) <- solutions$regimen
aucvec_100 <- sort(aucvec, decreasing = TRUE)[1:50]
Ramfjord <- c("03", "09", "12", "19", "25", "28")
# half_mouth <- c("02", "03", "04", "05", "06", "07", "08",
#                 "18", "19", "20", "21", "22", "23", "24")
solutions_interested <- 
  c(names(aucvec_100), "030814192530", paste(Ramfjord, collapse = ""))
solutions_interested_sens <- solutions_interested_auc <- 
  do.call(rbind, 
          replicate(n = n_recal, simplify = FALSE,
                    setNames(object = rep(0, length(solutions_interested)), 
                             nm = solutions_interested)))
cat("Recalculate all the solutions:\n")
for(recal_ind in 1:n_recal){
  for(sol in solutions_interested){
    current_node = str_sub(sol, start = (1:6)*2-1, end = (1:6)*2)
    cat("Current node:", current_node, "\n")
    partial_mouth <- str_c(prefix, 
                           rep(current_node, each = length(suffix)), 
                           rep(suffix, times = length(current_node)))
    res <- 
      predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = training_data, 
                  input_teeth = partial_mouth, n_cores = n.cores, type ="results", cl=cl)
    metrics_vec <- 
      calculate_metrics2(truth = training_data_raw[["MS"]]*1, 
                         pred = res$pred_MS, metrics = mtrcs, 
                         threshold = "auto")
    solutions_interested_auc[recal_ind, sol] <- metrics_vec["AUC"]
    solutions_interested_sens[recal_ind, sol] <- metrics_vec["Sensitivity"]
  }
}
colSDs <- function(x, na.rm = FALSE) {
  apply(x, 2, sd, na.rm = na.rm)
}
solutions_interested_auc_sd <- colSDs(solutions_interested_auc)
solutions_interested_auc <- colMeans(solutions_interested_auc)
solutions_interested_sens_sd <- colSDs(solutions_interested_sens)
solutions_interested_sens <- colMeans(solutions_interested_sens)
saveRDS(object = solutions_interested_auc, file = str_glue("{workdir}/res/EA_aucvec_{Date}.rds"))
saveRDS(object = solutions_interested_auc_sd, file = str_glue("{workdir}/res/EA_aucsd_{Date}.rds"))
saveRDS(object = solutions_interested_sens, file = str_glue("{workdir}/res/EA_senvec_{Date}.rds"))
saveRDS(object = solutions_interested_sens_sd, file = str_glue("{workdir}/res/EA_sensd_{Date}.rds"))

stopCluster(cl)