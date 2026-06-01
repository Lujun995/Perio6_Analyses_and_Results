library(tidyverse)
workdir = ".."
n.cores = 100
top_n_nodes = 20
prefix = "ohx"
suffix = c(str_c("pc", c("d", "m", "s", "p", "l", "a"), sep = ""), 
           str_c("la", c("d", "m", "s", "p", "l", "a"), sep = "")) # pcd, ..., lad, ...
training_data_raw <- readRDS(str_glue("{workdir}/res/nhanes_complete.rds"))
training_data <- training_data_raw %>%
  mutate(across(.cols = starts_with("ohx"),
                .fns = \(x) case_match(x, "0"~0, "1"~1, "2"~2, "3"~3,
                                       "4"~4, "5"~5, "6/6+"~6))) %>%
  column_to_rownames(var = "seqn") %>%
  select(starts_with("ohx")) %>%
  as.matrix()
mtrcs = c("Prevalence", "AUC", "Sensitivity","Specificity", "F1-score", 
          "Diagnostic Odds Ratio")
cal_mdl <- readRDS(str_glue("{workdir}/res/nhanes_cal_genmdl_ord.rds"))
ppd_mdl <- readRDS(str_glue("{workdir}/res/nhanes_ppd_genmdl_ord.rds"))

n_teeth = 5
Date = "0711"

lines <- readLines(str_glue("{workdir}/res/EA_{Date}.out"))
solutions <- list()
i = 1
while (i <= length(lines)) {
  line <- lines[i]
  if(str_detect(line, "Current node:")) {
    node_match <- str_match(line, "Current node:\\s*(.*)$")[, 2]
    node_vals <- str_split(str_trim(node_match), "\\s+")[[1]]
    entry <- as.list(node_vals)
    names(entry) <- paste0("V", seq_along(entry))
    
    if (i + 1 <= length(lines) && str_detect(lines[i + 1], "^Current auc:")) {
      auc_val <- as.numeric(str_match(lines[i + 1], "Current auc:\\s*([0-9\\.]+)")[, 2])
      entry$auc <- auc_val
      i <- i + 1  # Skip the auc line in next iteration
    } else {
      entry$auc <- NA
    }
    solutions[[length(solutions) + 1]] <- entry
  }
  i <- i + 1
}
all_solutions <- tidyr::drop_na(bind_rows(solutions))

save.image(file = str_glue("{workdir}/res/EA_{Date}.RData"))
saveRDS(object = all_solutions, file = str_glue("{workdir}/res/EA_record_{Date}.rds"))
saveRDS(object = top_n_nodes, file = str_glue("{workdir}/res/EA_ntop_{Date}.rds"))

solutions <- all_solutions %>% 
  group_by(across(starts_with("V"))) %>%
  summarise(repeated_n=n(), ave_auc= mean(auc)) %>%
  ungroup() %>%
  arrange(desc(ave_auc)) %>%
  rowwise() %>%
  mutate(regimen = str_c(c_across(starts_with("V")), collapse = "")) %>%
  ungroup() 
print(solutions, n=top_n_nodes)
saveRDS(object = solutions, file = str_glue("{workdir}/res/EA_solutions_{Date}.rds"))
