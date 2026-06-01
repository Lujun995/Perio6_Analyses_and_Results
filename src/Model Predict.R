source("predict.gen.dependencies.R") # assume the working directory is where Model Predict.R exists
source("case_def.R")
source("predict.gen.R")
n.cores = 5

# read in the dataset
data <- read_csv("...")
data_normalized <- data %>%
  rename_with(tolower) %>%
  rename("seqn" = "studyid") %>%
  filter(if_all(starts_with("ohx"), ~ .x >= 0)) %>%
  mutate(across(.cols = starts_with("ohx"),
                .fns = \(x) ifelse(x > 6, 6, x)))
data_normalized <- case_def(d = data_normalized) %>%
  mutate(MS = (mod_12|sev_12))


all_tooth <- sprintf("%02d", c(seq(2,15), seq(18,31))) 
Ramfjord <- c("03", "09", "12", "19", "25", "28")
half_mouth <- c("02", "03", "04", "05", "06", "07", "08",
                "18", "19", "20", "21", "22", "23", "24")
best_teeth_1 <- c("03", "14", "18", "25", "30", "31")
suffix = c(str_c("pc", c("d", "m", "s", "p", "l", "a"), sep = ""), 
           str_c("la", c("d", "m", "s", "p", "l", "a"), sep = "")) # pcd, ..., lad, ...
all_ohx <- str_c("ohx", 
                 rep(all_tooth, each = length(suffix)), 
                 rep(suffix, times = length(all_tooth)))
input_teeth <- str_c("ohx", 
                     rep(best_teeth_1, each = length(suffix)), 
                     rep(suffix, times = length(best_teeth_1)))

cal_dist <- data_normalized %>% 
  pivot_longer(cols = starts_with("ohx"), 
               names_to = "tooth_id", 
               values_to = "value") %>%
  mutate(value = factor(value, levels = 0:6)) %>%
  mutate(tooth = str_extract(tooth_id, "\\d\\d"),
         test = str_extract(tooth_id, "pc|la"),
         site = str_extract(tooth_id, "[dspaml]$")) %>%
  mutate(site = case_match(site, "d" ~ "DB" , "s" ~ "MB", "m" ~ "B",
                           "p" ~ "DL", "a" ~ "ML", "l" ~ "L"),
         site = factor(site, levels = c("DB", "B", "MB", "DL", "L", "ML"))) %>%
  mutate(from = "Temple-C") %>%
  filter(test == "la") %>%
  ggplot(aes(x = tooth, fill = value)) +
  geom_bar(color = NA, position = position_stack(reverse = FALSE)) + 
  geom_vline(xintercept = length(all_tooth)/2+0.5, 
             color = "grey20", linetype = "dashed") + 
  facet_grid(.~site) + 
  theme_bw() + 
  scale_fill_brewer(name = "CAL/mm", palette = "YlGnBu") +
  labs(x = "Tooth Number", y = "Number of Subjects") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
        #element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        panel.grid = element_blank())
cal_dist
ggsave(plot = cal_dist, filename = "../res/cal_dist.jpeg", 
       width = 12, height = 2)



# function parameters
cal_mdl <- readRDS("Data/Models/nhanes_cal_genmdl_ord.rds")
ppd_mdl <- readRDS("Data/Models/nhanes_ppd_genmdl_ord.rds")

# For Temple U data
data_normalized_mat <- data_normalized %>%
  column_to_rownames(var = "seqn") %>%
  select(all_of(all_ohx)) %>%
  as.matrix()
new_data <- data_normalized[1:100, ]

# multi-thread version
t <- proc.time()
res <- predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = new_data, 
                   input_teeth = input_teeth, 
                   n_cores = n.cores, type ="results")
print(proc.time() - t) # 35 sec
data_normalized$seqn <- as.character(data_normalized$seqn)
res <- res %>%
  left_join(x= ., y = data_normalized, by = "seqn") %>%
  select(seqn, pred_MS, MS, pred_sev_12, sev_12) %>%
  rename("truth_MS" = "MS", "truth_sev_12" = "sev_12") %>%
  print(n = 100)
calculate_metrics2(truth = res$truth_MS, 
                   pred = res$pred_MS,
                   threshold = 0.2)
calculate_metrics2(truth = res$truth_MS, 
                   pred = res$pred_MS,
                   threshold = "auto")


# # single-thread version when the sample size is small
# t <- proc.time()
# res <- predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = new_data, 
#                    input_teeth = input_teeth,
#                    n_cores = 1, type ="results")
# print(proc.time() - t) # 42 sec. single-thread version may be faster when n is small
# res %>%
#   left_join(x= ., y = data_complete, by = "seqn") %>%
#   select(seqn, pred_MS, MS, pred_sev_12, sev_12) %>%
#   rename("truth_MS" = "MS", "truth_sev_12" = "sev_12") %>%
#   print()

# # single person version. Even if n_core = 5 by default, it will run using a single thread
# t <- proc.time()
# res <- predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, 
#                    new_data = new_data[1, , drop = FALSE], 
#                    input_teeth = input_teeth, 
#                    type ="results")
# print(proc.time() - t) # 4.7 sec
# res %>%
#   left_join(x= ., y = data_complete, by = "seqn") %>%
#   select(seqn, pred_MS, MS, pred_sev_12, sev_12) %>%
#   rename("truth_MS" = "MS", "truth_sev_12" = "sev_12") %>%
#   print()

# # missing data in the input teeth
# if(file.exists("../res/aric_extended.rds")){
#   data_extended <- readRDS("../res/aric_extended.rds")
# } else {
#   data_extended <- readRDS("../res/nhanes_extended.rds")
# }
# new_data <- data_extended %>%
#   mutate(across(.cols = starts_with("ohx"),
#                 .fns = \(x) case_match(x, "0"~0, "1"~1, "2"~2, "3"~3,
#                                        "4"~4, "5"~5, "6/6+"~6))) %>%
#   column_to_rownames(var = "seqn") %>%
#   select(all_of(all_ohx)) %>%
#   as.matrix()
# new_data <- new_data[1:10, ]
# 
# # some test:
# res_temp <- predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = new_data,
#                         input_teeth = input_teeth,
#                         n_cores = 5, type ="results")
# res_temp <- predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = new_data, 
#                         input_teeth = input_teeth, 
#                         n_cores = 1, type ="data")
# temp1 <- as.matrix(res_temp[1:10, input_teeth]) # Unoberved values at input teeth were imputed
# temp2 <- new_data[, input_teeth]
# # Observed values at input teeth were not imputed
# all.equal(c(temp1)[!is.na(c(temp2))], c(temp2)[!is.na(c(temp2))]) # TRUE
# # Values at non-input teeth were always generated
# for(varr in sample(setdiff(colnames(new_data), input_teeth), 10)){
#   cat(varr, ":\n")
#   matrix(res_temp[[varr]] == rep(new_data[,varr], 50), 
#          ncol = 10, byrow = TRUE) %>%
#     apply(X = ., FUN = table, simplify = FALSE, useNA = "ifany", MARGIN = 2) %>%
#     print()
# }


# # complete cases should be faster than missing input teeth cases
# new_data <- data_extended %>%
#   mutate(across(.cols = starts_with("ohx"),
#                 .fns = \(x) case_match(x, "0"~0, "1"~1, "2"~2, "3"~3,
#                                        "4"~4, "5"~5, "6/6+"~6))) %>%
#   column_to_rownames(var = "seqn") %>%
#   select(all_of(all_ohx)) %>%
#   as.matrix()
# new_data <- new_data[1:10, ]
# t <- proc.time()
# res_temp <- predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = new_data, 
#                         input_teeth = input_teeth, 
#                         n_cores = 1, type ="data")
# proc.time() - t # 50s
# new_data <- data_complete %>%
#   mutate(across(.cols = starts_with("ohx"),
#                 .fns = \(x) case_match(x, "0"~0, "1"~1, "2"~2, "3"~3,
#                                        "4"~4, "5"~5, "6/6+"~6))) %>%
#   column_to_rownames(var = "seqn") %>%
#   select(all_of(all_ohx)) %>%
#   as.matrix()
# new_data <- new_data[1:10, ]
# t <- proc.time()
# res_temp <- predict.gen(cal_mdl = cal_mdl, ppd_mdl = ppd_mdl, new_data = new_data, 
#                         input_teeth = input_teeth, 
#                         n_cores = 1, type ="data")
# proc.time() - t # 40s ~ 50/28*22 = 39
  

# }



