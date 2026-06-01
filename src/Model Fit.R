source("../src/fit.gen.R")

# data prepossessing
all_tooth <- sprintf("%02d", c(seq(2,15), seq(18,31))) 
suffix = c(str_c("pc", c("d", "m", "s", "p", "l", "a"), sep = ""), 
           str_c("la", c("d", "m", "s", "p", "l", "a"), sep = "")) # pcd, ..., lad, ...
all_ohx <- str_c("ohx", 
                 rep(all_tooth, each = length(suffix)), 
                 rep(suffix, times = length(all_tooth)))
nhanes_complete <- readRDS("../res/nhanes_complete.rds")
training_data <- nhanes_complete %>%
  mutate(across(.cols = starts_with("ohx"), 
                .fns = \(x) case_match(x, "0"~0, "1"~1, "2"~2, "3"~3, 
                                       "4"~4, "5"~5, "6/6+"~6))) %>%
  column_to_rownames(var = "seqn") %>%
  select(all_of(all_ohx)) %>%
  as.matrix()

# set upper limit and check for leaping levels
# need more data normalization here for model training and prediction:
# the upper limit of ohx measurements:
# the upper limits should be determined in the training dataset 
# and applied to test dataset because they need to use the same set of parameters
# overall:
#   consider using 99.9% as the limit. 
#   threshold = ceiling(quantile(0.999)) will result in less than 0.1% of samples 
#   in the highest class
upper_limit <- ceiling(apply(training_data, MARGIN = 2, 
                             FUN = quantile, prob = 0.999))
min(upper_limit)
rowSums(t(training_data) >= upper_limit)
table(rowSums(t(training_data) >= upper_limit))
# however, for interdental measures:
#   upper_limit of CAL should be at least 6/6+
#   upper_limit of PPD should be at least 5/5+ 
#   (both of them can be at most 6/6+ because of our previous normalization)
upper_limit[upper_limit < 6 & str_detect(names(upper_limit), "ohx\\d\\dla[dspa]")]
upper_limit[str_detect(names(upper_limit), "ohx\\d\\dla[dspa]")] <- 6
upper_limit[upper_limit < 5 & str_detect(names(upper_limit), "ohx\\d\\dpc[dspa]")]
upper_limit[str_detect(names(upper_limit), "ohx\\d\\dpc[dspa]")] <- 5
rowSums(t(training_data) >= upper_limit)
# then, check if there are "jumping" levels:
#   for example, we have observations of 4 and 6, but not 5.
#   for ordinal regression, we have to set 6 to 5/5+ in this case.
#   if we only have observations 0-4, but not 5 or 6, it is okay
for(i in 1:6){
  # check unused upper_limit
  var_index <- (rowSums(t(training_data) >= upper_limit) == 0)
  if(sum(var_index) == 0) break
  if(i == 1) stopifnot(identical(names(upper_limit), names(var_index)))
  # print(var_index[var_index])
  upper_limit[var_index] <- upper_limit[var_index] - 1
}
for(i in 1:3){
  # check leaping levels
  var_index <- (rowSums(t(training_data) == pmax(upper_limit - i, 0)) == 0)
  if(sum(var_index) > 0) cat("var has a leaping level:", names(var_index[var_index]))
  upper_limit[var_index] <- upper_limit[var_index] - 1
}
# apply the upper_limit and normalize the data
mean(training_data != t(pmin(t(training_data), upper_limit))) 
# employing this upper limit only changes 0.065% of observations
training_data <- t(pmin(t(training_data), upper_limit))
# While we use different limits for different varaibles, note that, in fit.gen,
# all the variables will consistently have 0:ohx_max levels. For the unused levels,
# the probability is just 0.


# fit the model
if(!file.exists("../res/nhanes_ohx_genmdl_ord.rds")){
  gen_mdl <- fit.gen(data = training_data)
  saveRDS(gen_mdl, file = "../res/nhanes_ohx_genmdl_ord.rds")
}
if(!file.exists("../res/nhanes_cal_genmdl_ord.rds")){
  gen_mdl <- fit.gen(data = training_data[, str_subset(all_ohx, "ohx\\d\\dla[dspaml]")])
  saveRDS(gen_mdl, file = "../res/nhanes_cal_genmdl_ord.rds")
}
if(!file.exists("../res/nhanes_ppd_genmdl_ord.rds")){
  gen_mdl <- fit.gen(data = training_data[, str_subset(all_ohx, "ohx\\d\\dpc[dspaml]")])
  saveRDS(gen_mdl, file = "../res/nhanes_ppd_genmdl_ord.rds")
}

