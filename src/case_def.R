case_def <- function(
    d, id = "seqn", var_name = "^(ohx)(\\d\\d)(pc|la)([dmspla])$", 
    tooth_num_digit = 4:5, measurement_digit = 6:7, site_digit = 8, 
    PPD = "pc", CAL = "la", n_teeth_total = 28,
    Buccal = "m", Mesiobuccal = "s", Distobuccal = "d",
    Lingual = "l", Mesiolingual = "a", Distolingual = "p",
    detectable_threshold = 0){
    # this function will return a dataframe with these new columns:
    # ps_18, ps_12: periodontitis status based on 2018 EFP/AAP or 2012 CDC/AAP case definitions
    # mod_18, mod_12: moderate periodontitis cases based on the two case definitions
    # sev_18, sev_12: severe periodontitis cases based on the two case definitions
    # n_NA_tooth_cal: number of missing teeth (with all CAL and PD sites are missing)
    # n_NA_site_cal: number of missing CAL sites across present teeth (excluding lost teeth)
    
    # 2018 EFP/AAP Case Definitions (Tonetti et al., 2018, J Clin Periodontol):
    #   *  Pseudo-code for 2018 EFP/AAP Case Definition
    # +  interdental CAL is detectable (for tooth "x"):
    #   -  max_across_sites_at_one_tooth(DF_CAL, MF_CAL, DL_CAL, ML_CAL) >= 1 
    # -  comment: Any of CALs interdental site is detectable at tooth "x"
    # +  Interdental CAL is detectable at ≥2 non-adjacent teeth:
    #   -  count_across_teeth(interdental CAL is detectable) >= 3 | 
    #   (count_across_teeth(interdental CAL is detectable) == 2 & teeth are non-adjacent)
    # -  comment: if we have at least 3 teeth with detectable interdental CAL, at least two teeth of them have to be non-adjacent
    # +  Buccal/Oral CAL ≥ 3 mm with PPD > 3 mm at ≥2 teeth
    # -  count_across_teeth(B_CAL >= 3 & B_PD >= 4) >= 2
    # +  a periodontitis case: Interdental CAL is detectable at ≥2 non-adjacent teeth OR Buccal or oral CAL ≥3 mm with pocketing >3 mm is detectable at ≥2 teeth 
    # -  Interdental CAL is detectable at ≥2 non-adjacent teeth | Buccal/Oral CAL ≥ 3 mm with PPD > 3 mm at ≥2 teeth
    # +  moderate periodontitis cases:
    #   -  a periodontitis case & max_across_teeth(max_across_sites(DF_CAL, MF_CAL, DL_CAL, ML_CAL))  == 3 or 4
    # +  severe periodontitis cases: 
    #   -  a periodontitis case & max_across_teeth(max_across_sites(DF_CAL, MF_CAL, DL_CAL, ML_CAL)) >= 5
    # 
    # 2012 CDC/AAP Case Definition (Eke et al., 2012, J Periodontol):
    #   *  Pseudo-code for 2012 CDC/AAP Case Definition
    # +  Mild periodontitis: ≥2 interproximal sites with AL ≥3 mm, and ≥2 interproximal sites with PD ≥4 mm (not on same tooth) or one site with PD ≥5 mm
    # -  (count_across_teeth(max_across_sites_at_one_tooth(DF_CAL, MF_CAL, DL_CAL, ML_CAL) >= 3) >= 2) &
    #   ( (count_across_teeth(max_across_sites_at_one_tooth(DF_PD, MF_PD, DL_PD, ML_PD) >= 4) >= 2) |
    #       (count_across_teeth(max_across_sites_at_one_tooth(DF_PD, MF_PD, DL_PD, ML_PD) >= 5) >= 1) )
    # -  comments: while the case definition requires "interproximal sites" to achieve a threshold, we still need to first consider at the tooth level because the definition requires the evidence coming from at least two teeth 
    # +  Moderate periodontitis: ≥2 interproximal sites with AL ≥4 mm (not on same tooth), or ≥2 interproximal sites with PD ≥5 mm (not on same tooth)
    # -  (count_across_teeth(max_across_sites_at_one_tooth(DF_CAL, MF_CAL, DL_CAL, ML_CAL) >= 4) >= 2) |
    #   (count_across_teeth(max_across_sites_at_one_tooth(DF_PD, MF_PD, DL_PD, ML_PD) >= 5) >= 2)
    # +  Severe periodontitis: ≥2 interproximal sites with AL ≥6 mm (not on same tooth) and ≥1 interproximal site with PD ≥5 mm
    # -  (count_across_teeth(max_across_sites_at_one_tooth(DF_CAL, MF_CAL, DL_CAL, ML_CAL) >=6 ) >= 2) &
    #   (count_across_teeth(max_across_sites_at_one_tooth(DF_PD, MF_PD, DL_PD, ML_PD) >= 5 ) >= 1)
    
    interdental <- c(Mesiobuccal, Distobuccal, Mesiolingual, Distolingual)
    
    #only the pc and la columns
    num_teeth <- names(d) %>% 
      str_subset(pattern = var_name) %>% 
      str_sub(min(tooth_num_digit), max(tooth_num_digit)) %>% 
      unique() %>% length()
    if(num_teeth != n_teeth_total)
      stop("Number of teeth in the dataset is not equal to n_teeth_total")
    num_site <- names(d) %>% 
      str_detect(pattern = var_name) %>%
      sum()
    if(num_site != 2*6*num_teeth)
      stop("Number of CAL and PD measurements in the dataset is not equal to 2*6*n_teeth_total")
    
    d_melt <- d %>% 
      select( all_of(id), matches(var_name) ) %>%
      pivot_longer(cols = matches(var_name),
                   names_to = "measurement_id",
                   values_to = "value",
                   values_drop_na = FALSE) %>%
      mutate(tooth= str_sub(measurement_id, min(tooth_num_digit), max(tooth_num_digit)),
             test = str_sub(measurement_id, min(measurement_digit), max(measurement_digit)),
             site = str_sub(measurement_id, min(site_digit), max(site_digit))) %>%
      mutate(is_na = is.na(value), 
             value = replace_na(value, -1),
             is_ident = site %in% interdental, is_B = (site == Buccal))
    d_melt_tooth <- d_melt %>% 
      group_by(pick(all_of(id)), tooth, test) %>% 
      summarise(n_NA_site = sum(is_na), 
                imax = max(value[is_ident]), # NA values will be -1. Avoid warnings
                B = value[is_B],
                .groups = "drop")
    # for 2018 case definition, almost all teeth across all subject in NHANES 
    # have detectable (>0) interdental CALs
    # d_melt_tooth %>% filter(test == "la") %>% group_by(tooth) %>%
    #   summarize(freq = sum(imax > 0)/sum(imax >= 0)) %>%
    #   arrange(desc(freq)) %>% head()
    # d %>% select(seqn, matches("ohx02la[dspa]"), matches("ohx31la[dspa]")) %>%
    #   slice_sample(n=20)
    
    d_indi <- d_melt_tooth %>%
      mutate(is_missing = (n_NA_site == 6)) %>%
      # n_NA_site can be 4-6 when imax == -1. 
      # In this case, imax should be 0 (not consider the tooth -> 0 due to missing) 
      mutate(test = case_match(test, CAL ~ "cal", PPD ~ "ppd")) %>% # for convenience
      pivot_wider(id_cols = all_of(c(id, "tooth")), 
                  values_from = all_of(c("n_NA_site", "imax", "B", "is_missing")),
                  names_from = test) %>%
      group_by(pick(all_of(id))) %>%
      summarise(n_NA_site_cal = sum(n_NA_site_cal), 
                n_NA_tooth_cal = sum(is_missing_cal),
                # 2018 case definition
                n_det_cal = sum(imax_cal>detectable_threshold),
                # v_det_cal = str_c(collapse = ",", 
                #                   tooth[imax_cal>detectable_threshold]),
                adj = if_else(
                  sum(imax_cal>detectable_threshold) == 2 &
                    abs(as.numeric(tooth[imax_cal>detectable_threshold])[1] - 
                          as.numeric(tooth[imax_cal>detectable_threshold])[2]) == 1,
                  "adj", "NA"),
                n_B_CALge3_PDge4 = sum(B_cal >=3 & B_ppd >= 4),
                imax_cal_max = max(imax_cal),
                # 2012 case definition
                n_i_CALge3 = sum(imax_cal >= 3),
                n_i_CALge4 = sum(imax_cal >= 4),
                n_i_CALge6 = sum(imax_cal >= 6),
                
                n_i_PPDge4 = sum(imax_ppd >= 4),
                n_i_PPDge5 = sum(imax_ppd >= 5),
                .groups = "drop"
      ) %>%
      mutate(ps_18 = ((n_det_cal >= 3)| 
                        (n_det_cal == 2 & adj != "adj")|
                        (n_B_CALge3_PDge4 >= 2)),
             mod_18 = (ps_18 & (imax_cal_max >= 3 & imax_cal_max <= 4)),
             sev_18 = (ps_18 & (imax_cal_max >= 5)),
             mil_12 = ((n_i_CALge3 >= 2) & (n_i_PPDge4 >= 2 | n_i_PPDge5 >= 1)),
             mod_12 = ((n_i_CALge4 >= 2) | (n_i_PPDge5 >= 2)),
             sev_12 = ((n_i_CALge6 >= 2) & (n_i_PPDge5 >= 1)),
             ps_12 = mil_12 | mod_12 | sev_12,
      ) %>%
      select(all_of(id),
             ps_12, mil_12, mod_12, sev_12, 
             ps_18, mod_18, sev_18,
             n_NA_site_cal, n_NA_tooth_cal) %>%
      mutate(n_NA_site_cal = n_NA_site_cal - n_NA_tooth_cal*6) # n_NA_site_cal
    
    d <- left_join(d, d_indi, by = id)
    
    return(d)
  }