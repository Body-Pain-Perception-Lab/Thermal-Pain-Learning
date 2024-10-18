# function to retrieve the data from OSF. Options for getting the data is either to get the preprocessed data i.e. rerun = FALSE
# alternatively to get the raw data and preprocess it choose rerun = TRUE. The returned list contains 3 elements, the preprocessed TPL dataframe, the TPL threshold dataframe and
# lastly a dataframe of participants that were removed from the analysis given the extreme amount of missing trials.
get_data <- function(rerun = FALSE) {
  if (rerun == FALSE) {
    # load the files from the Analysis folder
    return(list(
      data = read.csv(here::here("Analysis", "Cleaned-data", "TPL_df.csv")),
      thresholds = read.csv(here::here("Analysis", "Cleaned-data", "Thresholds.csv")),
      removers = read.csv(here::here("Analysis", "Cleaned-data", "removed.csv"))
    ))
  } else {
    
    # # This is the directory to where the raw behavioral data is stored both the experimental TPL and the first QST thresholding
    # 
    datapath <- here::here("matlab", "csv_files2", "Data")
    datafiles <- list.files(path = datapath, pattern = "*.tsv$", recursive = T)

    # initializing dataframes
    qst <- data.frame(NULL)
    data <- data.frame(NULL)
    # reading in the dataframes:
    for (i in 1:length(datafiles)) {
      if (str_detect(datafiles[i], "qst")) {
        qst1 <- read.delim(here(datapath, datafiles[i]), sep = "\t")
        qst <- rbind(qst, qst1)
      } else if (str_detect(datafiles[i], "tpl")) {
        painlearning1 <- read.delim(here(datapath, datafiles[i]), sep = "\t")
        data <- rbind(data, painlearning1)
      }
    }
    # converting factors to factors:
    columns <- c("cohort", "stim", "predResp", "predAcc", "block_type", "sequence", "desired_prob", "blockNumber", "cue", "location", "id")
    data[, columns] <- lapply(data[, columns], as.factor)
    

    # function to determine participants that should be excluded. Where the argument prob is the proportion of missed trials in %
    make_remover <- function(data, prob) {
      # count the number of times a participant missed a response:
      g <- data %>%
        filter(is.na(predResp)) %>%
        dplyr::group_by(id) %>%
        dplyr::summarize(n())
      
      # find all the times a participant had at least one missing values in their VAS-ratings
      q <- data %>% filter(is.na(vasResp_1) | is.na(vasResp_2) | is.na(vasResp_3))
      
      # now counting the times where all VAS-responses were missing
      q$sort <- ifelse(is.na(q$vasResp_1) & is.na(q$vasResp_2) & is.na(q$vasResp_3), 0, 1)
      q <- q %>% 
        filter(sort == 1) %>%
        dplyr::group_by(id) %>%
        dplyr::summarize(n())
      
      # now finding the participants that have missing trials above the threshold
      g <- data.frame(g) %>% filter(n.. > 306 * prob)
      q <- data.frame(q) %>% filter(n.. > 306 * prob)
      
      
      removers <- rbind(g, q)
    }
    # find the participants above 10% missing trials (either VAS or predictions):
    removers <- make_remover(data, 0.1)
    data <- data %>% filter(id %in% removers$id == F)
    thresh <- qst %>% filter(id %in% removers$id == F)
    
    thresh = thresh %>% mutate(age = ifelse(age == 0, 27,age), gender = ifelse(gender == 0, 2, gender))
    
    
    # renaming factors
    levels(data$cue) <- c("low-tone", "high-tone")
    levels(data$stim) <- c("cold", "warm", "TGI", "NaN")
    levels(data$predResp) <- c("cold", "warm", "NaN")
    levels(data$predAcc) <- c("Wrong", "Right", "TGI", "Wrong")
    
    
    # making a column that gives the expectedness of the stimuli, that is given the underlying probability and the cue was the outcome to be expected (E), unexpected(UE) or neutral(N).
    # all stimuli in 0.5 blocks are neutral and all TGI trials are coded as TGI. Here its important to note that the probabilities in desired_prob reflect the probability of a
    # high-tone to give a cold stimulus
    data$expected <- ifelse(data$desired_prob == "0.5" & data$stim != "TGI", a <- "N",
                            ifelse(data$desired_prob == "0.18" & data$cue == "low-tone" & data$stim == "cold", a <- "E",
                                   ifelse(data$desired_prob == "0.18" & data$cue == "high-tone" & data$stim == "warm", a <- "E",
                                          ifelse(data$desired_prob == "0.82" & data$cue == "low-tone" & data$stim == "warm", a <- "E",
                                                 ifelse(data$desired_prob == "0.82" & data$cue == "high-tone" & data$stim == "cold", a <- "E", a <- "UE")
                                          )
                                   )
                            )
    )
    
    data$expected <- ifelse(data$stim == "TGI", a <- "TGI", b <- data$expected)
    
    data$expected <- as.factor(data$expected)
    
    write.csv(removers %>% mutate(id = as.numeric(as.character(str_sub(id,-4,-1))),
                                  n.. = as.numeric(n..)),
              here::here("matlab", "removers2.csv"))
    
    return(list(data = data, thresholds = thresh, removers = removers))
  }
}


get_models <- function(osftoken){
  # get the preanalyzed data from osf:
  if (!dir.exists(here::here("Analysis","Models"))) {
    osfr::osf_auth(token = osf_token)

    TPL <- osfr::osf_retrieve_node("https://osf.io/q5z39/")

    TPL %>%
      osfr::osf_ls_files(pattern = "Models") %>%
      osfr::osf_download(path = here::here("Analysis"), recurse = TRUE, conflicts = "overwrite", progress = TRUE)
  }
  
}


# function to check simulated residudals of a fitted model:
residual_check <- function(model) {
  model_sim <- DHARMa::simulateResiduals(model)
  return(plot(model_sim))
}




# convinient function to make p-values to desired format
make_pvalue <- function(p_value) {
  if(is.na(p_value)){
    return("NA")
  }
  if (p_value > 0.05) {
    p_value <- round(p_value, 2)
    p <- paste("p = ", p_value)
  }
  if(p_value == 0.05){
    p = "p = .05 "
  }
  if (p_value < 0.05) {
    p <- "p < .05 "
  }
  if (p_value < 0.01) {
    p <- "p < .01"
  }
  if (p_value < 0.001) {
    p <- "p < .001"
  }
  if (p_value < 0.0001) {
    p <- "p < .0001"
  }
  
  return(p)
}

# get summary statistics from a generalized linear mixed effects model.
# Arguments: model is the model to get statistics on. Coefficients is the index of the coefficients when using summary(model) (that is excluding the intercept and until the coefficients arguemnt),
# round is the rounding of the statistics.

summary_stats_zoib_new <- function(model, coefficients, round, part, intercept = FALSE) {
  
  dfdataframe = parameters::model_parameters(model)
  if(part == "mu"){
    part = "conditional"
  }
  
  if("Component" %in% colnames(data.frame(dfdataframe))){
    dd = data.frame(dfdataframe) %>% filter(Component == part)
  }else{
    dd = data.frame(dfdataframe)
  }
  if(intercept){
    dd = dd[1:(coefficients+1),]
  }else{
    dd = dd[2:(coefficients+1),]
  }
  
  beta = round(dd$Coefficient,round)
  se = round(dd$SE, round)
  CI_low = round(dd$CI_low, round)
  CI_high = round(dd$CI_high, round)
  stat = round(dd$t, round)
  p = round(dd$p, round)
  df = round(dd$df_error, round)
  
  
  
  return(list(beta = beta, se = se, stat = stat,CI_low = CI_low,CI_high = CI_high,df = df, p = p))
  
}



# function to combine the dataframe of behavior aquired from MATLAB with the fitted HGF trajectories fitted in the TAPAS tool-box in MATLAB
get_hgf <- function(data, rerun) {
  
  
  if (rerun) {
    # run the main matlab script
    matlabr::run_matlab_script(here::here("matlab", "pain_main.m"))
  }
  
  
    # read in the trajectories
    trajfeel2 <- read.csv(here::here("matlab", "created files", "binary_response_feel_22.csv"), header = F)
    data$stim <- as.factor(data$stim)
    data$id <- as.factor(data$id)
    data$stim <- relevel(data$stim, ref = "cold")
    
    data1 <- data %>% filter(id %in% trajfeel2$V33)
    traj1 <- trajfeel2 %>% filter(V33 %in% data$id)
    traj1[, 34] <- rep(1:306, length(unique(traj1$V33)))
    
    # renaming all the columns from MATLAB
    traj1 <- traj1 %>% dplyr::rename(
      mu1 = V1, mu2 = V2,
      mu3 = V3, sa1 = V4,
      sa2 = V5, sa3 = V6,
      mu1hat = V7, mu2hat = V8,
      mu3hat = V9, sa1hat = V10,
      sa2hat = V11, sa3hat = V12,
      weight_v_1 = V13, weight_v_2 = V14,
      weight_v_3 = V15, weight_w_1 = V16,
      weight_w_2 = V17, volatility_pe1 = V18,
      volatility_pe2 = V19, volatility_pe3 = V20,
      update_ud_1 = V21, update_ud_2 = V22,
      update_ud_3 = V23, pweights_on_pe1 = V24,
      pweights_on_pe2 = V25, pweights_on_pe3 = V26,
      pwpe1 = V27, pwpe2 = V28, pwpe3 = V29,
      learn_rate1 = V30, learn_rate2 = V31,
      learn_rate3 = V32, id = V33,
      trial = V34
    )
    
    traj1$id <- as.factor(traj1$id)
    df <- inner_join(data1, traj1, by = c("id", "trial"))
    
    df$burnbeta <- df$vasResp_3 / 100
    df$coldbeta <- df$vasResp_1 / 100
    df$warmbeta <- df$vasResp_2 / 100
    
    return(df)
   
  }
    


parameter_model_recovery = function(){

  # run the main matlab script
  matlabr::run_matlab_script(here::here("matlab", "model_parameter_recovery.m"))
  
}


make_new_reporting = function(stats, number, Z, inc_df = F){
  
  if(Z){
    z_stat = stats$stat[number]
    p_value = stats$p[number]
    ci = list(low = as.numeric(stats$beta[number])-2*as.numeric(stats$std[number]), high = as.numeric(stats$beta[number])+2*as.numeric(stats$std[number]))
    beta_stat = stats$beta[number]
    if(inc_df){
      text = paste0("$\\", "beta", "$"," = ",beta_stat,", 95% CI = [", ci$low, " ; ", ci$high,"]",", Z = ", z_stat, ", ", make_pvalue(p_value))
    }
    
    return(text)
    
  }else if(!Z){
    beta_stat = stats$beta[number]
    t_stat = stats$stat[number]
    p_value = stats$p[number]
    ci = list(low = stats$CI_low[number], high = stats$CI_high[number])
    df = stats$df[number]
    if(inc_df){
      text = paste0("$\\", "beta", "$"," = ",beta_stat, ", 95% CI = [", ci$low, "; ", ci$high,"]",", t(",df,") = ",t_stat,", ",make_pvalue(p_value))
    }else{
      text = paste0("$\\", "beta", "$"," = ",beta_stat, ", 95% CI = [", ci$low, "; ", ci$high,"]",", ",make_pvalue(p_value))
      
    }
    
    return(text)
    
  }
  
  return("error")
}




get_all_stats_and_plots = function(){
  
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load("renv", "here", "knitr")
  required_packages = c("rmarkdown","tidyverse","lmerTest","brms","ggeffects", "gamlss", "DHARMa",
                        "glmmTMB", "cowplot", "scales", "matlabr","osfr","parameters")
  lapply(required_packages, library, character.only = TRUE)
  detach("package:renv", unload = T)
  
  source(here::here("scripts", "utils.R"))
  

  # extract data to rerun the analysis
  data <- get_data(rerun = F)
  df <- data$data
  trajfeel2 <- get_hgf(df, rerun = F)
  
  if(!dir.exists(here::here("Analysis","Models"))){
    # Retrieve the OSF authentication token if the file exists
    osf_file_path <- here::here("osf","osf.txt")
    if (file.exists(osf_file_path)) {
      osf_token <- read_lines(osf_file_path)[1]
    } else {
      stop("OSF token file does not exist!")
    }
    
    
  }
  
  models = list.files(here::here("Analysis","Models"), full.names = T)
  
  lapply(models, function(model) {
    load(model, envir = .GlobalEnv)  # Loads each file into the global environment
  })
  
  
  ###################################
  # starting with stuff for the manuscript on the methods:
  nf <- data$thresholds %>%
    group_by(gender, id) %>%
    summarize(n = n()) %>%
    filter(gender == 1) %>%
    mutate(females = length(gender))
  
  temps <- data$data %>%
    group_by(id) %>%
    summarize(maxw = max(targetT_1, na.rm = T), maxc = min(targetT_2, na.rm = T)) %>%
    ungroup() %>%
    summarize(meanw = mean(maxw), sdw = sd(maxw), meanc = mean(maxc), sdc = sd(maxc))
  
  nr <- length(unique(data[[3]]$id))
  np <- length(unique(data[[1]]$id))
  
  
  age <- data$data %>%
    filter(age != 0) %>%
    summarize(maxage = max(age, na.rm = T), minage = min(age, na.rm = T), meanage = mean(age, na.rm = T),sdage = sd(age, na.rm = T))
  
  
  ex_ids = read.csv(here::here("matlab","excluded_ids.csv"))
  
  qq = data$thresholds %>% filter(!id %in% ex_ids$excluded_ids)
  
  MRI_meanage_sdage = round(qq  %>% summarize(meanage = mean(age), sdage = sd(age)),1)
  
  MRI_nf = qq %>%
    group_by(gender, id) %>%
    summarize(n = n()) %>%
    filter(gender == 1) %>%
    mutate(females = length(gender))
  

  ################################################
  #starting with how the percived quaility of the TGI shapes whether you are correct:
  
  
  df$predAcc2 <- lead(df$predAcc, 1)
  
  df_percieved_TGI_acc <- df %>%
    group_by(stim, id, trial) %>%
    summarize(coolness = vasResp_1 / (vasResp_1 + vasResp_2), warmness = vasResp_2 / (vasResp_1 + vasResp_2)) %>%
    inner_join(df) %>%
    mutate(stim = as.factor(stim), id = as.factor(id), desired_prob = as.factor(desired_prob)) %>%
    dplyr::select(predAcc2, coolness, cue, desired_prob, trial, id) %>%
    filter(stim == "TGI", desired_prob != "0.5") %>%
    mutate(predAcc2 = ifelse(predAcc2 == "Right", 1, 0), cue = as.factor(cue)) %>%
    drop_na() %>%
    mutate(contingency = as.factor(ifelse(desired_prob == "0.82" & cue == "high-tone" | desired_prob == "0.18" & cue == "low-tone", "Predicts cold", "Predicts warm")))
  
  
  df_percieved_TGI_acc$contingency <- relevel(df_percieved_TGI_acc$contingency, ref = "Predicts warm")
  
  model_percieved_TGI_acc_contingency_warm <- gamlss::gamlss(predAcc2 ~ coolness * contingency + trial + random(id),
                                                             family = BI(mu.link = "logit"),
                                                             control = gamlss.control(n.cyc = 100, trace = T),
                                                             data = df_percieved_TGI_acc
  )
  
  # extrating estimates to give probabilities in the manuscript
  
  int_percieved_TGI_acc_contingency_cold <- summary_stats_zoib_new(model_percieved_TGI_acc_contingency_cold, coefficients = 1, round = 2, part = "mu", intercept = TRUE)
  
  stats_percieved_TGI_acc_contingency_cold <- summary_stats_zoib_new(model_percieved_TGI_acc_contingency_cold, coefficients = 4, round = 2, part = "mu")
  
  int_percieved_TGI_acc_contingency_warm <- summary_stats_zoib_new(model_percieved_TGI_acc_contingency_warm, coefficients = 1, round = 2, part = "mu", intercept = TRUE)
  
  stats_percieved_TGI_acc_contingency_warm <- summary_stats_zoib_new(model_percieved_TGI_acc_contingency_warm, coefficients = 4, round = 2, part = "mu")
  
  base::save(list = "stats_percieved_TGI_acc_contingency_warm", file = here::here("Analysis","Temp", "reporting_statistics_stats_percieved_TGI_acc_contingency_warm.RData"))
  
  
  # what does the effects mean?
  #calculating what this means in how much the probability of answering correct changes going from having the "right" percept of the TGI (given the contingency) compared to precieving the opposite.
  
  #because the effects are opposite we take the absolute value and calculate the average effect of a one unit increase in coolness (i.e. going from responding completely cold to completely warm)
  mean_percieved_TGI_acc_contingency <- mean(abs(c(stats_percieved_TGI_acc_contingency_warm$beta[1], stats_percieved_TGI_acc_contingency_cold$beta[1])))
  
  #we then get the uncertainty associated with these values:
  mean_percieved_TGI_acc_contingency_plus_2sd <- mean(abs(c(
    stats_percieved_TGI_acc_contingency_warm$beta[1] - 2 * stats_percieved_TGI_acc_contingency_warm$se[1],
    stats_percieved_TGI_acc_contingency_cold$beta[1] + 2 * stats_percieved_TGI_acc_contingency_cold$se[1]
  )))
  
  mean_percieved_TGI_acc_contingency_minus_2sd <- mean(abs(c(
    stats_percieved_TGI_acc_contingency_warm$beta[1] + 2 * stats_percieved_TGI_acc_contingency_warm$se[1],
    stats_percieved_TGI_acc_contingency_cold$beta[1] - 2 * stats_percieved_TGI_acc_contingency_cold$se[1]
  )))
  
  
  #average intercept
  mean_acc_TGI_acc <- mean(abs(c(int_percieved_TGI_acc_contingency_cold$beta[1], int_percieved_TGI_acc_contingency_warm$beta[1])))
  
  #average effect of coolness
  effect_TGI_acc <- brms::inv_logit_scaled(mean_acc_TGI_acc + mean_percieved_TGI_acc_contingency)- brms::inv_logit_scaled(mean_acc_TGI_acc)
  
  #and 95% confidence intervals
  effect_TGI_acc_95ci <- c(
    brms::inv_logit_scaled(mean_acc_TGI_acc + mean_percieved_TGI_acc_contingency_plus_2sd) - brms::inv_logit_scaled(mean_acc_TGI_acc),
    brms::inv_logit_scaled(mean_acc_TGI_acc + mean_percieved_TGI_acc_contingency_minus_2sd) - brms::inv_logit_scaled(mean_acc_TGI_acc)
  )
  
  
  save(
    list = c(
      "effect_TGI_acc_95ci", "effect_TGI_acc"
    ),
    file = here::here("Analysis","Temp", "effect_TGI_acc.RData")
  )
  
  
  ## for plots
  percieved_TGI_acc_contingency_cold_95 <- ggeffects::ggpredict(model_percieved_TGI_acc_contingency_warm,
                                                                terms = c("coolness[all]", "contingency"),
                                                                type = "random")
  
  percieved_TGI_acc_contingency_cold_80 <- ggeffects::ggpredict(model_percieved_TGI_acc_contingency_warm,
                                                                terms = c("coolness[all]", "contingency"),
                                                                type = "random",
                                                                ci_level = 0.8)
  
  percieved_TGI_acc_contingency_cold_50 <- ggeffects::ggpredict(model_percieved_TGI_acc_contingency_warm,
                                                                terms = c("coolness[all]", "contingency"),
                                                                type = "random",
                                                                ci_level = 0.5)
  
  
  # save extra model results for plotting
  
  
  save(
    list = c(
      "percieved_TGI_acc_contingency_cold_95", "percieved_TGI_acc_contingency_cold_80", "percieved_TGI_acc_contingency_cold_50"
    ),
    file = here::here("Analysis", "Workspace", "percieved_TGI.RData")
  )
  
  
  ################################################
  # Now the three way interaction of binary responses:
  
  ## wrangling the data
  df_expectation <- df %>%
    mutate(coldbeta = vasResp_1 / 100, warmbeta = vasResp_2 / 100) %>%
    dplyr::select(predResp, coldbeta, vasRT_1, vasRT_2, trial, id, stim, warmbeta) %>%
    drop_na() %>%
    filter(predResp != "NaN") %>%
    mutate(stim = as.factor(stim), id = as.factor(id), trial = scale(trial)[,1]) %>%
    pivot_longer(cols = c("coldbeta", "warmbeta"), names_to = "Ratingscale") %>%
    pivot_longer(cols = c("vasRT_1", "vasRT_2"), names_to = "Ratingscale_rt", values_to = "RT") %>%
    mutate(RateCon = as.factor(ifelse(stim == "cold" &
                                        Ratingscale == "coldbeta" | stim == "warm" &
                                        Ratingscale == "warmbeta", "Con", "Incon"))) %>%
    filter(stim != "TGI") %>%
    mutate(predResp = as.factor(predResp)) %>%
    filter(Ratingscale == "coldbeta" & Ratingscale_rt == "vasRT_1" | Ratingscale == "warmbeta" & Ratingscale_rt == "vasRT_2")
  
  
  df_expectation$predResp <- relevel(df_expectation$predResp, ref = "cold")
  df_expectation$stim <- relevel(df_expectation$stim, ref = "cold")
  df_expectation$RateCon <- relevel(df_expectation$RateCon, ref = "Con")
  
  inncous_expectation_cold_con <- update(model_INNOUOUS_expectation_warm_incon, data = df_expectation)
  
  stat_INNOUOUS_expectation_cold_con_mu <- summary_stats_zoib_new(inncous_expectation_cold_con, coefficients = 8, round = 2, part = "mu")
  
  base::save(list = "stat_INNOUOUS_expectation_cold_con_mu", file = here::here("Analysis","Temp", "reporting_statistics_stat_INNOUOUS_expectation_cold_con_mu.RData"))
  
  ## now with the other factors as references
  df_expectation$predResp <- relevel(df_expectation$predResp, ref = "cold")
  df_expectation$stim <- relevel(df_expectation$stim, ref = "cold")
  df_expectation$RateCon <- relevel(df_expectation$RateCon, ref = "Incon")
  
  inncous_expectation_cold_incon <- update(model_INNOUOUS_expectation_warm_incon, data = df_expectation)
  
  stat_INNOUOUS_expectation_cold_incon_mu <- summary_stats_zoib_new(inncous_expectation_cold_incon, coefficients = 8, round = 2, part = "mu")
  
  
  base::save(list = "stat_INNOUOUS_expectation_cold_incon_mu", file = here::here("Analysis","Temp", "reporting_statistics_stat_INNOUOUS_expectation_cold_incon_mu.RData"))
  
  ## now with the other factors as references
  df_expectation$predResp <- relevel(df_expectation$predResp, ref = "warm")
  df_expectation$stim <- relevel(df_expectation$stim, ref = "warm")
  df_expectation$RateCon <- relevel(df_expectation$RateCon, ref = "Con")
  
  model_INNOUOUS_expectation_warm_con <- update(model_INNOUOUS_expectation_warm_incon, data = df_expectation)
  
  stat_INNOUOUS_expectation_warm_con_mu <- summary_stats_zoib_new(model_INNOUOUS_expectation_warm_con, coefficients = 8, round = 2, part = "mu")
  
  
  base::save(list = "stat_INNOUOUS_expectation_warm_con_mu", file = here::here("Analysis","Temp", "reporting_statistics_stat_INNOUOUS_expectation_warm_con_mu.RData"))
  
  
  ################################################
  # Then the three way interaction of the continuous belief to cold:
  
  
  hgf_expectation_tgi <- trajfeel2 %>%
    filter(stim != "NaN") %>%
    mutate(coldbeta = vasResp_1 / 100, warmbeta = vasResp_2 / 100, burnbeta = vasResp_3 / 100) %>%
    dplyr::select(stim, id, sa1hat, trial, coldbeta, warmbeta, burnbeta, cue, mu1hat) %>%
    drop_na() %>%
    mutate(stim = as.factor(stim), id = as.factor(id)) %>%
    pivot_longer(cols = c("coldbeta", "warmbeta", "burnbeta"), names_to = "Ratingscale") %>%
    mutate(belief_to_cold = ifelse(cue == "high-tone", mu1hat, 1 - mu1hat))%>% 
    mutate(trial = scale(trial)[,1], Ratingscale = as.factor(Ratingscale))
  
  
  # Full model with burning ratings
  hgf_expectation_tgi$Ratingscale <- relevel(hgf_expectation_tgi$Ratingscale, ref = "warmbeta")
  hgf_expectation_tgi$stim <- relevel(hgf_expectation_tgi$stim, ref = "warm")
  
  model_figure4_warm_warm <-  gamlss(value ~ belief_to_cold * Ratingscale * stim + trial + random(id),
                                     nu.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                     tau.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                     sigma.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                     data = hgf_expectation_tgi,
                                     family = BEINF(mu.link = "logit", sigma.link = "logit", nu.link = "logit", tau.link = "logit"),
                                     control = gamlss.control(n.cyc = 100, trace = T)
  )
  

  stat_expectation_HGF_warm_warm <- summary_stats_zoib_new(model_figure4_warm_warm, coefficients = 18, round = 2, part = "mu")
  
  base::save(list = "stat_expectation_HGF_warm_warm", file = here::here("Analysis","Temp", "reporting_statistics_stat_expectation_HGF_warm_warm.RData"))
  
  
  #marginal effects for the manuscript
  
  figure4a_95  <- ggeffects::ggpredict(model_figure4_warm_warm, terms = c("belief_to_cold[all]","stim","Ratingscale"), type = "random")
  
  figure4a_80  <- ggeffects::ggpredict(model_figure4_warm_warm, terms = c("belief_to_cold[all]","stim","Ratingscale"), type = "random",ci.lvl = 0.8)
  
  figure4a_50  <- ggeffects::ggpredict(model_figure4_warm_warm, terms = c("belief_to_cold[all]","stim","Ratingscale"), type = "random",ci.lvl = 0.5)
  
  
  # save extra model results for plotting
  #
  save(
    list = c(
      "figure4a_95", "figure4a_80", "figure4a_50"
    ),
    file = here::here("Analysis", "Workspace", "figure4a.RData")
  )
  

  
  hgf_expectation_tgi$Ratingscale <- relevel(hgf_expectation_tgi$Ratingscale, ref = "coldbeta")
  hgf_expectation_tgi$stim <- relevel(hgf_expectation_tgi$stim, ref = "TGI")
  
  model_figure4_full_TGI_cold <-  gamlss(value ~ belief_to_cold * Ratingscale * stim + trial + random(id),
                                         nu.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                         tau.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                         sigma.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                         data = hgf_expectation_tgi,
                                         family = BEINF(mu.link = "logit", sigma.link = "logit", nu.link = "logit", tau.link = "logit"),
                                         control = gamlss.control(n.cyc = 100, trace = T)
  )
  
  stat_expectation_HGF_TGI_cold <- summary_stats_zoib_new(model_figure4_full_TGI_cold, coefficients = 18, round = 2, part = "mu")
  
  base::save(list = "stat_expectation_HGF_TGI_cold", file = here::here("Analysis","Temp", "reporting_statistics_stat_expectation_HGF_TGI_cold.RData"))
  
  
  
  hgf_expectation_tgi$Ratingscale <- relevel(hgf_expectation_tgi$Ratingscale, ref = "warmbeta")
  hgf_expectation_tgi$stim <- relevel(hgf_expectation_tgi$stim, ref = "TGI")
  
  model_figure4_full_TGI_warm <-  gamlss(value ~ belief_to_cold * Ratingscale * stim + trial + random(id),
                                         nu.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                         tau.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                         sigma.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                         data = hgf_expectation_tgi,
                                         family = BEINF(mu.link = "logit", sigma.link = "logit", nu.link = "logit", tau.link = "logit"),
                                         control = gamlss.control(n.cyc = 100, trace = T)
  )
  
  
  stat_expectation_HGF_TGI_warm <- summary_stats_zoib_new(model_figure4_full_TGI_warm, coefficients = 18, round = 2, part = "mu")
  
  base::save(list = "stat_expectation_HGF_TGI_warm", file = here::here("Analysis","Temp", "reporting_statistics_stat_expectation_HGF_TGI_warm.RData"))
  
  
  
  hgf_expectation_tgi$Ratingscale <- relevel(hgf_expectation_tgi$Ratingscale, ref = "burnbeta")
  hgf_expectation_tgi$stim <- relevel(hgf_expectation_tgi$stim, ref = "TGI")
  
  model_figure4_full_TGI_burn <-  gamlss(value ~ belief_to_cold * Ratingscale * stim + trial + random(id),
                                         nu.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                         tau.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                         sigma.formula = ~ belief_to_cold * Ratingscale * stim+ trial + random(id),
                                         data = hgf_expectation_tgi,
                                         family = BEINF(mu.link = "logit", sigma.link = "logit", nu.link = "logit", tau.link = "logit"),
                                         control = gamlss.control(n.cyc = 100, trace = T)
  )
  
  
  stat_expectation_HGF_TGI_burn <- summary_stats_zoib_new(model_figure4_full_TGI_burn, coefficients = 18, round = 2, part = "mu")
  
  base::save(list = "stat_expectation_HGF_TGI_burn", file = here::here("Analysis","Temp", "reporting_statistics_stat_expectation_HGF_TGI_burn.RData"))
  
  
  reporting = list.files(here::here("Analysis","Temp"), full.names = T)
  for(file in reporting){
    load(file)
    
  }
  
  
  # Get a list of all variable names that start with "stats_" and "stat_"
  stat_vars <- ls(pattern = "^stat_")
  stats_vars <- ls(pattern = "^stats_")
  effect_vars <- ls(pattern = "^effect_")
  # # Combine the variable names into a single vector
  vars_to_save <- c(effect_vars, stat_vars, stats_vars, "age", "nf", "nr", "np", "temps", "MRI_nf","MRI_meanage_sdage")
  
  # # Save the variables to a file
  
  
  
  base::save(list = vars_to_save, file = here::here("Manuscripts", "Workspace", "reporting_statistics.RData"))

    
}
