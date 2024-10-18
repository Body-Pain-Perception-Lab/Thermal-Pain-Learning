#font size and style for all plots for everything
font = "sans"
font_size = 12
font_size_small = 8
axis_width = 1

text = ggplot2::theme(text = ggplot2::element_text(family = font, size = font_size))
theme = cowplot::theme_cowplot()

#patch theme
patchtheme = ggplot2::theme(plot.tag = ggplot2::element_text(size = font_size+4,       
                                                             family = "sans",     
                                                             face = "bold",            
                                                             hjust = 0.5,              
                                                             vjust = 0.5))


# plotting the parameter recovery that was run in MATLAB
parameter_recovery <- function(model) {
  
  
  make_pr_plot = function(data, correlation, name,xlim, ylim,xlabel, ylabel){
    
    plot <- data %>% ggplot(aes(x = data[,paste0(name,"_simulated")], y = data[,paste0(name,"_recovered")])) +
      geom_point() +
      theme_classic() +
      ggtext::geom_richtext(aes(x = xlabel, y = ylabel, label = paste("r = ", round(correlation$estimate, 2), " ", "[", round(correlation$conf.int[1], 2), ";", round(correlation$conf.int[2], 2), "] <br>", make_pvalue(correlation$p.value))), size = 3,stat = "unique") +
      ylab(bquote(Recovered ~ .(as.name(name)))) +
      xlab(bquote(Simulated ~ .(as.name(name))))+
      coord_cartesian(xlim = xlim, ylim = ylim)+theme+text+
      theme(legend.position = "center",
            legend.box = "horizontal",
            legend.direction="horizontal",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            legend.spacing.x = unit(0.35, "cm"),
            plot.title = element_text(hjust = 0.5),
            legend.title = element_blank(),
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width))
    
    return(plot)
  }
  
  if (model == "hgf_2level") {
    # read the data
    #pr_hgf <- data.frame(t(read.csv(here::here("Analysis", "Computational_modeling", "Parameter_recovery", "hgf2_parameter_recovery.csv"), header = F)))
    pr_hgf <- data.frame(t(read.csv(here::here("matlab", "created files", "hgf2_parameter_recovery.csv"), header = F)))
    
    # rename columns
    pr_hgf <- pr_hgf %>% dplyr::rename(
      omega_simulated = X1, omega_recovered = X2,
      zeta_simulated = X3, zeta_recovered = X4
    )
    
    # making correlation analysis of the simulated and recovered parameters
    cortest_beta <- cor.test(pr_hgf$zeta_simulated, pr_hgf$zeta_recovered)
    cortest_omega <- cor.test(pr_hgf$omega_simulated, pr_hgf$omega_recovered)
    
    
    zeta = make_pr_plot(data = pr_hgf,correlation = cortest_beta, name = "zeta",xlim = c(0,15), ylim = c(0,30),xlabel = 10, ylabel = 30)
    omega2 = make_pr_plot(data = pr_hgf,correlation = cortest_omega, name = "omega",xlim = c(-12,2), ylim = c(-12,1),xlabel = -8, ylabel = -1.5)
    
    
    
    # plotting it in the same plot using patchwork
    plot <- zeta + omega2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
    
    plot
    return(plot)
  }
  
  # first the HGF model
  if (model == "hgf") {
    # read the data
    #pr_hgf <- data.frame(t(read.csv(here::here("Analysis", "Computational_modeling", "Parameter_recovery", "hgf_parameter_recovery.csv"), header = F)))
    pr_hgf <- data.frame(t(read.csv(here::here("matlab", "created files", "hgf_parameter_recovery.csv"), header = F)))
    
    # rename columns
    pr_hgf <- pr_hgf %>% dplyr::rename(
      omega_simulated = X1, omega_recovered = X2,
      theta_simulated = X3, theta_recovered = X4,
      kappa_simulated = X5, kappa_recovered = X6,
      zeta_simulated = X7, zeta_recovered = X8
    )
    
    # making correlation analysis of the simulated and recovered parameters
    cortest_beta <- cor.test(pr_hgf$zeta_simulated, pr_hgf$zeta_recovered)
    cortest_omega <- cor.test(pr_hgf$omega_simulated, pr_hgf$omega_recovered)
    
    
    cortest_omega31 <- cor.test(pr_hgf$theta_simulated, pr_hgf$theta_recovered)
    cortest_kappa <- cor.test(pr_hgf$kappa_simulated, pr_hgf$kappa_recovered)
    
    
    zeta = make_pr_plot(data = pr_hgf,correlation = cortest_beta, name = "zeta",xlim = c(0,15), ylim = c(0,30),xlabel = 10, ylabel = 25)
    omega2 = make_pr_plot(data = pr_hgf,correlation = cortest_omega, name = "omega",xlim = c(-12,2), ylim = c(-12,2),xlabel = -8, ylabel = -1.5)
    theta = make_pr_plot(data = pr_hgf,correlation = cortest_omega31, name = "theta",xlim = c(-12,4), ylim = c(-12,4),xlabel = -8, ylabel = -1)
    kappa = make_pr_plot(data = pr_hgf,correlation = cortest_kappa, name = "kappa",xlim = c(0,2), ylim = c(0,3),xlabel = 1.5, ylabel = 2.2)
    
    
    
    # plotting it in the same plot using patchwork
    plot <- zeta + omega2 + theta + kappa + plot_layout(guides = "collect") & theme(legend.position = "bottom")
    
    
    return(plot)
  }
  # for the rescorla wagner learning model:
  if (model == "rw") {
    # reading the data
    #pr_rw <- data.frame(t(read.csv(here::here("Analysis", "Computational_modeling", "Parameter_recovery", "rw_parameter_recovery.csv"), header = F)))
    pr_rw <- data.frame(t(read.csv(here::here("matlab", "created files", "rw_parameter_recovery.csv"), header = F)))
    
    
    #renaming parameters
    pr_rw <- pr_rw %>% dplyr::rename(alpha_simulated = X1,
                                     alpha_recovered = X2,
                                     zeta_simulated = X3,
                                     zeta_recovered = X4)
    
    #making correlation test
    cortest_zeta <- cor.test(pr_rw$zeta_simulated, pr_rw$zeta_recovered)
    cortest_alpha <- cor.test(pr_rw$alpha_simulated, pr_rw$alpha_recovered)
    
    
    
    zeta = make_pr_plot(data = pr_rw,correlation = cortest_zeta, name = "zeta",xlim = c(0,15), ylim = c(0,30),xlabel = 5, ylabel = 25)
    alpha = make_pr_plot(data = pr_rw,correlation = cortest_alpha, name = "alpha",xlim = c(0,1), ylim = c(0,1),xlabel = 0.40, ylabel = 0.75)
    
    #combining plots with patchwork
    plot <- zeta + alpha + plot_layout(guides = "collect") & theme(legend.position = "bottom")
    return(plot)
  }
  # for the sutton k1 learning model:
  if (model == "su1") {
    #read the data
    #pr_su1 <- data.frame(t(read.csv(here::here("Analysis", "Computational_modeling", "Parameter_recovery", "su1_parameter_recovery.csv"), header = F)))
    pr_su1 <- data.frame(t(read.csv(here::here("matlab", "created files", "su1_parameter_recovery.csv"), header = F)))
    
    #rename the variables
    pr_su1 <- pr_su1 %>% dplyr::rename(zeta_simulated = X3,
                                       zeta_recovered = X4,
                                       mu_simulated = X1,
                                       mu_recovered = X2)
    
    #correlation test
    cortest_zeta <- cor.test(pr_su1$zeta_simulated, pr_su1$zeta_recovered)
    cortest_mu <- cor.test(pr_su1$mu_simulated, pr_su1$mu_recovered)
    
    zeta = make_pr_plot(data = pr_su1,correlation = cortest_zeta, name = "zeta",xlim = c(0,15), ylim = c(0,30),xlabel = 7, ylabel = 23)
    mu = make_pr_plot(data = pr_su1,correlation = cortest_mu, name = "mu",xlim = c(0,10), ylim = c(0,20),xlabel = 3, ylabel = 10)
    
    #combining the plots
    plot <- zeta + mu + plot_layout(guides = "collect") & theme(legend.position = "bottom")
    
    return(plot)
  }
  if (model == "ph") {
    # reading the data
    #pr_rw <- data.frame(t(read.csv(here::here("Analysis", "Computational_modeling", "Parameter_recovery", "rw_parameter_recovery.csv"), header = F)))
    pr_ph <- data.frame(t(read.csv(here::here("matlab", "created files", "ph_parameter_recovery.csv"), header = F)))
    
    
    #renaming parameters
    pr_ph <- pr_ph %>% dplyr::rename(S_simulated = X1,
                                     S_recovered = X2,
                                     zeta_simulated = X3,
                                     zeta_recovered = X4)
    
    #making correlation test
    cortest_zeta <- cor.test(pr_ph$zeta_simulated, pr_ph$zeta_recovered)
    cortest_S <- cor.test(pr_ph$S_simulated, pr_ph$S_recovered)
    
    
    
    zeta = make_pr_plot(data = pr_ph,correlation = cortest_zeta, name = "zeta",xlim = c(0,15), ylim = c(0,30),xlabel = 5, ylabel = 25)
    S = make_pr_plot(data = pr_ph,correlation = cortest_S, name = "S",xlim = c(0,1), ylim = c(0,1),xlabel = 0.40, ylabel = 0.75)
    
    #combining plots with patchwork
    plot <- zeta + S + plot_layout(guides = "collect") & theme(legend.position = "bottom")
    return(plot)
  }
  if (model == "ph1") {
    # reading the data
    #pr_rw <- data.frame(t(read.csv(here::here("Analysis", "Computational_modeling", "Parameter_recovery", "rw_parameter_recovery.csv"), header = F)))
    pr_ph <- data.frame(t(read.csv(here::here("matlab", "created files", "modph_parameter_recovery.csv"), header = F)))
    
    
    #renaming parameters
    pr_ph <- pr_ph %>% dplyr::rename(kappa_simulated = X1,
                                     kappa_recovered = X2,
                                     eta_simulated = X3,
                                     eta_recovered = X4,
                                     zeta_simulated = X5,
                                     zeta_recovered = X6)
    
    #making correlation test
    cortest_zeta <- cor.test(pr_ph$zeta_simulated, pr_ph$zeta_recovered)
    cortest_kappa <- cor.test(pr_ph$kappa_simulated, pr_ph$kappa_recovered)
    cortest_eta <- cor.test(pr_ph$eta_simulated, pr_ph$eta_recovered)
    
    
    
    zeta = make_pr_plot(data = pr_ph,correlation = cortest_zeta, name = "zeta",xlim = c(0,15), ylim = c(0,30),xlabel = 5, ylabel = 25)
    
    kappa = make_pr_plot(data = pr_ph,correlation = cortest_kappa, name = "kappa",xlim = c(0,1), ylim = c(0,1),xlabel = 0.40, ylabel = 0.75)
    eta = make_pr_plot(data = pr_ph,correlation = cortest_eta, name = "eta",xlim = c(0,1), ylim = c(0,1),xlabel = 0.40, ylabel = 0.75)
    
    
    #combining plots with patchwork
    plot <- (zeta + eta)/kappa + plot_layout(guides = "collect") & theme(legend.position = "bottom")
    return(plot)
  }
  
  
  
}


#function to analyse the model_recovery done in the MATLAB script.
get_model_recovery <- function() {
  
  #reading the model_recovery data (this is log likelihoods of models where the data simulating model was the HGF)
  mr_hgf <- data.frame(read.csv(here::here("matlab", "created files", "hgf_data_model_recovery.csv"), header = F))[-1, ]
  
  #renaming columns so they match the model fitted
  mr_hgf <- mr_hgf %>% dplyr::rename(rw = V1,
                                     hgf = V2,
                                     su1 = V3,
                                     ph = V4)
  
  #find the highest log likelihood per simulation
  mr_hgf$lc <- colnames(mr_hgf)[apply(mr_hgf, 1, which.max)]
  
  #count the number of times a model "won" i.e had the highest log likelihood.
  hgf <- mr_hgf %>%
    group_by(lc) %>%
    summarize(hgf = n())
  
  #we can even plot the parameter space that we simulated and see which model "won" where
  zeta <- data.frame(zeta = t(read.csv(here::here("matlab", "created files", "simulated_zeta.csv"), header = F)))
  omega <- data.frame(omega = t(read.csv(here::here("matlab", "created files", "simulated_omega2.csv"), header = F)))
  
  mr_hgf %>%
    mutate(zeta = zeta$zeta, omega = omega$omega) %>%
    ggplot(aes(x = omega, y = zeta, col = lc)) +
    geom_point() +
    theme_classic()
  
  
  #same for when the data was simulated from the rescorla wagner:
  mr_rw <- data.frame(read.csv(here::here("matlab", "created files", "rw_data_model_recovery.csv"), header = F))[-1, ]
  mr_rw <- mr_rw %>% dplyr::rename(rw = V1, hgf = V2, su1 = V3, ph = V4)
  mr_rw$lc <- colnames(mr_rw)[apply(mr_rw, 1, which.max)]
  rw <- mr_rw %>%
    group_by(lc) %>%
    summarize(rw = n())
  
  #and can plot this as well
  zeta <- data.frame(zeta = t(read.csv(here::here("matlab", "created files", "simulated_zeta.csv"), header = F)))
  alpha <- data.frame(alpha = t(read.csv(here::here("matlab", "created files", "simulated_alpha.csv"), header = F)))
  
  mr_rw %>%
    mutate(zeta = zeta$zeta, alpha = alpha$alpha) %>%
    ggplot(aes(x = alpha, y = zeta, col = lc)) +
    geom_point() +
    theme_classic()
  
  
  
  #same when the data simulating model was the sutton K1
  mr_su1 <- data.frame(read.csv(here::here("matlab", "created files", "su1_data_model_recovery.csv"), header = F))[-1, ]
  mr_su1 <- mr_su1 %>% dplyr::rename(rw = V1, hgf = V2, su1 = V3, ph = V4)
  mr_su1$lc <- colnames(mr_su1)[apply(mr_su1, 1, which.max)]
  su1 <- mr_su1 %>%
    group_by(lc) %>%
    summarize(su1 = n())
  
  #and plotting
  
  zeta <- data.frame(zeta = t(read.csv(here::here("matlab", "created files", "simulated_zeta.csv"), header = F)))
  mu <- data.frame(mu = t(read.csv(here::here("matlab", "created files", "simulated_mu.csv"), header = F)))
  
  mr_su1 %>%
    mutate(zeta = zeta$zeta, mu = mu$mu) %>%
    ggplot(aes(x = mu, y = zeta, col = lc)) +
    geom_point() +
    theme_classic()
  
  
  #same when the data simulating model was the pearce hall
  mr_ph <- data.frame(read.csv(here::here("matlab", "created files", "ph_data_model_recovery.csv"), header = F))[-1, ]
  mr_ph <- mr_ph %>% dplyr::rename(rw = V1, hgf = V2, su1 = V3, ph = V4)
  mr_ph$lc <- colnames(mr_ph)[apply(mr_ph, 1, which.max)]
  ph <- mr_ph %>%
    group_by(lc) %>%
    summarize(ph = n())
  
  #and plotting
  
  zeta <- data.frame(zeta = t(read.csv(here::here("matlab", "created files", "simulated_zeta.csv"), header = F)))
  S <- data.frame(mu = t(read.csv(here::here("matlab", "created files", "simulated_S.csv"), header = F)))
  
  mr_ph %>%
    mutate(zeta = zeta$zeta, S = S$mu) %>%
    ggplot(aes(x = S, y = zeta, col = lc)) +
    geom_point() +
    theme_classic()
  
  
  #the rescorla wagner model never "won" the model selection when the data_simulating model was sutton K1
  
  if (nrow(su1) == 2) {
    su1 <- rbind(su1, c("ph", 0), c("rw",0)) %>% arrange(lc)
  }
  
  if (nrow(ph) == 2) {
    ph <- rbind(ph, c("rw", 0),c("hgf",0)) %>% arrange(lc)
  }
  
  if (nrow(rw) == 3) {
    rw <- rbind(rw,c("ph",0)) %>% arrange(lc)
  }
  hgf = hgf %>% arrange(lc)
  rw = rw %>% arrange(lc)
  
  #combine the number of times each model won on each of the simulated data models 
  table <- cbind(hgf,ph, rw, su1)
  #remove duplicated columns such that we end up with a 3x3 matrix of simulated vs fitted models
  table <- table[, !duplicated(colnames(table))]
  table <- table %>% dplyr::rename(" " = lc)
  
  #renaming
  row.names(table) <- c("Recovered HGF", "Recovered RW", "Recovered SU1","Recovered PH")
  table$` ` <- NULL
  #renaming
  colnames(table) <- c("HGF", "RW", "Sutton", "pearce hall")
  #making the actual table
  table <- table %>% rownames_to_column(var = "  ")
  table$Sutton <- as.numeric(table$Sutton)
  ft <- flextable(table)
  ft <- add_header_row(ft, values = c(" ", "Simulated", " ", " "), colwidths = c(2, 1, 1,1))
  
  # Adjusting font size and column widths
  ft <- ft %>%
    fontsize(size = 10, part = "all") %>%
    width(j = 1:4, width = 1.5)
  return(ft)
}

#function to do model comparison given the log likelihoods of the actual model fitted to the real data
get_model_comparison <- function() {
  #read the data
  modelcompar <- data.frame(read.csv(here::here("matlab", "created files", "model_comparison2.csv"), header = F))
  
  
  #some renamning and tidying
  modelcompar$model <- c("Hierarchical \nGaussian\nFilter(2-level)", "Rescorla\nWagner", "Sutton\nK1", "Pearce\nhall")
  
  names(modelcompar)[1] <- "Model frequencies"
  names(modelcompar)[2] <- "Exceedance probability"
  
  #here the 267 is the number of participants in this analysis
  modelcompar$`Model frequencies` <- modelcompar$`Model frequencies` * 267
  modelcompar$`Exceedance probability` <- modelcompar$`Exceedance probability` * 267
  
  df <- reshape2::melt(modelcompar, id.vars = "model")
  
  model_compar = ggplot(df, aes(x = model, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(name = "Model Frequency", sec.axis = sec_axis(~ . / 267, name = "Exceedance probability")) +
    scale_fill_manual(values = c("#c44e52", "#5f9e6e")) +
    theme(axis.text.y = element_text(color = "red")) +
    guides(fill = guide_legend(override.aes = list(size = 8)))+
    theme(axis.title.x = element_blank())+
    theme+text+
    theme(legend.position = c(.95, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=font_size_small),
          axis.text=element_text(size=font_size_small),
          axis.title=element_text(size=font_size),
          axis.line=element_line(size=axis_width,colour = "black"),
          axis.ticks=element_line(size=axis_width),
          axis.title.y = element_text(color = "#c44e52"),
          axis.text.y = element_text(color = "#c44e52"),
          axis.title.y.right = element_text(color = "#5f9e6e"),
          axis.text.y.right = element_text(color = "#5f9e6e"))
  
  
  
  return(model_compar)
  
  
}



#function take takes the models used in the study and makes a table of all regression coefficients. (These are saved in xlsx files)
get_main_tables <- function(model, round = 2) {
  
  #retrieve the fixed effects of the model
  
  if("Component" %in% names(data.frame(parameters::model_parameters(model)))){
    
    fixedeffecs <- parameters::model_parameters(model) %>%
      mutate(CI = NULL, CI_low = NULL, CI_high = NULL, df_error = NULL) %>%
      dplyr::rename(parameter = Component) %>%
      dplyr::select(parameter, everything()) %>% 
      mutate(parameter = ifelse(str_detect(parameter, "conditional"), "μ", ifelse(str_detect(parameter, "sigma"), "σ", ifelse(str_detect(parameter, "tau"), "τ", "ν"))))
    #renaming
    names(fixedeffecs) <- c("parameter","contrast", "\u03B2", "SE", "t", "p")
    #formular for the model (i.e. the math)
    formular <- as.character(formula(model))
  }else{
    
    fixedeffecs <- parameters::model_parameters(model) %>%
      mutate(CI = NULL, CI_low = NULL, CI_high = NULL, df_error = NULL, Component = "conditional") %>%
      dplyr::rename(parameter = Component) %>%
      dplyr::select(parameter, everything()) %>% 
      mutate(parameter = ifelse(str_detect(parameter, "conditional"), "μ", ifelse(str_detect(parameter, "sigma"), "σ", ifelse(str_detect(parameter, "tau"), "τ", "ν"))))
    #renaming
    names(fixedeffecs) <- c("parameter","contrast", "\u03B2", "SE", "t", "p")
    #formular for the model (i.e. the math)
    formular <- as.character(formula(model))
    
    
  }
  #get family
  family = family(model)[2]
  link = model$mu.link
  
  if(family == "Beta Inflated"){
    family = "ZOIB"
  }
  
  if(family != "ZOIB"){
    coefnames <- fixedeffecs$contrast
    
    
    name_mapping <- c(
      "expectedN" = "Expectation(Neutral)",
      "expectedUE" = "Expectation(Unpredicted)",
      "expectedTGI" = "Expectation(TGI)",
      "coolness" = "Perceived Coldness",
      "contingencyPredicts warm" = "Contingency predicts(Warm)",
      "trial" = "Trialnumber",
      "coolness:contingencyPredicts warm" = "Perceived Coldness:\nContingency_predicts(Warm)",
      "sa1hat" = "prediction uncertainty"
    )
    
    # Create a new vector of row names by applying the name mapping to the current row names
    new_coefnames <- sapply(coefnames, function(x) {
      if (x %in% names(name_mapping)) {
        return(name_mapping[[x]])
      } else {
        return(x)
      }
    })
    
    #add the new names
    fixedeffecs$contrast <- new_coefnames
    
  }
  
  
  #big renaming of columns: for all cases of reference for stimulus
  # reference of TGI
  if (sum(grepl("TGI", fixedeffecs$contrast)) == 0) {
    fixedeffecs <- fixedeffecs %>% mutate(contrast = case_when(
      contrast == "stimwarm" ~ "Stimulus(Warm)",
      contrast == "stimcold" ~ "Stimulus(Cold)",
      contrast == "trial" ~ "Trialnumber",
      contrast == "sa2" ~ "Estimation uncertainty",
      contrast == "sa2:stimcold" ~ "Estimation uncertainty:\nStimulus(Cold)",
      contrast == "sa2:stimwarm" ~ "Estimation uncertainty:\nStimulus(Warm)",
      contrast == "Ratingscaleburnbeta" ~ "Ratingscale(Burning)",
      contrast == "Ratingscalewarmbeta" ~ "Ratingscale(Warm)",
      contrast == "belief_to_cold:Ratingscaleburnbeta" ~ "belief to cold:\nRatingscale(Burning)",
      contrast == "belief_to_cold:Ratingscalewarmbeta" ~ "belief to cold:\nRatingscale(Warm)",
      contrast == "(Intercept)" ~ "Intercept",
      contrast == "RateConCon" ~ "RatingScale(Factual)",
      contrast == "predRespcold" ~ "Prediction(Cold)",
      contrast == "belief_to_cold:stimwarm" ~ "belief to cold:\nStimulus(Warm)",
      contrast == "RateConCon:stimcold" ~ "RatingScale (Factual):\nStimulus(Cold)",
      contrast == "RateConCon:predRespcold" ~ "RatingScale (Factual):\nPrediction(Cold)",
      contrast == "stimcold:predRespcold" ~ "Stimulus (Cold):\nPrediction(Cold)",
      contrast == "RateConCon:stimcold:predRespcold" ~ "RatingScale(Factual):\n Stimulus (Cold):\n Prediction(Cold)",
      TRUE ~ contrast
    ))
    references <- c(
      "stimulus" = "TGI",
      "RatingScale" = "Counterfactual",
      "Prediction" = "Warm"
    )
  }
  # reference of cold
  if (sum(grepl("cold", fixedeffecs$contrast)) == 0) {
    fixedeffecs <- fixedeffecs %>% mutate(contrast = case_when(
      contrast == "stimwarm" ~ "stimulus(Warm)",
      contrast == "stimTGI" ~ "stimulus(TGI)",
      contrast == "trial" ~ "trialnumber",
      contrast == "Ratingscaleburnbeta" ~ "Ratingscale(Burning)",
      contrast == "Ratingscalewarmbeta" ~ "Ratingscale(Warm)",
      contrast == "belief_to_cold:Ratingscaleburnbeta" ~ "belief to cold:\nRatingscale(Burning)",
      contrast == "belief_to_cold:Ratingscalewarmbeta" ~ "belief to cold:\nRatingscale(Warm)",
      contrast == "(Intercept)" ~ "Intercept",
      contrast == "RateConCon" ~ "Rating Scale(Factual)",
      contrast == "predRespcold" ~ "Prediction(Cold)",
      contrast == "belief_to_cold:stimwarm" ~ "belief to cold:\nStimulus(Warm)",
      contrast == "RateConCon:stimcold" ~ "RatingScale (Factual):\nStimulus (Cold)",
      contrast == "RateConCon:predRespcold" ~ "RatingScale (Factual):\nPrediction(Cold)",
      contrast == "stimcold:predRespcold" ~ "Stimulus(Cold):\nPrediction(Cold)",
      contrast == "RateConCon:stimcold:predRespcold" ~ "RatingScale(Factual):\nStimulus(Cold):\nPrediction(Cold)",
      TRUE ~ contrast
    ))
    references <- c(
      "stimulus" = "cold",
      "RatingScale" = "Counterfactual",
      "Prediction" = "Warm"
    )
  }
  # reference of warm
  if (sum(grepl("warm", fixedeffecs$contrast)) == 0) {
    fixedeffecs <- fixedeffecs %>% mutate(contrast = case_when(
      contrast == "stimcold" ~ "stimulus(Cold)",
      contrast == "stimTGI" ~ "stimulus (TGI)",
      contrast == "trial" ~ "trialnumber",
      contrast == "Ratingscaleburnbeta" ~ "Ratingscale(Burning)",
      contrast == "Ratingscalewarmbeta" ~ "Ratingscale(Warm)",
      contrast == "belief_to_cold:Ratingscaleburnbeta" ~ "belief to cold:\nRatingscale(Burning)",
      contrast == "belief_to_cold:Ratingscalewarmbeta" ~ "belief to cold:\nRatingscale(Warm)",
      contrast == "(Intercept)" ~ "Intercept",
      contrast == "RateConCon" ~ "Rating Scale(Factual)",
      contrast == "predRespcold" ~ "Prediction(Cold)",
      contrast == "RateConCon:stimcold" ~ "RatingScale(Factual):\nStimulus(Cold)",
      contrast == "RateConCon:predRespcold" ~ "RatingScale(Factual):\nPrediction (Cold)",
      contrast == "stimcold:predRespcold" ~ "Stimulus(Cold):\nPrediction(Cold)",
      contrast == "RateConCon:stimcold:predRespcold" ~ "RatingScale(Factual):\nStimulus(Cold):\nPrediction(Cold)",
      TRUE ~ contrast
    ))
    references <- c(
      "stimulus" = "warm",
      "RatingScale" = "Counterfactual",
      "Prediction" = "Warm"
    )
  }
  
  
  if(formular[2] == "value" & grepl("stim",formular[3]) & !grepl("sa2",formular[3])& !grepl("RateCon",formular[3])){
    formular[2] = "Rating(Coldbeta | Warmbeta | Burningbeta)"
  }
  
  if(formular[2] == "value" & grepl("stim",formular[3]) & grepl("sa2",formular[3])){
    formular[2] = "Burningbeta"
  }
  
  if(grepl("sa2", formular[3])){
    formular[3] = str_replace(formular[3],"sa2","Est. uncertainty")
  }
  
  if(grepl("sa1hat", formular[3])){
    formular[3] = str_replace(formular[3],"sa1hat","Pred. uncertainty")
  }
  
  #formating and rounding the numeric values:
  fixedeffecs[, 3:6] <- apply(fixedeffecs[, 3:6], 2, function(x) formatC(x, format = "g", digits = round))
  #the table
  ft <- flextable(fixedeffecs) %>%
    add_header_row(values = paste0(formular[2], formular[1], formular[3], ", ", family, "(link = ",link,")"), colwidths = c(ncol(fixedeffecs))) %>%
    #add_header_lines(values = title) %>%
    width(j = c(1, 3:ncol(fixedeffecs)), width = 1) %>%
    width(j = 2, width = 1.8) %>%
    fontsize(size = 10, part = "all") %>%
    theme_vanilla() %>%
    align(i = 1:2, j = NULL, align = "center", part = "header")
  
  
  return(ft)
}

