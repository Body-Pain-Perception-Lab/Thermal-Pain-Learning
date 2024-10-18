logis_sigmoid = function(expect,zeta){
  return(expect^zeta / (expect^zeta + (1-expect)^zeta))
}


make_conting = function(){
  # Source additional scripts for plotting and utilities
  source(here::here("scripts","plots.R"))
  source(here::here("scripts","utils.R"))
  
  # Retrieve the OSF authentication token if the file exists
  osf_file_path <- here::here("osf","osf.txt")
  if (file.exists(osf_file_path)) {
    osf_token <- read_lines(osf_file_path)[1]
  } else {
    stop("OSF token file does not exist!")
  }
  
  # Get data
  data <- get_data(osf_token = osf_token, rerun = FALSE)
  df <- data$data
  
  # Get HGF trajectories
  trajfeel2 = get_hgf(df, rerun = F, osf_token = osf_token)
  df = trajfeel2 %>% filter(predRT > 0.1)%>% mutate(id = as.factor(id),stim =as.factor(stim))
  
  df$u = ifelse(df$cue == "low-tone" & df$stim == "TGI", a <- 1-(df$vasResp_1/(df$vasResp_1+df$vasResp_2)),
                ifelse(df$cue == "high-tone" & df$stim == "TGI", a <- (df$vasResp_1/(df$vasResp_1+df$vasResp_2)),
                       ifelse(df$cue == "low-tone" & df$stim == "cold", a <- 0, 
                              ifelse(df$cue == "high-tone" & df$stim == "warm", a <- 0,
                                     ifelse(df$cue == "low-tone" & df$stim == "warm", a <- 1,
                                            ifelse(df$cue == "high-tone" & df$stim == "cold", a <- 1, a <- 0.1))))))
  
  
  ids = df %>% select(u,id) %>% drop_na() %>% group_by(id) %>% summarize(n = n()) %>% filter(n == 306)
  
  df = df %>% filter(id %in% c(300,281)) 
  
  
  u1 = df %>% filter(id == 300) %>% .$u
  cue1 = df %>% filter(id == 300) %>% mutate(cue = ifelse(cue == "high-tone",1,ifelse(cue == "low-tone",0,NA))) %>% .$cue
  u2 = df %>% filter(id == 281) %>% .$u
  cue2 = df %>% filter(id == 281) %>% mutate(cue = ifelse(cue == "high-tone",1,ifelse(cue == "low-tone",0,NA))) %>% .$cue
  
  
  write.csv(u1, here::here("Shiny","Data","contingency1.csv"))
  write.csv(u2, here::here("Shiny","Data","contingency2.csv"))
  return(list(u1,cue1,u2,cue2))
}


get_hgf = function(u,e_0,zeta,omega,pi2_0){
  
  trial = length(u)
  mu2 = array(NA, trial)
  pi2 = array(NA, trial)
  sa2 = array(NA, trial)
  expect = array(NA, trial)
  
  mu2[1] = brms::logit_scaled(e_0)
  pi2[1] = pi2_0
  
  for(i in 1:trial){
    
    expect[i] = brms::inv_logit_scaled(mu2[i])
    sa2[i] = 1/pi2[i]
    pi2[i+1] = (1/(1/pi2[i] + (omega))) + 1/(1/(expect[i] * (1 - expect[i])))
    
    mu2[i+1] = mu2[i]+(1/pi2[i])*(u[i]-expect[i])
    
  }
  
  resp = rbinom(length(expect),1,logis_sigmoid(expect,zeta))
  data = data.frame(expect = expect, resp = resp) %>% mutate(omega = omega,pi2_0 = pi2_0, e_0 = e_0, zeta = zeta, model = "hgf") %>% 
    mutate(trial = 1:n(), u = u, lr = sa2)
  
  return(data)
}
get_rw = function(u,e_0,zeta,alpha){
  trial = length(u)
  expect2 = array(NA, trial)
  expect = array(NA, trial)
  
  expect[1] = (e_0)
  
  for(i in 1:trial){
    expect2[i] = expect[i]
    expect[i+1] = expect[i]+(alpha)*(u[i]-expect[i])
  }
  resp = rbinom(length(expect2),1,logis_sigmoid(expect2,zeta))
  data = data.frame(expect = expect2, resp = resp) %>% mutate(alpha = alpha, e_0 = e_0, zeta = zeta, model = "rw") %>% 
    mutate(trial = 1:n(), u = u,lr = alpha)
  
  
  return(data)
}
get_ph = function(u,e_0,zeta,S,a_0){
  trial = length(u)
  expect = array(NA, trial)
  expect2 = array(NA, trial)
  
  expect[1] = (e_0)
  
  for(i in 1:trial){
    expect2[i] = expect[i]
    if(i == 1){
      expect[i+1] = expect[i]+(S)*(u[i]-expect[i]) * a_0
    }else{
      expect[i+1] = expect[i]+(S)*(u[i]-expect[i]) * abs((u[i-1] - expect[i-1]))
    }
  }
  resp = rbinom(length(expect2),1,logis_sigmoid(expect2,zeta))
  data = data.frame(expect = expect2, resp = resp) %>% mutate(S = S, a_0 = a_0, e_0 = e_0, zeta = zeta, model = "ph") %>% 
    mutate(trial = 1:n(), u = u, lr = S*abs((u - expect)))
  
  return(data)
}
get_su1 = function(u,e_0,zeta,mu,Rhat,h_0){
  
  trial = length(u)
  expect2 = array(NA, trial)
  expect = array(NA, trial)
  be = array(NA,trial)
  al = array(NA,trial)
  h = array(NA,trial)
  
  expect[1] = (e_0)
  be[1] = log(Rhat)
  h[1] = h_0
  
  for(i in 1:trial){
    
    expect2[i] = expect[i]
    
    be[i+1] = be[i] + mu * (u[i]-expect[i])*h[i]
    al[i] = exp(be[i]) / (Rhat + exp(be[i]))
    h[i+1] = (h[i] + al[i] * (u[i]-expect[i])) * max((1-al[i]),0)
    
    
    expect[i+1] = expect[i] + al[i] * (u[i]-expect[i]);
    
  }
  
  
  resp = rbinom(length(expect2),1,logis_sigmoid(expect2,zeta))
  data = data.frame(expect = expect2, resp = resp) %>% mutate(mu = mu, Rhat = Rhat, h_0 = h_0, e_0 = e_0, zeta = zeta, model = "su1") %>% 
    mutate(trial = 1:n(), u = u,lr = al)
  
  return(data)
}




make_plot = function(data, model_name){
  
  p1 = ggplot(data, aes(x = trial, y = expect)) +
    geom_line() +
    geom_point(aes(y = resp)) +
    labs(x = "Trial", y = "Value", title = paste0(model_name," Plot")) +
    theme_minimal()
  
  lr = ggplot(data, aes(x = trial, y = lr)) +
    geom_line() +
    labs(x = "Trial", y = "Value", title = ("Learning rate")) +
    theme_minimal()
  
  return(p1/lr)
  
}
