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
                                                             vjust = 0.5),
                            text=ggplot2::element_text(family=font))
#colors
breaks = 4
color_cold = "#99CCED" # "#5495CB"
color_warm = "#FB6A4A" # "#F90047"
color_tgi =  "#8549A7" # "#800CEA"
color_grey = "#B3B3B3"


#function to generate plot 1
plot1 = function(df){
  
  #contingency plot (B) in figure 1
  plot_contingency = function(data1){
    #colors
    col = c(color_cold, color_warm, color_tgi)
    
    #recode data
    data1$line = 1
    data1$stim <- factor(data1$stim, levels = c("cold", "warm", "TGI"))
    data1 = data1 %>% dplyr::rename(Stimulus = stim)
    data1$Stimulus = as.factor(data1$Stimulus)
    data1 = data1 %>% filter(Stimulus != "NaN")
    data1 = droplevels(data1)
    levels(data1$Stimulus) = c("Cold","Warm","TGI")
    
    
    #taking ID 312 as the example contigency space for display
    data = data1 %>% filter(id == 312)
    
    data$desired_prob = as.numeric(as.character(data$desired_prob))
    data$desired_prob = data$desired_prob*100
    
    #getting values for where to put lines that indicate stimulus (i.e. red = warm, blue = cold, purple = TGI)
    data$upper = ifelse(data$desired_prob < 50 , a <- data$desired_prob-5, ifelse(data$desired_prob > 50, b <- data$desired_prob+5, 55))
    
    data$lower = ifelse(data$desired_prob < 50 , a <- data$desired_prob-0.6, ifelse(data$desired_prob > 50, b <- data$desired_prob+0.6,data$desired_prob))
    
    size = 1
    
    #get the high and low tone icons
    img <- png::readPNG(here::here("Figures","Pictures", "low-tone.PNG"))
    lowtone <- rasterGrob(img, interpolate=TRUE)
    img <- png::readPNG(here::here("Figures","Pictures", "high-tone.PNG"))
    hightone <- rasterGrob(img, interpolate=TRUE)
    
    #coordinates
    lowtone_cord = -20
    stim_cord = 0.46 #placement of middle legend
    
    #plotting
    plot_contingency = data1 %>% filter(id == 312)%>% 
      ggplot(aes(x= trial, y = as.numeric(as.character(desired_prob))*100, group = line))+
      geom_line(aes(x = trial, y = 100-as.numeric(as.character(desired_prob))*100, linetype = "dotted"), linewidth = 1)+
      geom_segment(data = data %>% filter(cue == "high-tone"), aes(x = trial, y = lower, col = Stimulus, xend = trial, yend = upper), linewidth = 0.5)+
      geom_segment(data = data %>% filter(cue == "low-tone"), aes(x = trial, y = 100-lower, col = Stimulus, xend = trial, yend = 100-upper), linewidth = 0.5)+
      geom_line(aes(linetype = "Solid"), size = 1)+
      xlab("Trial")+
      ylab(" ")+
      scale_y_continuous(lim = c(0,100),breaks=seq(0,100,by = 25), labels = c("0","25","50","75","100"))+
      scale_x_continuous(lim = c(0,306),breaks=pretty_breaks(n = breaks))+
      scale_color_manual(name = "Stimulus", values = col, labels = c("Cold","Warm","TGI"))+
      annotation_custom(lowtone, xmin = lowtone_cord, xmax = lowtone_cord+60, ymin = 95, ymax = 105)+
      annotation_custom(hightone, xmin = lowtone_cord+50, xmax = lowtone_cord+110, ymin = 95, ymax = 105)+
      geom_segment(aes(x = lowtone_cord+90, xend = lowtone_cord+100, y = 99, yend = 99), linewidth = 1)+
      geom_segment(aes(x = lowtone_cord+40, xend = lowtone_cord+55, y = 99, yend = 99), linewidth = 1, linetype = "dotted")+
      guides(linetype = "none",color=guide_legend(direction='horizontal',key.spacing.x = unit(0.15,"cm")))+ #spaceing of legends
      theme+text+
      theme(legend.position = c(stim_cord,0.96),
            legend.key.width = unit(0.85, 'cm'), #size of bars
            plot.title = element_text(hjust = 0.5), 
            legend.title = element_blank(),
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width))
    
    plot_contingency
    
    return(plot_contingency)
  }
  
  # prediction reaction time plot (D) in figure 1.
  plot_rts_new = function(data1){
    
    #getting colors for expected, neutral, unexpected and TGI stimulus
    col = c("#9cd8df", "#019090","#006072", color_tgi)
    
    #recoding and naming
    data1$expected = as.factor(data1$expected)
    
    levels(data1$expected) = c("P","N","UP","TGI")
    
    data1 = data1 %>% dplyr::rename(Expectation = expected)
    
    data1$Expectation = as.factor(data1$Expectation)
    
    #calculating the group level differences for prediction reaction time (predRT2) for each condition and id
    ff = data1 %>% filter(Expectation != "NA") %>% 
      mutate(predRT2 = lead(predRT)) %>% 
      group_by(id,Expectation) %>% 
      summarize(predRT2mean = mean(predRT2, na.rm = T),
                predRT2sd = sd(predRT2, na.rm = T))
    
    #plotting
    rts_expected = ff %>% ggplot(aes(x = Expectation, y = predRT2mean, fill = Expectation))+
      gghalves::geom_half_point(aes(col = Expectation), shape = 20, show.legend = FALSE, side = "r", alpha = 0.25,
                      transformation = position_jitter(width = 0.05, height = 0, seed = 123))+
      geom_half_violin(side = "l", show.legend = FALSE, alpha = .60)+
      geom_half_boxplot(notch = TRUE, width = 0.15, alpha = .6, outlier.shape = NA, show.legend = FALSE, side = "l")+
      scale_fill_manual(values = col)+
      scale_color_manual(values = col)+
      ylab(label = "RT (s)")+
      xlab(label = " ")+
      scale_y_continuous(lim = c(0,2),breaks=pretty_breaks(n = breaks), labels = comma)+
      theme+text+
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
    
    
    return(rts_expected)
    
  }
  
  #Accuracy plot / error rates (C) in figure 1. 
  plot_acc = function(data1){
    
    #getting colors
    col = c("#9cd8df", "#019090", "#006072", color_tgi)
    
    #renaming and recoding
    data1 = data1 %>% dplyr::rename(Expectation = expected)
    data1$Expectation = as.factor(data1$Expectation)
    
    
    levels(data1$Expectation) = c("P","N","TGI","UP")
    
    
    #getting participant and condition-wise accuracy
    acc = data1 %>%  filter(Expectation != "TGI") %>% 
      group_by(id,Expectation, predAcc, .drop = FALSE) %>% 
      summarize(n = n()) %>% 
      filter(predAcc != "TGI" & Expectation != "TGI")
    
    #initializing a list to store the results
    placeholder = list()
    
    #going through the data frame by 2 such that we calculate the procent of correct answers in each of the different conditions
    for (i in seq(1,nrow(acc),by=2)){
      value=acc[i,4]/(acc[i+1,4]+acc[i,4])
      placeholder = rbind(placeholder,value,value)
    }
    
    #adding it to the dataframe
    acc$procent = placeholder$n
    
    #taking every second row of the dataframe as it contains all the information
    acc = acc[seq(1,nrow(acc), by = 2),]
    
    #given that we plot error-rates we get this by taking one minus the accuracy 
    acc$procent = 1-acc$procent
    
    #plotting
    acc_expected = acc %>% ungroup() %>%
      ggplot(aes(x = Expectation, y = procent, fill = Expectation))+
      geom_half_violin(side = "l", show.legend = FALSE, alpha = .60) +
      geom_half_point(aes(col = Expectation), shape = 20, show.legend = FALSE, side = "r", alpha = 0.25,
                      transformation = position_jitter(width = 0.05, height = 0, seed = 123))+
      geom_half_boxplot(notch = TRUE, width = 0.15, alpha = .6, outlier.shape = NA, show.legend = FALSE, side = "l")+
      scale_fill_manual(values = col)+
      scale_color_manual(values = col)+
      ylab(label = "Error rate (%)")+
      xlab(label = " ")+
      labs(color = " ", fill = " ")+
      coord_cartesian(ylim = c(0,1))+
      scale_y_continuous(breaks=seq(0,1,by = 0.25), labels = c("0","25","50","75","1"))+
      theme+text+
      theme(legend.position = "center",
            legend.box = "horizontal",
            legend.direction="horizontal",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            legend.spacing.x = unit(0.35, "cm"),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    
    return(acc_expected)
  }
  
  
  #making the plots
  plot_contingency1 = plot_contingency(df)
  rts_expected1 = plot_rts_new(df)
  acc_expected1 = plot_acc(df)
  
  
  #making the layout
  layout <- "
AAAA
BBBB
CCDD
"
  
  #making a placeholder before adding (A) in figure 1
  plotspacer = ggplot(data.frame())+theme_void()
  
  #adding all elements together
  plot1 = plotspacer+
    plot_contingency1+
    acc_expected1+rts_expected1+
    plot_layout(design = layout)+
    plot_annotation(tag_levels = "A")&
    patchtheme
  
  
  #extracing the y-lab for (B) in figure 1. 
  imfile = here::here("Figures","Pictures", "ylab_fig1.PNG")
  
  
  #adding the y-lab to (B)
  plot1 = ggdraw() + 
    draw_plot(plot1) + 
    draw_image(magick::image_read(imfile),
               scale = 0.22,
               y = 0.065,
               x = -0.44)
  
  #adding the experimental paradigm. (A) in figure 1
  plot1 = ggdraw() + 
    draw_plot(plot1) + 
    draw_image(magick::image_read(here::here("Figures","Pictures","paradigm.PNG")),
               y = 0.34,
               x = 0.03,
               scale = 0.88)
  
  
  plot1 <- plot1 & xlab(NULL)
  
  # Use the tag label as an x-axis label
  plot1 = wrap_elements(panel = plot1) +
    labs(tag = "Trial type") +
    theme(
      plot.tag = element_text(size = rel(1)),
      plot.tag.position = c(0.55,0.03)
    )
  
  
  
  #saving the image such that it can be read in the manuscript markdown!
  ggsave(filename = here::here("Figures", "figure1.png"),
         plot = plot1,
         width = 7.2,
         height = 7.2,
         units = "in",
         dpi = 600)
  
  return(plot1)
}

#function to generate plot 2
plot2 = function(data){
  
  #main effect of stimulus on burning ratings figure 2A.
  burn_main_effect = function(data){
    
    #recoding and renaming
    data = data %>% filter(stim != "NaN")
    data$stim = factor(data$stim, levels = c("TGI", "warm", "cold"))
    col = c(color_tgi, color_warm, color_cold)
    
    #putting the data in a long format and extracting only the burning ratings (vasResp_3)
    data = data %>% 
      pivot_longer(cols = c(vasResp_1,vasResp_2,vasResp_3)) %>% 
      filter(name == "vasResp_3")
    
    #extracting all participant and stimulus wise mean ratings for burning.
    ff = data %>% filter(expected != "NA",stim != "NaN") %>% 
      group_by(id,stim, name) %>% 
      summarize(value = mean(value, na.rm = T)) %>% 
      mutate(name = as.factor(name))
    
    #renaming the levels of the factors
    levels(ff$name) = c("Burning","","")
    levels(ff$stim) = c("TGI","Warm","Cold")
    
    #plotting
    maineffect_stim = ff %>% ggplot(aes(x = stim, y = value,fill = stim))+
      geom_half_point(aes(col = stim), 
                      shape = 20, show.legend = FALSE, side = "l", alpha = 0.25,
                      transformation = position_jitter(width = 0.05, height = 0, seed = 123))+
      geom_half_violin(side = "r", show.legend = FALSE, alpha = .60, scale = 'width')  +
      geom_half_boxplot(notch = TRUE, width = 0.15, alpha = .6, outlier.shape = NA, show.legend = FALSE, side = "r")+
      scale_fill_manual(values = col)+
      scale_color_manual(values = col)+
      ylab("Burning Rating")+
      xlab(label = "Stimulus")+
      labs(color = "Stimulus", fill = "Stimulus")+
      scale_fill_manual(values = col)+
      scale_color_manual(values = col)+
      coord_flip()+theme+text+
      theme(legend.position = "center",
            legend.box = "horizontal",
            legend.direction="horizontal",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            legend.spacing.x = unit(0.35, "cm"),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text.y = element_blank())
    
    maineffect_stim
    
    return(maineffect_stim)
  }
  
  #three way interaction effect of expectations figure 2B.
  expectation_plot2 = function(df){
    #colors
    cols = c("#C51C1C", "#13216B", "#5363AC", "#FC6767")
    col = c(color_cold, color_warm, color_tgi)
    
    #first plotting the veridical ratings (cold stimulus on cold ratings and warm stimulus on warm ratings) across predictions (cold or warm)
    expect_plot_con = df %>% pivot_longer(c(vasResp_1,vasResp_2,vasResp_3)) %>% #filter(value != 0,value != 100) %>% 
      filter(stim != "TGI", predResp != "NaN", name != "vasResp_3")%>% 
      group_by(stim, name, predResp) %>% 
      summarize(mean = mean(value, na.rm =T), se = sd(value, na.rm =T)/sqrt(n())) %>%
      #making the congruency / veridical rating to stimulus
      mutate(RateCon = ifelse(stim == "cold" & name == "vasResp_1" | stim == "warm" & name == "vasResp_2"  , "Con", "Incon")) %>%
      mutate(stim = ifelse(stim == "cold","Cold",ifelse(stim == "warm","Warm",NA))) %>% 
      mutate(predResp = ifelse(predResp == "cold","Cold",ifelse(predResp == "warm","Warm",NA))) %>% 
      #taking only the congurrent (veridical)
      filter(RateCon == "Con") %>% 
      #plotting
      ggplot()+
      ylab(label = "VAS Rating")+
      labs(fill = "Prediction")+
      geom_bar(aes(x = stim,y = mean, fill = predResp), show.legend = TRUE,
               stat = "identity", position = "dodge",size=1, width = .75, color = "black", alpha = .6)+
      geom_errorbar(aes(x = stim, y = mean, group = predResp, ymin = mean-2*se, ymax = mean+2*se),
                    width = 0, position = position_dodge(0.75), size = 1)+
      scale_fill_manual(values = c(col[1],col[2]))+
      coord_cartesian(ylim = c(40,50))+
      scale_y_continuous(breaks=pretty_breaks(n = breaks), labels = comma) +
      #ggtitle("Factual\nRatings")+
      annotate("text",x = 1.6, y = 49, label = "Factual\nRatings", family = font, size = (font_size_small-1)/2) +
      xlab("Cold  Warm")+
      theme+text+
      theme(legend.position = "top",
            legend.box = "horizontal",
            legend.direction="horizontal",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            legend.spacing.x = unit(0.35, "cm"),
            axis.title.x = element_text(vjust = -0.4),
            text=element_text(family=font),
            #plot.title = element_text(hjust = 0.5, size = font_size,vjust = -27), #change legend key height
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_blank(),
            plot.margin = margin(r = 0  # Right margin
            ))
    
    
    
    
    #same for the incongurrent or counterfactual case i.e. cold stimulus on warm ratings and warm stimulus on cold ratings!
    
    expect_plot_incon = df %>% pivot_longer(c(vasResp_1,vasResp_2,vasResp_3)) %>% 
      filter(stim != "TGI", predResp != "NaN", name != "vasResp_3")%>% 
      group_by(stim, name, predResp) %>% summarize(mean = mean(value, na.rm =T), se = sd(value, na.rm =T)/sqrt(n())) %>%
      mutate(RateCon = ifelse(stim == "cold" & name == "vasResp_1" |stim == "warm" & name == "vasResp_2"  , "Con", "Incon")) %>%
      mutate(stim = ifelse(stim == "cold","Cold",ifelse(stim == "warm","Warm",NA))) %>% 
      mutate(predResp = ifelse(predResp == "cold","Cold",ifelse(predResp == "warm","Warm",NA))) %>% 
      #only inconcurrent (counter factual)
      filter(RateCon == "Incon") %>% 
      ggplot()+
      ylab(label = " ")+
      labs(fill = "Prediction")+
      geom_bar(aes(x = stim,y = mean, fill = predResp), show.legend = TRUE,
               stat = "identity", position = "dodge",size=1, width = .75, color = "black", alpha = .6)+
      geom_errorbar(aes(x = stim, y = mean, group = predResp, ymin = mean-2*se, ymax = mean+2*se),
                    width = 0, position = position_dodge(0.75), size = 1)+
      scale_fill_manual(values = c(col[1],col[2]))+
      coord_cartesian(ylim = c(10,20))+
      scale_y_continuous(breaks=pretty_breaks(n = breaks), labels = comma) +
      #ggtitle("Counterfactual\nRatings")+
      annotate("text",x = 1.6, y = 19, label = "Counterfactual\nRatings", family = font, size = (font_size_small-1)/2) +
      xlab("Cold  Warm")+
      theme+text+
      theme(legend.position = "top",
            legend.box = "horizontal",
            legend.direction="horizontal",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            legend.spacing.x = unit(0.35, "cm"),
            text=element_text(family=font),
            axis.title.x = element_text(vjust = -0.4),
            #plot.title = element_text(hjust = 0.5, size = font_size,vjust = -27), #change legend key height
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_blank(),
            plot.margin = margin(r = 0  # Right margin
            ))
    
    #combing the two and returning them
    this_plot <- expect_plot_con + expect_plot_incon + plot_layout(guides = "collect") & theme(legend.position = 'top')
    
    return(this_plot)
    
  }
  
  #Individual differences plot figure 2C
  individual_dif = function(df){
    
    #extracting TGI trials and summarizing across mean and standard errors for the cold and warm ratings across participants
    test = df %>% filter(stim == "TGI") %>% 
      pivot_longer(c(vasResp_1, vasResp_2,vasResp_3)) %>% 
      filter(name != "vasResp_3") %>% 
      group_by(id,name) %>% 
      summarize(n = n(),mean = mean(value, na.rm = T), se = sd(value, na.rm =T )/sqrt(n)) %>% 
      pivot_wider(names_from = name, values_from = c(mean,se))
    
    #getting the difference in mean rating of cold and warm
    test$color = test$mean_vasResp_1-test$mean_vasResp_2
    
    #getting the standard error on this difference
    test$color_se = sqrt(test$se_vasResp_1^2+test$se_vasResp_2^2)
    
    #pivoting the mean responses into a long format such that its easier to plot
    test = test %>% pivot_longer(c(mean_vasResp_1, mean_vasResp_2), names_to = "mean", values_to = "Rating")
    
    #making a color column that depends on whether the 95% confidence interval crosses 0 or not.
    test$col = ifelse(test$color-2*test$color_se > 0, "Cold",
                      ifelse(test$color+2*test$color_se < 0, "Warm", "Ambiguous"))
    
    #naming
    test$col = factor(test$col, levels = c("Cold", "Ambiguous", "Warm"))
    
    #colors for the points
    cols = c(color_cold, color_grey, color_warm)
    
    #ordering the ids to make the plot look nice
    test1 = test[order(test$color),]
    #making new id column
    test1$id2 = sort(rep(seq(1,267,by=1),2))
    
    #plotting (looks werid when not combined)
    feel_plot = test1 %>% ggplot(aes(x = id2, y = color)) +
      coord_flip() +
      geom_pointrange(aes(x = id2, y = -color, ymin = -color-2*color_se, ymax = -color+2*color_se, col = col), alpha = 0.3, size = 0.1) +
      geom_hline(yintercept = 0, size = axis_width, linetype = 2) +
      scale_color_manual("Perceived quality\nof the TGI", values = cols) +
      scale_y_continuous("Mostly Cold              Mostly Warm", breaks = pretty_breaks(n = 6), labels = comma) + 
      scale_x_continuous("Participant")+
      labs(color = "Perceived quality\nof the TGI")+
      theme+text+
      guides(color = guide_legend(override.aes = list(linetype = "blank", size = 1)))+
      theme(legend.position = c(0.8, 0.83),
            axis.title.y = element_text(margin = margin(r = -59)),
            legend.box = "vertical",
            legend.direction="vertical",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.title = element_text(size=font_size_small+2), #change legend title font size
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    
    return(feel_plot)
  }
  
  
  #three way / two way interaction plot figure 2D
  threewayint2 = function(df){
    
    #colors
    col = c("snow3","slategray")
    
    #labels for cues.
    labels <- c(
      "low-tone" = "Low Tone",
      "high-tone" = "High Tone"
    )
    
    #load the necessary for the plot
    
    base::load(here::here("Analysis", "Workspace", "percieved_TGI.RData"))
    
    
    num_colors <- 5000
    
    # Create a color vector using colorRampPalette
    color_vector <- colorRampPalette(c("blue","blue","white","white","red","red"))(num_colors)
    
    
    percieved_TGI_acc_contingency_cold_95 = data.frame(percieved_TGI_acc_contingency_cold_95) %>% 
      mutate(x = 1-x)
    
    percieved_TGI_acc_contingency_cold_80 = data.frame(percieved_TGI_acc_contingency_cold_80) %>% 
      mutate(x = 1-x)
    
    percieved_TGI_acc_contingency_cold_50 = data.frame(percieved_TGI_acc_contingency_cold_50) %>% 
      mutate(x = 1-x)
    
    levels(percieved_TGI_acc_contingency_cold_95$group) = c("Predicts cold", "Predicts warm")
    #plotting the results
    predacc_cool =percieved_TGI_acc_contingency_cold_95 %>% 
      ggplot(aes(x = x, y = predicted*100, fill = group))+
      geom_line(col = "black")+
      geom_ribbon(data = percieved_TGI_acc_contingency_cold_95, aes(ymin = conf.low*100, ymax = conf.high*100), alpha = 0.15, colour = NA)+
      geom_ribbon(data = percieved_TGI_acc_contingency_cold_80, aes(ymin = conf.low*100, ymax = conf.high*100), alpha = 0.25, colour = NA)+
      geom_ribbon(data = percieved_TGI_acc_contingency_cold_50, aes(ymin = conf.low*100, ymax = conf.high*100), alpha = 0.35, colour = NA)+
      geom_line(data = data.frame(x = rep(seq(0, 1, length.out = 101),2), group = c(rep("Predicts cold",101),rep("Predicts warm",101))),aes(x = x, y = 37, color = x), linewidth = 4)+
      coord_cartesian(ylim = c(37,70), xlim = c(0,1))+
      scale_color_gradientn(colors = color_vector, guide = "none")+
      #scale_colour_gradient(low = "red", high = "blue", guide = "none")+
      scale_fill_manual(" ",values = c(col[1],col[2]), labels = c("Contingency predicts Cold","Contingency predicts Warm"))+
      theme(axis.ticks.x = element_blank(),panel.spacing.x=unit(1.5, "lines"))+
      scale_x_continuous("Perceived quality of the TGI", breaks =c(0,1), labels = c("Cold","Warm"))+
      scale_y_continuous("P(Correct on following trial) (%)", breaks=pretty_breaks(n = breaks), labels = comma)+
      theme+text+
      theme(legend.position = c(0.5,0.9),
            legend.box = "vertical",
            legend.direction="vertical",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    
    predacc_cool
    
    return(predacc_cool)
    
  }
  
  
  #layout for figure 2
  layout <- "
AB
CD
"
  #getting the plots for figure 2
  main = burn_main_effect(data)
  expectation = expectation_plot2(data)
  individual = individual_dif(data)
  threeway = threewayint2(df)
  
  
  #expectation <- this_plot
  
  # Use the tag label as an x-axis label
  
  #making figure 2 with the layout
  plot = main+
    expectation+
    individual+
    threeway+
    plot_layout(design = layout)+plot_annotation(tag_levels = list(c("A", "B", " ", "C", "D"))) &
    patchtheme
  
  
  image <- image_read(here("Figures","Pictures", "Stimuli.PNG"))
  rotated_image <- image_rotate(image, 90)  # Rotate the image by 45 degrees
  image <- image_read(here("Figures","Pictures", "Stimuli_ver.png"))
  # rotated_image <- image_rotate(image, 90)  # Rotate the image by 45 degrees
  
  # for main effect of stim
  
  plot1 = ggdraw() + 
    draw_plot(plot) + 
    draw_image(image,
               y = -0.364,
               x = -0.422,
               height = 2.2,
               scale = 0.122)
  
  
  image <- image_read(here("Figures","Pictures", "cold_warm_pair.PNG"))
  
  plot2 = ggdraw() + 
    draw_plot(plot1) + 
    draw_image(image,
               y = 0.025,
               x = 0.167,
               scale = 0.09)
  
  plot2
  
  
  plot3 = ggdraw() + 
    draw_plot(plot2) + 
    draw_image(image,
               y = 0.025,
               x = 0.392,
               scale = 0.09)
  
  plot3
  
  
  
  plot1 = wrap_elements(panel = plot3) +
    labs(tag = "Stimulus") +
    theme(
      plot.tag = element_text(size = font_size),
      plot.tag.position = c(0.78,0.47)
    )
  plot1
  #sacing it for the manuscript markdown
  ggsave(filename = here::here("Figures", "figure2.png"), plot = plot1, width = 7.2, height = 7.2, units = "in", dpi = 600)
  #returning it for the plot markdown.
  return(plot1)
  
  
  
  this_plot = wrap_elements((expect_plot_con + expect_plot_incon) +plot_annotation(title = 'A'))
  
  
}

#function to generate plot 3
plot3 = function(data){
  
  #plot of HGF trajectories figure 3B
  plot_trajectories = function(trajfeel2){
    
    #extracing a participants data
    df = trajfeel2 %>% filter(id == 594)
    
    #getting the cue-stimulus contingency i.e. the input sequence of the task
    df$u = ifelse(df$cue == "low-tone" & df$stim == "TGI", a <- 1-(df$vasResp_1/(df$vasResp_1+df$vasResp_2)),
                  ifelse(df$cue == "high-tone" & df$stim == "TGI", a <- (df$vasResp_1/(df$vasResp_1+df$vasResp_2)),
                         ifelse(df$cue == "low-tone" & df$stim == "cold", a <- 0, 
                                ifelse(df$cue == "high-tone" & df$stim == "warm", a <- 0,
                                       ifelse(df$cue == "low-tone" & df$stim == "warm", a <- 1,
                                              ifelse(df$cue == "high-tone" & df$stim == "cold", a <- 1, a <- 0.1))))))
    
    #getting the participants' reponses in the contingency space
    df$predResp2 = ifelse(df$cue == "low-tone" & df$predResp == "cold", a <- 0, 
                          ifelse(df$cue == "high-tone" & df$predResp == "warm", a <- 0,
                                 ifelse(df$cue == "low-tone" & df$predResp == "warm", a <- 1,
                                        ifelse(df$cue == "high-tone" & df$predResp == "cold", a <- 1, a <- NaN))))
    
    
    #extracting the first level latent parameters of the model.
    q1 = df %>% dplyr::select(desired_prob,u,predResp2,trial,mu1hat,sa1hat,mu2,sa2) %>% 
      mutate(level = 1, mu2 = NA, sa2 = NA)
    
    #extracting the second level latent parameters of the model
    q2 = df %>% dplyr::select(desired_prob,u,predResp2,trial,mu1hat,sa1hat,mu2,sa2) %>% 
      mutate(level = 2, desired_prob = NA, mu1hat = NA, u = NA, predResp2 = NA)
    
    #combing the two
    q3 = rbind(q1,q2)
    q3$level = as.factor(q3$level)
    
    #blue, #purple #red
    col = c(color_cold, color_tgi, color_warm)
    
    #green, grey, purple, yellow, 
    col = c("#8dcad1", "#9FADBD", color_tgi, "#f3d8a5")
    
    
    
    #different sizing arguments
    p_size = 0.5
    
    p_adj = 0.02
    #adjust the input sequence such that responses and inputs do not overlap
    q3$u = ifelse(q3$u == 0 , a <- q3$u-p_adj, b <- q3$u+p_adj)
    #same for responses
    q3$predResp2 = ifelse(q3$predResp2 == 0 , a <- q3$predResp2+p_adj,b <- q3$predResp2-p_adj)
    
    
    
    
    #plotting
    trajectories = q3 %>% mutate(level = factor(level, labels = c("Prediction","Estimation")),
                                 level = relevel(level, ref = "Estimation"),
                                 lower1 = (mu1hat-sa1hat)*100,
                                 upper1 = (mu1hat+sa1hat)*100,
                                 mu1hat = mu1hat*100,
                                 predResp2 = predResp2*100,
                                 u = u*100,
                                 lower2 = mu2-sa2,
                                 upper2 = mu2+sa2,
                                 desired_prob = as.numeric(as.character(desired_prob))*100) %>% 
      dplyr::rename(Trial = trial) %>% 
      ggplot(aes())+
      
      ggh4x::facet_wrap2(~level,nrow = 2, scales = "free",
                  strip = strip_themed(background_x = list(element_rect(fill = "#9db2d4"),
                                                           element_rect(fill = "#df9ea0"))))+
      geom_point(aes(x = Trial, y = predResp2+0.05, col = "#8dcad1"), size = p_size, show.legend = FALSE)+
      geom_point(aes(x = Trial, y = u-0.05, col = "#f3d8a5"), size = p_size, show.legend = FALSE)+
      geom_line(aes(x = Trial, y = mu1hat, col = "a"), show.legend = FALSE)+
      geom_ribbon(aes(x = Trial, ymax = upper1, ymin = lower1, alpha = 0.75, fill = "#9db2d4"), show.legend = FALSE)+
      geom_ribbon(aes(x = Trial, ymax = upper2, ymin = lower2, alpha = 0.75, fill  = "#df9ea0"), show.legend = FALSE)+
      geom_line(aes(x = Trial, y = desired_prob), size = 1)+
      geom_line(aes(x = Trial, y = mu2, col = "d"), show.legend = FALSE)+
      scale_color_manual(values = c("#8dcad1","#f3d8a5","#c44e52", "#366dc7"))+
      scale_fill_manual(values = c("#df9ea0","#9db2d4"))+
      ylab(expression("x"[1]~"(%)                    x"[2]~"(a.u.)"))+
      theme+text+
      coord_cartesian(clip = "off")+
      theme(legend.position = c(0.5,0.9),
            legend.box = "vertical",
            legend.direction="vertical",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
            strip.background = element_rect(fill = "white", colour = "black", linewidth = axis_width)
      )
    
    
    
    trajectories
    
    
    return(trajectories)
  }
  
  #barplot of the reaction times as a function of the prediction uncertainty figure 3D
  RT_uncertainty = function(trajeel2){
    
    #we bin the prediction uncertainty to make a barplot
    trajfeel2$bins = cut(trajfeel2$sa1hat, breaks = seq(0,0.25, by = 0.025))
    
    #calculate the mean and standard error for the bins
    qq = trajfeel2 %>% group_by(bins) %>% 
      summarize(n = n(), mean = mean(predRT,na.rm = T),
                se = sd(predRT,na.rm = T)/sqrt(n))
    
    #plot
    rts = qq %>% mutate(filler = as.numeric(qq$bins))%>% 
      ggplot(aes(x = bins, y = mean, fill = filler))+
      geom_errorbar(aes(y = mean,ymin = mean-2*se, ymax=mean+2*se), width = 0.2, col = "black")+
      geom_bar(stat = "summary",position = "dodge", fun = "mean")+
      coord_cartesian(ylim=c(0.75,1.11))+
      scale_y_continuous("RT (s)", breaks = seq(0.8,1.1,by = 0.1), labels = c("0.8","0.9","1.0","1.1"))+
      scale_x_discrete(name = "Prediction Uncertainty",
                       labels = c("Low","","","","","","","","High"),
                       limits = c("(0,0.025]","(0.025,0.05]","(0.05,0.075]","(0.075,0.1]","(0.1,0.125]","(0.125,0.15]","(0.175,0.2]","(0.2,0.225]","(0.225,0.25]"))+
      scale_fill_gradient(low = "#9cd8df", high = "#006072")+
      theme+text+
      theme(legend.position = "none",
            legend.box = "vertical",
            legend.direction="vertical",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
            strip.background = element_rect(fill = "white", colour = "white"))
    
    
    
    return(rts)
  }
  
  #barplot of the accuracy / errorrates as a function of the prediction uncertainty figure 3C
  errorrates_uncertainty = function(trajfeel2){
    
    #refactoring and making the bins
    trajfeel2$predAcc = as.character(trajfeel2$predAcc)
    trajfeel2$id = as.integer(trajfeel2$id)
    trajfeel2$bins = cut(trajfeel2$sa1hat, breaks = seq(0,0.25, by = 0.025))
    
    #first we get the number of trials in each bin for each participant in both correct and incorrect trials
    popper = trajfeel2 %>% filter(predAcc != "TGI", is.na(predAcc) == FALSE) %>% 
      mutate(id = as.factor(id), bins = as.factor(bins), predAcc = as.factor(predAcc)) %>%
      group_by(id,bins,predAcc) %>% 
      summarize(n = n())
    
    #next we get all combinations of these as some cases in the above would have all correct or incorrect trials in a bin.
    popper1 = trajfeel2 %>% filter(predAcc != "TGI", is.na(predAcc) == FALSE) %>%
      mutate(id = as.factor(id), bins = as.factor(bins), predAcc = as.factor(predAcc)) %>%
      ungroup() %>% 
      tidyr::expand(id,bins,predAcc)
    
    #now we can join the two too get a full data frame of possible options and fill the ones missing with 0 as this is here were a participant
    #could have all incorrect trials in a bin and therefore 0 correct.
    popper = full_join(popper,popper1) %>% 
      replace_na(list(n = 0))
    
    #ordering such that we can loop through them to get a procent of correct respones
    popper = popper[with(popper, order(id,bins)),]
    
    #initalize empty list
    placeholder = list()
    #loop through every other row and calculate the accuracy procent
    for (i in seq(1,nrow(popper),by=2)){
      value=popper[i,4]/(popper[i+1,4]+popper[i,4])
      placeholder = rbind(placeholder,value,value)
    }
    
    #now we add this accuracy procent to the original data frame, but given that we plot error rates we say 1-accuracy
    popper$procent = 1-placeholder$n
    
    #now extracting every other row as this contains all the information needed
    popper = popper[seq(1,nrow(popper), by = 2),]
    
    #getting the group level estimates by bin. Here multipling by 100 gives the procents instead of de decimals
    means = popper %>% group_by(bins) %>% 
      summarize(n = n(),meanp = mean(procent,na.rm= T)*100, sep = (sd(procent,na.rm= T)/(sqrt(n)))*100)
    
    #plotting
    errors  = means %>%  ggplot(aes(x = bins, y = meanp, fill = as.numeric(bins)))+
      geom_errorbar(aes(y = meanp,ymin = meanp-2*sep, ymax=meanp+2*sep), width = 0.2, col = "black")+
      geom_bar(stat = "summary",position = "dodge", fun = "mean")+
      #geom_point(aes(x = bins, y = procent*100),show.legend = FALSE,shape = 1,position = position_nudge(x = -0.05,y=0))+
      #ggdist::stat_slab(scale = 0.8, alpha = 0.8, show.legend = FALSE, interval_size = 0, point_size = 0)+
      ylab("Error rate (%)")+
      scale_x_discrete(name = " ", labels = c("Low uncertainty"," "," "," ", ""," "," "," ","High uncertainty"), limits = c("(0,0.025]","(0.025,0.05]","(0.05,0.075]","(0.075,0.1]","(0.1,0.125]","(0.125,0.15]","(0.175,0.2]","(0.2,0.225]","(0.225,0.25]"))+
      scale_fill_gradient(low = "#9cd8df", high = "#006072")+
      coord_cartesian(ylim=c(15,51))+
      theme+text+
      theme(legend.position = "none",
            legend.box = "vertical",
            legend.direction="vertical",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.title = element_blank(), #change legend title font size
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.background = element_rect(fill = "white", colour = "white"))
    
    
    return(errors)
  }
  
  
  #getting the plots:
  trajectory = plot_trajectories(data)
  RT = RT_uncertainty(data)
  errorrates = errorrates_uncertainty(data)
  #the layout for figure 3
  layout = c(patchwork::area(0,0,4,2),
             patchwork::area(0,3,2,4),
             patchwork::area(3,3,3,4),
             patchwork::area(4,3,4,4))
  
  #making a spacer for figure 3A
  plotspacer = data.frame() %>% 
    ggplot() + 
    theme_void()+theme(plot.tag.position  = c(0.05, 1))
  
  #adding the plots together
  plot1 = plotspacer+
    trajectory+
    errorrates+
    RT+
    plot_layout(design = layout)+
    plot_annotation(tag_levels = 'A')&
    patchtheme
  
  
  #draw on figure 3A
  plot1 = ggdraw() + 
    draw_plot(plot1) + 
    draw_image(magick::image_read(here::here("Figures","Pictures","figure3_180124_merged.tif")),
               scale = 0.94,
               x = -0.265,
               y = -0.009)
  plot1
  
  #save the plot for the markdown manuscript
  ggsave(filename = here::here("Figures", "figure3.png"), plot = plot1, width = 7.2, height = 7.2, units = "in", dpi = 600)
  
  #return the plot for the markdown plot
  return(plot1)
}


plot4_v2 = function(data){
  
  #Figure 4B
  main_sa2 = function(data2){
    #getting colors
    col = c(color_cold, color_warm,color_tgi)
    model = read.csv(paste0(here::here("Analysis","Plotting"),"/summarystats_umti.csv"))
    
    parameters = data.frame(coef = model[1:7,2:3]) %>% rownames_to_column("param") %>% drop_na() %>% 
      rename(est = coef.Estimate,
             se = coef.Std..Error)
    
    sa2 <- seq(-1, 4, by = 0.1)  # Continuous variable
    stim <- factor(rep(c("TGI", "Cold", "Warm"), each = length(sa2)))  # 3-level factor
    stim <- relevel(stim, ref = "TGI")
    trial <- 0  # Constant value for 'trial'
    sa2_expanded <- rep(sa2, times = 3)
    data <- data.frame(sa2 = sa2_expanded, stim = stim, trial = trial)
    design_matrix <- model.matrix(~ sa2 * stim + trial, data = data)
    
    
    yhat_mean = brms::inv_logit_scaled(design_matrix %*% parameters$est)
    yhat_q5 = brms::inv_logit_scaled(design_matrix %*% (parameters$est - 1.96  * parameters$se))
    yhat_q95 = brms::inv_logit_scaled(design_matrix %*% (parameters$est + 1.96 * parameters$se))
    yhat_q20 = brms::inv_logit_scaled(design_matrix %*% (parameters$est - 1.3 * parameters$se))
    yhat_q80 = brms::inv_logit_scaled(design_matrix %*% (parameters$est + 1.3 * parameters$se))
    yhat_q50l = brms::inv_logit_scaled(design_matrix %*% (parameters$est - 0.68 * parameters$se))
    yhat_q50h = brms::inv_logit_scaled(design_matrix %*% (parameters$est + 0.68 * parameters$se))
    
    
    mean = data.frame(pred = yhat_mean*100, sa2 = data$sa2, stim = data$stim)%>% 
      mutate(Rating = "Burning")
    
    q95 = data.frame(conf.high = yhat_q95*100, conf.low = yhat_q5*100, sa2 = data$sa2, stim = data$stim)%>% 
      mutate(Rating = "Burning")
    
    q80 = data.frame(conf.high = yhat_q80*100, conf.low = yhat_q20*100, sa2 = data$sa2, stim = data$stim)%>% 
      mutate(Rating = "Burning")
    
    q50 = data.frame(conf.high = yhat_q50h*100, conf.low = yhat_q50l*100, sa2 = data$sa2, stim = data$stim)%>% 
      mutate(Rating = "Burning")
    
    mean$stim = factor(mean$stim, levels = c("Cold","Warm","TGI"))
    q95$stim = factor(q95$stim, levels = c("Cold","Warm","TGI"))
    q80$stim = factor(q80$stim, levels = c("Cold","Warm","TGI"))
    q50$stim = factor(q50$stim, levels = c("Cold","Warm","TGI"))
    
    
    #plotting. Note for visual purposes we have limited the graph to be between 0 and 8.
    sa2_main = mean %>% 
      ggplot()+
      geom_line(aes(x = sa2, y = pred, group = stim),col = "black")+
      geom_ribbon(data = data.frame(q95), aes(x = sa2,ymin = conf.low, ymax = conf.high, fill = stim), alpha = 0.15)+
      geom_ribbon(data = data.frame(q80), aes(x = sa2,ymin = conf.low, ymax = conf.high, fill = stim), alpha = 0.25)+
      geom_ribbon(data = data.frame(q50), aes(x = sa2,ymin = conf.low, ymax = conf.high, fill = stim), alpha = 0.35)+
      scale_fill_manual(name = "Stimulus",values = col, labels = c("Cold","Warm","TGI"))+
      facet_wrap(~Rating, labeller = label_both)+
      scale_y_continuous("Predicted VAS", breaks=seq(15,35,by = 10), labels = seq(15,35,by = 10))+
      scale_x_continuous("Estimation Uncertainty", breaks =seq(-1,4,by = 5), labels = c("Low","High"))+
      theme+text+
      coord_cartesian(ylim = c(15,35),clip = "off")+
      theme(legend.position = "top",
            legend.box = "horizontal",
            legend.direction="horizontal",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(linewidth=axis_width),
            axis.ticks=element_line(linewidth=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
            strip.background = element_rect(fill = "white", colour = "black",linewidth = axis_width))
    
    sa2_main
    
    return(sa2_main)
  }
  
  #Figure 4A
  belief_tgi = function(data2){
    #colors
    col = c(color_cold, color_warm,color_tgi)
    
    #loading the data for the plot
    load(here::here("Analysis", "Workspace", "figure4a.RData"))
    
    gg_model1 = data.frame(figure4a_95)
    
    levels(gg_model1$facet) = c("Cold","Burning","Warm")
    
    gg_model1$group <- factor(gg_model1$group, levels = c("cold","warm","TGI"))
    
    gg_model2 = data.frame(figure4a_80)
    levels(gg_model2$facet) = c("Cold","Burning","Warm")
    gg_model2$group <- factor(gg_model2$group, levels = c("cold","warm","TGI"))
    
    
    gg_model3 = data.frame(figure4a_50)
    levels(gg_model3$facet) = c("Cold","Burning","Warm")
    gg_model3$group <- factor(gg_model3$group, levels = c("cold","warm","TGI"))
    
    
    num_colors <- 5000
    
    # Create a color vector using colorRampPalette
    color_vector <- colorRampPalette(c("blue","blue","white","white","red","red"))(num_colors)
    
    gg_model1 = data.frame(gg_model1) %>% mutate(x = 1-x)
    gg_model2 = data.frame(gg_model2) %>% mutate(x = 1-x)
    gg_model3 = data.frame(gg_model3) %>% mutate(x = 1-x)
    
    
    Cold = gg_model1 %>% 
      filter(facet == "Cold") %>%
      dplyr::rename(Ratingscale = facet) %>% 
      ggplot(aes(x = x, y = predicted,group = group))+
      geom_line(col = "black")+
      geom_ribbon(data = gg_model1 %>% filter(facet == "Cold")%>% dplyr::rename(Rating = facet), aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.15)+
      geom_ribbon(data = gg_model2 %>% filter(facet == "Cold")%>% dplyr::rename(Rating = facet), aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25)+
      geom_ribbon(data = gg_model3 %>% filter(facet == "Cold")%>% dplyr::rename(Rating = facet), aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.35)+
      geom_line(data = data.frame(x = seq(0, 1, length.out = num_colors), group = NA),aes(x = x, y = 0.24, color = x), linewidth = 3)+
      scale_fill_manual(name = "Stimulus",values = col, labels = c("Cold","Warm","TGI"))+
      scale_color_gradientn(colors = color_vector, guide = "none")+
      theme_classic()+
      facet_wrap2(~Rating, labeller = label_both)+
      scale_y_continuous("Predicted VAS", breaks=seq(0.3,0.6,by = 0.1), labels = seq(30,60,by = 10))+
      #scale_x_continuous("Prediction Uncertainty",breaks=pretty_breaks(n = 4), labels = comma,lim = c(0,1))+
      scale_x_continuous("Prediction Uncertainty", breaks =c(0,0.5,1), labels = c("Low","High","Low"))+
      guides(fill = guide_legend(title = "Stimulus"))+
      theme+text+
      coord_cartesian(ylim = c(0.25,0.60),clip = "off")+
      theme(legend.position = "top",
            legend.box = "horizontal",
            legend.direction="horizontal",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(linewidth=axis_width),
            axis.ticks=element_line(linewidth=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
            strip.background = element_rect(fill = "white", colour = "black",linewidth = axis_width))
    # plot.margin = margin(r = 0  # Right margin
    #                      )) # Left margin)
    
    
    Warm = gg_model1 %>% 
      filter(facet == "Warm") %>%
      dplyr::rename(Ratingscale = facet) %>% 
      ggplot(aes(x = x, y = predicted,group = group))+
      geom_line(col = "black")+
      geom_ribbon(data = gg_model1 %>% filter(facet == "Warm")%>% dplyr::rename(Rating = facet), aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.15)+
      geom_ribbon(data = gg_model2 %>% filter(facet == "Warm")%>% dplyr::rename(Rating = facet), aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25)+
      geom_ribbon(data = gg_model3 %>% filter(facet == "Warm")%>% dplyr::rename(Rating = facet), aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.35)+
      geom_line(data = data.frame(x = seq(0, 1, length.out = num_colors), group = NA),aes(x = x, y = 0.24, color = x), linewidth = 3)+
      scale_fill_manual(name = "Stimulus",values = col, labels = c("Cold","Warm","TGI"))+
      scale_color_gradientn(colors = color_vector, guide = "none")+
      theme_classic()+
      facet_wrap2(~Rating, labeller = label_both)+
      scale_y_continuous(" ", breaks=seq(0.3,0.6,by = 0.1), labels = seq(30,60,by = 10))+
      #scale_x_continuous("Prediction Uncertainty",breaks=pretty_breaks(n = 4), labels = comma,lim = c(0,1))+
      scale_x_continuous("Prediction Uncertainty", breaks =c(0,0.5,1), labels = c("Low","High","Low"))+
      guides(fill = guide_legend(title = "Stimulus"))+
      theme+text+
      coord_cartesian(ylim = c(0.25,0.60), clip = "off")+
      theme(legend.position = "top",
            legend.box = "horizontal",
            legend.direction="horizontal",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(linewidth=axis_width),
            axis.ticks=element_line(linewidth=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
            strip.background = element_rect(fill = "white", colour = "black",linewidth = axis_width))
    # plot.margin = margin(r = 0  # Right margin
    # ))
    
    
    
    
    Burning = gg_model1 %>% 
      filter(facet == "Burning") %>%
      dplyr::rename(Ratingscale = facet) %>% 
      ggplot(aes(x = x, y = predicted,group = group))+
      geom_line(col = "black")+
      geom_ribbon(data = gg_model1 %>% filter(facet == "Burning")%>% dplyr::rename(Rating = facet), aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.15)+
      geom_ribbon(data = gg_model2 %>% filter(facet == "Burning")%>% dplyr::rename(Rating = facet), aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.25)+
      geom_ribbon(data = gg_model3 %>% filter(facet == "Burning")%>% dplyr::rename(Rating = facet), aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.35)+
      geom_line(data = data.frame(x = seq(0, 1, length.out = num_colors), group = NA),aes(x = x, y = 0.195, color = x), linewidth = 3)+
      scale_fill_manual(name = "Stimulus",values = col, labels = c("Cold","Warm","TGI"))+
      # scale_colour_gradient2(low = "red",
      #                        mid = "white",
      #                        high = "blue",
      #                        midpoint = 0.5,
      #                        guide = "none")+
      scale_color_gradientn(colors = color_vector, guide = "none")+
      theme_classic()+
      facet_wrap2(~Rating, labeller = label_both)+
      scale_y_continuous(" ",breaks=pretty_breaks(n = 4), labels = number_format(accuracy = 0.01))+
      #scale_x_continuous("Prediction Uncertainty",breaks=pretty_breaks(n = 4), labels = comma,lim = c(0,1))+
      scale_x_continuous("", breaks =c(0,0.5,1), labels = c("Low","High","Low"))+
      guides(fill = guide_legend(title = "Stimulus"))+
      theme+text+
      coord_cartesian(ylim = c(0.1975,0.4),clip = "off")+
      theme(legend.position = "top",
            legend.box = "horizontal",
            legend.direction="horizontal",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(linewidth=axis_width),
            axis.ticks=element_line(linewidth=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
            strip.background = element_rect(fill = "white", colour = "black",linewidth = axis_width),
            plot.margin = margin(r = 0  # Right margin
            ))
    
    
    
    
    
    return(list(Cold,Warm,Burning))
  }
  
  #Figure 4C
  responder_plot_v2 = function(data2){
    
    responders <- data2 %>%
      filter(stim != "NaN") %>% 
      rename(VAS_Cold = vasResp_1,VAS_Warm = vasResp_2,VAS_Burn = vasResp_3,Stimulus = stim) %>% 
      pivot_longer(cols = c("VAS_Cold","VAS_Warm","VAS_Burn")) %>% 
      filter(name == "VAS_Burn") %>% 
      group_by(id, Stimulus) %>%
      dplyr::summarize(meanburn = mean(value, na.rm = T), seburn = sd(value,na.rm = T)/sqrt(n()))
    
    # want it in a wide format instead of long
    dfa <- responders %>% pivot_wider(names_from = Stimulus, values_from = c("seburn","meanburn"))
    
    # define the responders as a continus variable that is defined as the burning rating on TGI minus the average burning rating on the cold and warm stimulus
    dfa$max <- ifelse(dfa$meanburn_cold > dfa$meanburn_warm, dfa$meanburn_cold, dfa$meanburn_warm)
    dfa$color <- ifelse(dfa$meanburn_cold > dfa$meanburn_warm, "Cold", "Warm")
    
    #get the uncertainty on this:
    dfa$max_se <- ifelse(dfa$meanburn_cold > dfa$meanburn_warm, dfa$seburn_cold, dfa$seburn_warm)
    
    # Mean of the difference which is the responsivity index
    dfa$mean_continous_responsiveness <- dfa$meanburn_TGI - dfa$max
    
    # and the uncertainty  
    dfa$se_continous_responsiveness <- sqrt((dfa$seburn_TGI)^2 + (dfa$max_se)^2)
    
    dfa = dfa %>% mutate(responder = ifelse(mean_continous_responsiveness-2*se_continous_responsiveness > 0,1,0)) %>% mutate(responder = as.factor(responder))
    
    plot = dfa %>% ungroup() %>% arrange(mean_continous_responsiveness) %>% arrange() %>% mutate(ids = 1:nrow(.)) %>% ggplot()+
      geom_vline(xintercept = 0, size = axis_width, linetype = 2)+
      geom_pointrange(aes(y = ids, x = mean_continous_responsiveness,xmin = mean_continous_responsiveness-2*se_continous_responsiveness, xmax = mean_continous_responsiveness+2*se_continous_responsiveness),col = color_grey, alpha = 0.5, size = 0.1)+
      theme_classic()+
      scale_x_continuous("TGI responsiveness", breaks = pretty_breaks(n = 6), labels = comma) + 
      scale_y_continuous("Participant")+
      labs(color = "Maximal Burning")+
      # scale_color_manual(values = c(color_cold,color_warm))+
      theme+text+
      guides(color = guide_legend(override.aes = list(linetype = "blank", size = 1)))+
      theme(legend.position = c(0.8, 0.5),
            legend.box = "vertical",
            legend.direction="vertical",
            legend.justification='center',
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'),
            plot.title = element_text(hjust = 0.5), #change legend key height
            legend.title = element_text(size=font_size_small+2), #change legend title font size
            legend.text = element_text(size=font_size_small),
            axis.text=element_text(size=font_size_small),
            axis.title=element_text(size=font_size),
            axis.line=element_line(size=axis_width),
            axis.ticks=element_line(size=axis_width),
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    plot
    return(plot)
  }
  
  
  #getting the plots
  sa2_belief = belief_tgi(data)
  sa2_main = main_sa2(data)
  
  #new
  responder = responder_plot_v2(data)
  
  
  combined_sa2_belief_main <- (sa2_belief[[1]] + sa2_belief[[2]]) + 
    plot_layout(widths = c(2, 2, 3)) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = 'top')
  
  sa2_main = sa2_main+theme(legend.position = "none")
  
  
  #combing the plots
  plot4 = (combined_sa2_belief_main+plot_layout(widths = c(2,2))) /
    ((sa2_main+responder)+plot_layout(widths = c(3,3)))+
    plot_annotation(tag_levels = list(c("A","","B","C")))&
    patchtheme
  
  
  #save figure for manuscript markdown
  ggsave(filename = here::here("Figures", "figure4.png"), plot = plot4, width = 7.2, height = 7.2, units = "in", dpi = 600)
  
  #return plot for plot markdown
  return(plot4)
  
}





draw_table_main = function(){
  #read csv file from SPM (this table was generated by going into matlab running the results of the SPM.Mat from the MT folder and then selecting the pos_beta contrast)
  #here a threshold of p = 0.001 was set together with a cluseter threshold of 150. The resulting table was then extractet which is the csv file below:
  #Note that the mask inside the mask folder in VBQ directory was also used.
  #the regions found was found using the SPM anatomi toolbox on the same data.
  
  #function to generate tables for the qMRI results
  make_table = function(data,map,sign,contrast){
    
    data = data[,c(7,5,9,10,12:15)]
    names(data) = c("p (FWE)","k","TFCE", "Z value","x","y","z","Region")
    
    
    data = data %>%
      dplyr::select(Region,k,`p (FWE)`,TFCE,`Z value`,x,y,z) %>% 
      mutate_at(vars(-Region), ~as.numeric(.)) %>% 
      #drop_na() %>% 
      mutate_if(is.numeric, ~round(.,3))
    
    data = data %>% mutate(`p (FWE)` = ifelse(`p (FWE)` == 0 , "<.001",`p (FWE)`))
    
    data$MAP = map
    data$`Con.` = contrast
    
    if(sign == "neg"){
      data$`Z value` = -data$`Z value` 
    }
    
    return(data)
    
  }
  
  
  tables = list.files(here::here("matlab","VBQ","TFCE","tables"), full.names = T, recursive = T)
  
  
  table_df = data.frame()
  
  for(table in tables){
    
    # extract the MAP
    parts <- unlist(strsplit(table, "/"))
    MAP <- "R2*"
    
    # Extract "pos" or "neg" from the file name
    file_name <- parts[length(parts)]
    sign <- sub("_.*", "", basename(file_name))
    
    contrast <- sub(".*_(.*)\\.xls$", "\\1", file_name)  # Extract the contrast
    
    if(contrast == "omega"){
      contrast = "ω"
    }else if(contrast == "zeta"){
      contrast  = "ζ"
    }else if(contrast == "umti"){
      contrast = "UMTI"
    }else if(contrast == "resp"){
      contrast == "Resp"
    }
    
    
    table = make_table(readxl::read_xls(table, skip = 1),MAP,sign = 1,contrast)
    
    table_df = rbind(table_df, table)
  }
  
  
  
  
  table_df[,6:8] = round(table_df[,6:8])
  table_df[4] = round(table_df[4])
  
  table_df[5] = round(table_df[5],2)
  
  table_df$`p (FWE)` = gsub("^0", "", as.character(table_df$`p (FWE)`))
  
  table_df = table_df %>% mutate(TFCE = as.character(TFCE), k  = as.character(k))
  
  
  #keep first two rows
  table_df = table_df %>% mutate(k = if_else(is.na(k) & lag(!is.na(k)), lag(k), k)) %>% drop_na()
  
  table_df = table_df %>% mutate(`Z value` = as.character(`Z value`)) %>% group_by(k) %>%
    mutate(
      Region = ifelse(duplicated(Region), " ", Region),
      `p (FWE)` = ifelse(duplicated(k), " ", `p (FWE)`),
      `Z value` = ifelse(duplicated(k), " ", `Z value`),
      TFCE = ifelse(duplicated(k), " ", TFCE),
      k = ifelse(duplicated(k), " ", k)
    )
  
  collapse_rows <- function(df) {
    df %>%
      mutate(across(everything(), as.character)) %>%
      group_by(grp = cumsum(k != " ")) %>%
      summarise(across(everything(), ~ paste(na.omit(.), collapse = "\n\n")), .groups = 'drop') %>%
      select(-grp)
  }

  table_omega = table_df %>% drop_na() %>% filter(`Con.` == "ω") %>% mutate(`Con.` = NULL, MAP = NULL)
  table_zeta = table_df %>% drop_na() %>% filter(`Con.` == "ζ") %>% mutate(`Con.` = NULL, MAP = NULL)
  
  table_omega = collapse_rows(table_omega)
  table_zeta = collapse_rows(table_zeta)
  
  two_tables_small_omega = flextable::flextable(table_omega) %>% 
    theme_vanilla() %>% 
    bold(i = 1, j = NULL, bold = TRUE, part = "header") %>% 
    width(j = 1, width = 2.6) %>%
    width(j = 2, width = 0.75) %>%
    width(j = 3, width = 0.75) %>%
    width(j = 4, width = 0.6) %>%
    width(j = 5, width = 0.75) %>%
    width(j = 6:8, width = 0.45) %>%
    #align(i = 2, j = NULL, align = "center", part = "header") %>%
    fontsize(size = 10, part = "all")
  two_tables_small_omega
  
  
  two_tables_small_zeta = flextable::flextable(table_zeta) %>% 
    theme_vanilla() %>% 
    bold(i = 1, j = NULL, bold = TRUE, part = "header") %>% 
    width(j = 1, width = 2.6) %>%
    width(j = 2, width = 0.75) %>%
    width(j = 3, width = 0.75) %>%
    width(j = 4, width = 0.6) %>%
    width(j = 5, width = 0.75) %>%
    width(j = 6:8, width = 0.45) %>%
    #align(i = 2, j = NULL, align = "center", part = "header") %>%
    fontsize(size = 10, part = "all")
  two_tables_small_zeta
  

  
  return(list(two_tables_small_omega,two_tables_small_zeta))
}




draw_table_supplementary = function(){
  #read csv file from SPM (this table was generated by going into matlab running the results of the SPM.Mat from the MT folder and then selecting the pos_beta contrast)
  #here a threshold of p = 0.001 was set together with a cluseter threshold of 150. The resulting table was then extractet which is the csv file below:
  #Note that the mask inside the mask folder in VBQ directory was also used.
  #the regions found was found using the SPM anatomi toolbox on the same data.

  #function to generate tables for the qMRI results
  make_table = function(data,map,sign,contrast){
    
    data = data[,c(3,5,10,12:15)]
    names(data) = c("p(FWE-corr)","k", "Z value","x","y","zmm","Region")
    
    
    data = data %>%
      dplyr::select(Region,k,`p(FWE-corr)`,`Z value`,x,y,zmm) %>% 
      mutate_at(vars(-Region), ~as.numeric(.)) %>% 
      mutate_if(is.numeric, ~round(.,3)) %>% drop_na()
    
    data = data %>% mutate(`p(FWE-corr)` = ifelse(`p(FWE-corr)` == 0 , "<.001",`p(FWE-corr)`))
    
    data$MAP = map
    data$Contrast = contrast
    
    if(sign == "neg"){
      data$`Z value` = -data$`Z value` 
    }
    
    return(data)
    
  }
  
  
  tables = list.files(here::here("matlab","VBQ","tables"), full.names = T, recursive = T)
  
  
  table_df = data.frame()
  
  for(table in tables){
    
    # extract the MAP
    parts <- unlist(strsplit(table, "/"))
    MAP <- parts[length(parts) - 1]
    
    # Extract "pos" or "neg" from the file name
    file_name <- parts[length(parts)]
    sign <- sub("_.*", "", basename(file_name))
  
    contrast <- sub(".*_(.*)\\.xls$", "\\1", file_name)  # Extract the contrast
    
    if(contrast == "omega"){
      contrast = "ω"
    }else if(contrast == "zeta"){
      contrast  = "ζ"
    }else if(contrast == "umti"){
      contrast = "UMTI"
    }else if(contrast == "resp"){
      contrast == "Resp"
    }
    
    
    table = make_table(readxl::read_xls(table, skip = 1),MAP,sign,contrast)
    
    table_df = rbind(table_df, table)
  }
  
  
  
  
  table_df[,5:7] = round(table_df[,5:7])
  
  table_df[4] = round(table_df[4],2)
  
  ft = flextable::flextable(table_df)
  
  library(flextable)
  
  #modifying the flextable to look nice
  table1 = ft %>%
    theme_vanilla() %>% 
    bold(i = 1, j = NULL, bold = TRUE, part = "header") %>% 
    width(j = 1, width = 1.45) %>%
    width(j = 2, width = 0.5) %>%
    width(j = 3, width = 1.2) %>%
    width(j = 4, width = 0.75) %>%
    width(j = 5:7, width = 0.4) %>%
    width(j = 8, width = 0.5) %>%
    width(j = 9, width = 0.85) %>%
    align(i = NULL, j = 3:4, align = "center", part = "body") %>%
    align(i = NULL, j = 9, align = "center", part = "body") %>%
    #align(i = 2, j = NULL, align = "center", part = "header") %>%
    fontsize(size = 10, part = "all")
  
  exportxlsx(table1, path = here::here("matlab","MPM cluster-based results","cluster-based-inference.xlsx"))
  
  
  return(table1)
}


