function plot_RW_mean(u1, sim)
subplot(3,3,1)
rw01 = pp_check_obs_rw(u1,sim,tapas_rw_binary_configsd01, tapas_rw_binary_optimal_config)
m1 = round(mean(rw01),3)
s1 = round(std(rw01),3)
title("0.1      "+"mean = "+ m1 + "              sd ="+ s1)
subplot(3,3,2)
rw05 = pp_check_obs_rw(u1,sim,tapas_rw_binary_configsd02, tapas_rw_binary_optimal_config)
m1 = round(mean(rw05),3)
s1 = round(std(rw05),3)
title("0.2      "+"mean = "+ m1 + "              sd ="+ s1)
subplot(3,3,3)
rw07 = pp_check_obs_rw(u1,sim,tapas_rw_binary_configsd03, tapas_rw_binary_optimal_config)
m1 = round(mean(rw07),3)
s1 = round(std(rw07),3)
title("0.3      "+"mean = "+ m1 + "              sd ="+ s1)
subplot(3,3,4)
rw1 = pp_check_obs_rw(u1,sim,tapas_rw_binary_configsd04, tapas_rw_binary_optimal_config)
m1 = round(mean(rw1),3)
s1 = round(std(rw1),3)
title("0.4      "+"mean = "+ m1 + "              sd ="+ s1)
subplot(3,3,5)

rw1 = pp_check_obs_rw(u1,sim,tapas_rw_binary_configsd05, tapas_rw_binary_optimal_config)
m1 = round(mean(rw1),3)
s1 = round(std(rw1),3)
title("0.5      "+"mean = "+ m1 + "              sd ="+ s1)



subplot(3,3,6)

rw1 = pp_check_obs_rw(u1,sim,tapas_rw_binary_configsd1, tapas_rw_binary_optimal_config)
m1 = round(mean(rw1),3)
s1 = round(std(rw1),3)
title("1      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,7)

rw1 = pp_check_obs_rw(u1,sim,tapas_rw_binary_configsd2, tapas_rw_binary_optimal_config)
m1 = round(mean(rw1),3)
s1 = round(std(rw1),3)
title("2      "+"mean = "+ m1 + "              sd ="+ s1)

