function plot_means_w2(u1, sim)
subplot(3,3,1)
omegawm1 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wm1, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm1),3)
s1 = round(std(omegawm1),3)
title("-1      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,2)
omegawm2 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wm2, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm2),3)
s1 = round(std(omegawm2),3)
title("-2      "+"mean = "+ m1 + "              sd ="+ s1)

subplot(3,3,3)
omegawm3 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wm3, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm3),3)
s1 = round(std(omegawm3),3)
title("-3      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,4)
omegawm4 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wm4, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm4),3)
s1 = round(std(omegawm4),3)
title("-4      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,5)
omegawm5 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wm5, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm5),3)
s1 = round(std(omegawm5),3)
title("-5      "+"mean = "+ m1 + "              sd ="+ s1)

subplot(3,3,6)

omegawm6 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wm6, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm6),3)
s1 = round(std(omegawm6),3)
title("-6      "+"mean = "+ m1 + "              sd ="+ s1)

subplot(3,3,7)

omegawm7 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wm7, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm7),3)
s1 = round(std(omegawm7),3)
title("-7      "+"mean = "+ m1 + "              sd ="+ s1)

subplot(3,3,8)
omegawm8 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wm8, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm8),3)
s1 = round(std(omegawm8),3)
title("-8      "+"mean = "+ m1 + "              sd ="+ s1)
