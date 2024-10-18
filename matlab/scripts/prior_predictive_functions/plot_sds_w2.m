function plot_means_w2(u1, sim)
subplot(3,3,1)
omegawm1 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wsd1, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm1),3)
s1 = round(std(omegawm1),3)
title("1      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,2)
omegawm2 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wsd2, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm2),3)
s1 = round(std(omegawm2),3)
title("2      "+"mean = "+ m1 + "              sd ="+ s1)

subplot(3,3,3)
omegawm3 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wsd3, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm3),3)
s1 = round(std(omegawm3),3)
title("3      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,4)
omegawm4 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wsd4, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm4),3)
s1 = round(std(omegawm4),3)
title("4      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,5)
omegawm5 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wsd5, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm5),3)
s1 = round(std(omegawm5),3)
title("5      "+"mean = "+ m1 + "              sd ="+ s1)

subplot(3,3,6)

omegawm6 = pp_check_obs(u1,sim,tapas_hgf_binary_pu_tgi_config_test_demo_wsd6, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(omegawm6),3)
s1 = round(std(omegawm6),3)
title("6      "+"mean = "+ m1 + "              sd ="+ s1)
