function plot_zeta(u1, sim)

subplot(3,3,1)
dn1 = pp_check(u1,sim,tapas_hgf_binary_pu_tgi_config_test,tapas_unitsq_sgm_config1, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(dn1),3)
s1 = round(std(dn1),3)
title("1      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,2)
dn2 = pp_check(u1,sim,tapas_hgf_binary_pu_tgi_config_test,tapas_unitsq_sgm_config2, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(dn2),3)
s1 = round(std(dn2),3)
title("2      "+"mean = "+ m1 + "              sd ="+ s1)

subplot(3,3,3)
dn3 = pp_check(u1,sim,tapas_hgf_binary_pu_tgi_config_test,tapas_unitsq_sgm_config3, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(dn3),3)
s1 = round(std(dn3),3)
title("3      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,4)
dn4 = pp_check(u1,sim,tapas_hgf_binary_pu_tgi_config_test,tapas_unitsq_sgm_config4, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(dn4),3)
s1 = round(std(dn4),3)
title("4      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,5)
dn5 = pp_check(u1,sim,tapas_hgf_binary_pu_tgi_config_test,tapas_unitsq_sgm_config5, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(dn5),3)
s1 = round(std(dn5),3)
title("5      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,6)
dn10 = pp_check(u1,sim,tapas_hgf_binary_pu_tgi_config_test,tapas_unitsq_sgm_config10, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(dn10),3)
s1 = round(std(dn10),3)
title("10      "+"mean = "+ m1 + "              sd ="+ s1)


subplot(3,3,7)
dn50 = pp_check(u1,sim,tapas_hgf_binary_pu_tgi_config_test,tapas_unitsq_sgm_config50, tapas_hgf_binary_pu_tgi_config_test_optimal)
m1 = round(mean(dn50),3)
s1 = round(std(dn50),3)
title("50      "+"mean = "+ m1 + "              sd ="+ s1)
