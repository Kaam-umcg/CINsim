library(CINsim)

sim_res <- Cinsim(pMisseg = 0.0025, g = 25, selection_mode = NULL, collect_fitness_score = TRUE)

cnvHeatmap(sim_res)

plot_cn(sim_res)

plot_cn(sim_res, final_g = FALSE)

plot_misseg_freq(sim_res)

#plot_cinsim_summary(sim_res)

start_pop <- makeKaryotypes(n = 10)

sim_res <- Cinsim(karyotypes = start_pop, pMisseg = 0.0025, g = 25, selection_mode = NULL)

plot_clonality(sim_res)

