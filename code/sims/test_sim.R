args$sim_name <- "bivariate_loss"
args$vimp_measure <- c("deviance", "accuracy", "auc")
args$cv <- 0
job_id <- 1
set.seed(current_seed)
system.time(sim_output <- replicate(5,
                                    run_sim_binary_bivariate_once(n = current_dynamic_args$n,
                                                                  j = current_dynamic_args$j,
                                                                  p = p,
                                                                  mu_0 = mu_0, mu_1 = mu_1, sigma = Sigma,
                                                                  truth = truth[current_dynamic_args$j],
                                                                  b = args$b, V = V,  
                                                                  risk_type = args$risk_type,
                                                                  learner_lib = learner_lib,
                                                                  type = args$vimp_measure,
                                                                  cv = args$cv),
                                    simplify = FALSE)
           )
sim_output <- lapply(as.list(1:length(sim_output)), function(x) tibble::add_column(sim_output[[x]], mc_id = x))
sim_output_tib <- do.call(rbind.data.frame, sim_output)
sim_output_tib

auc_cover <- sim_output_tib %>% 
  mutate(cover = cil <= truth & ciu >= truth) %>% 
  select(-n, -j, -se) %>% 
  filter(type == "auc")
auc_cover
auc_cover %>% 
  summarize(cover = mean(cover))
