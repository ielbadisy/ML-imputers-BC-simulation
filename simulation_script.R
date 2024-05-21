rm(list=ls())
source("utils.R")
#install.packages("devtools")
#devtools::install_github("ielbadisy/imputer")
pacman::p_load(dplyr, ggplot2, gower, miceadds, stringr, imputer, parallel, fastDummies, mice, survival, survAUC, purrr, pbapply)
set.seed = 1234
R = 500
ncores = 32

#---- Single imoutation methods
## grid of parameters 
grid_par_single <- expand.grid(N = 700,
                               pmiss = c(0.3),
                               method = c("complete", "knn", "cart", "famd", "missranger", "missforest", "misscforest"))

## simulation function
simcox_sim <- function(N, pmiss, method = method){
  dat <- gencox(n = N)
  truelogHR <- c(0.1, 0.3, 0.6)
  ampdat <- genNA(data= dat, pmiss = pmiss)
  suppressWarnings({impdat <- imputer::imputer(ampdat, method = method)})
  res_est <- evaluate_coxest_single(impdat, truelogHR)
  AUC <- coxperf_AUC(impdat)
  Cindex <- coxperf_Cindex(impdat)
  Gowerdist <- mean(gower_dist(dat, impdat))
  NRMSE <- missForest::mixError(ximp = impdat, xmis = ampdat, xtrue = dat)[1]
  PFC <- missForest::mixError(ximp = impdat, xmis = ampdat, xtrue = dat)[2]
  res <- cbind(res_est, AUC, Cindex, Gowerdist, NRMSE, PFC, method, pmiss, row.names = NULL)
  return(res)
  }

## simulation routine
simcox_sim2 <- purrr::possibly(.f = simcox_sim, otherwise = NULL)
out_single <- pbapply::pbreplicate(R, purrr::pmap(grid_par_single, simcox_sim2), simplify = FALSE, cl = ncores)
res_single <- dplyr::bind_rows(out_single, .id = "ID")

## single imputation results
write.csv(res_single, "res_single5.csv")

#---- Multiple imputation methods
## grid of parameters 
grid_par_multiple <- expand.grid(N = 700,
                                 pmiss = c(0.3),
                                 method = c("micecart", "micerf"))

## simulation function
simcox_sim_multiple <- function(N, pmiss, method = method){
  dat <- gencox(n = N)
  truelogHR <- c(0.1, 0.3, 0.6)
  ampdat <- genNA(data= dat, pmiss = pmiss)
  impdat_mids <- imputer::imputer(ampdat, method = method)
  impdat_list <- miceadds::mids2datlist(impdat_mids)
  dat_mice <- mice::complete(impdat_mids, "all")
  res_est <- evaluate_coxest_multiple(impdat_mids, truelogHR)
  AUC <- purrr::map(impdat_list, coxperf_AUC) %>% dplyr::bind_rows() %>% unlist() %>% as.numeric() %>% mean()
  Cindex <- purrr::map(impdat_list, coxperf_Cindex) %>% dplyr::bind_rows() %>%unlist() %>% as.numeric() %>% mean()
  Gowerdist <-  lapply(dat_mice, function(x) mean(gower_dist(x, dat))) %>% as.numeric() %>% mean()
  NRMSE <- NRMSE_mutiple(impdat_list, ampdat, dat) %>% dplyr::bind_rows() %>% colMeans()
  PFC <- PFC_mutiple(impdat_list, ampdat, dat) %>% dplyr::bind_rows() %>% colMeans()
  res <- cbind(res_est, AUC, Cindex, Gowerdist, NRMSE, PFC, method, pmiss, row.names = NULL)
  return(res)
  }

## multiple imputation results
simcox_sim_multiple2 <- purrr::possibly(.f = suppressWarnings({simcox_sim_multiple}), otherwise = NULL)
out_multiple <- pbapply::pbreplicate(R, purrr::pmap(grid_par_multiple, simcox_sim_multiple2), simplify = FALSE, cl = ncores)
res_multiple <- dplyr::bind_rows(out_multiple, .id = "ID")

write.csv(res_multiple, "res_multiple5.csv")

#-------- bind results
res <- rbind(res_single, res_multiple)
write.csv(res, "res_simulation.csv")

