#------------------------- utils

## generate compte data
gencox <- function(n = 700, maxTime = 7) {
  r <- 0.4 ; a <- 0.2 ; k <- 2 ;  rateC = 0.09 # ~ 30%
  ## x1 normal distribution
  x1 = rnorm(n)
  ## x2 standard normal derived from x1
  x2 = x1^2 + x1 + runif(n)
  ## x3  binary variable correlated to x2
  x3 <- rbinom(n, 1, 0.5) # binary variable
  
  ## latent event times
  U <- runif(n)
  
  ##  survival times
  exp.betaZ <- exp(0.1*x1 + 0.3*x2 +  0.6*x3)
  Tlat <- (1/r*(a*((1/(1-U))^(1/(a*exp.betaZ))-1))^(1/k))
  
  ## censoring times
  Ctimes <- rexp(n, rate = rateC)
  
  ## observed time is min of censored and true
  time <- pmin(Tlat, Ctimes)
  status <- as.numeric(Tlat <= Ctimes)
  
  # administrative censoring
  time <- ifelse(time > maxTime, maxTime, time)
  status <- ifelse(time > maxTime | status == 1, 1, 0)
  
  data <- data.frame(time, status, x1, x2, x3)
  x3 <- as.factor(x3)
  data <- data.frame(time, status, x1, x2, x3)
  
  ## Observed marginal cumulative hazard (Nelson-Altschuler-Aalen estimator)
  data$cumhaz <- mice::nelsonaalen(data, time, status)
  return(data)
  }


## generate NA under MAR (inspired from Shah et al., 2014a)
genNA <- function(data, pmiss){
  logistic <- function(x){
    exp(x) / (1 + exp(x))
  }
  predictions <- function(lp, n){
    trialn <- function(lptrial){
      sum(logistic(lptrial))
    }
    stepsize <- 32
    lptrial <- lp
    if (any(is.na(lptrial))){
      lp[is.na(lptrial)] <- mean(lptrial, na.rm = TRUE)
    }
    while(abs(trialn(lptrial) - n) > 1){
      if (trialn(lptrial) > n){
        lptrial <- lptrial - stepsize
      } else {
        lptrial <- lptrial + stepsize
      }
      stepsize <- stepsize / 2
    }
    as.logical(rbinom(logical(length(lp)), 1, logistic(lptrial)))
  }
  data$x2[predictions(0.1 * data$x1  + 0.1 * data$status + 0.1 * data$time, nrow(data) * pmiss/1.8)] <- NA
  data$x3[predictions(0.1 * data$x1  + 0.1 * data$status + 0.1 * data$time, nrow(data) * pmiss/1.8)] <- NA
  
  return(data)
  }


## evaluate post-imputation estimates accuracy 

### for single imputation methods
evaluate_coxest_single <- function(data, truelogHR){
  myformula <- as.formula(survival::Surv(time, status) ~ x1 + x2 + x3)
  dat <- data
  coefs <- as.data.frame(summary(survival::coxph(myformula, data = dat))$coef)

  out <- data.frame(covariates = row.names(coefs),
                    est = coefs$coef,
                    se_est = coefs$`se(coef)`,
                    lo95 = (coefs$coef + qnorm(0.025) * coefs$`se(coef)`),
                    hi95 = (coefs$coef + qnorm(0.975) * coefs$`se(coef)`)
                    )
  out$bias = coefs$coef - truelogHR
  out$bias_rel =  ((coefs$coef - truelogHR) / truelogHR) * 100
  out$MSE_estimate = (coefs$coef - truelogHR)^2
  out$cover <- truelogHR >= out$lo95 & truelogHR <= out$hi95
  data.frame(out)
  }

### for multiple imputation methods
evaluate_coxest_multiple <- function(mids, truelogHR){
  # pooled results
  fit <- with(mids, coxph(Surv(time, status) ~ x1 + x2 + x3))
  res <- summary(pool(fit))
  out <- data.frame(covariates = as.character(res$term),
                    est = res$estimate,
                    se_est = res$std.error,
                    lo95 = (res$estimate + qnorm(0.025) * res$std.error),
                    hi95 = (res$estimate + qnorm(0.975) * res$std.error)
                    )
  out$bias = out$est - truelogHR
  out$bias_rel =  ((out$est - truelogHR) / truelogHR)*100
  out$MSE_estimate = (out$est - truelogHR)^2
  out$cover <- truelogHR >= out$lo95 & truelogHR <= out$hi95
  data.frame(out)
}


## evaluate post-imputation predictive accuracy 

### AUC 9Chambless, L. E. and G. Diao, 2006)
coxperf_AUC <- function(data){
  myformula <- as.formula(survival::Surv(time, status) ~ x1 + x2 + x3)
  N = nrow(data)
  #Initialization
  index.train = sample(1:N,2/3*N)
  data.train = data[index.train,]
  data.test = data[-index.train,]
  dis_time = sort(data.train$time[data.train$status == 1])  #the default time points
  fitcox = survival::coxph(myformula, data = data.train, x = TRUE)
  predcox <- predict(fitcox)
  predcoxnew <- predict(fitcox, newdata=data.test)
  surv_obj <- survival::Surv(data.train$time, data.train$status)
  surv_obj_new <- survival::Surv(data.test$time, data.test$status)
  AUC = survAUC::AUC.cd(surv_obj, surv_obj_new, predcox, predcoxnew, dis_time)$iauc
  data.frame(AUC)
}

### Cindex
coxperf_Cindex <- function(data){
  (survival::coxph(Surv(time, status) ~ x1 + x2 + x3, data))$concordance[6]
  }

## substantive model-free performance metrics (adapted only for multiple imputation, see simcox_sim() for single imputation)

NRMSE_mutiple <- function(impdat_list, amputed_data, complete_data){
  NRMSE <- function(impdat_list){
    missForest::mixError(ximp = impdat_list, xmis = amputed_data, xtrue = complete_data)[1]
    }
  purrr::map(impdat_list, NRMSE)
  }

PFC_mutiple <- function(impdat_list, amputed_data, complete_data){
  PFC <- function(impdat_list){
    missForest::mixError(ximp = impdat_list, xmis = amputed_data, xtrue = complete_data)[2]
    }
  purrr::map(impdat_list, PFC)
  }

