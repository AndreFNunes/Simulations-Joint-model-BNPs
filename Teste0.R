##################################################################################
################
#                TESTE 0 REFERÊNCIA SEM BNP's
################
##################################################################################

simCode <- nimbleCode({
  for (i in 1:nsite) {
    for (j in 1:nvisit) {
      logit(prob[i,j]) <- beta0 + inprod(beta[1:p], x[i, 1:p]) + b[i] # True occupancy status
      logit(prob2[i,j]) <- beta02 + inprod(beta2[1:p], x[i, 1:p]) + b[i] # True occupancy status
      
      y1[i,j] ~ dbern(prob[i,j]) # Observed data
      y2[i,j] ~ dbern(prob2[i,j]) # Observed data
    }
    
    b[i] ~ dnorm(0, var = tau2_b)
  }
  
  tau2_b ~ dgamma(0.1,0.2)
  
  beta0 ~ dnorm(0, sd = 50)
  beta02 ~ dnorm(0, sd = 50)
  for(k in 1:p) {
    beta[k] ~ dnorm(0, sd = 50)
    beta2[k] ~ dnorm(0, sd = 50)
  }
  
})

p <- 3

simModel <- nimbleModel(code = simCode,
                        constants = list(nsite = 120, nvisit = 21, p = p, 
                                         x = X), 
                        inits = list(beta0 = 1, beta = rnorm(p, 0, 1),
                                     beta02 = -1, beta2 = rnorm(p, 0, 1),
                                     b = rnorm(120, 0, 2), tau2_b = 1)
)

CsimModel <- compileNimble(simModel)

CsimModel$setData(list(y1 = simulatedY1s))

CsimModel$setData(list(y2 = simulatedY2s))


gamma_list <- sprintf("b[%d]", 1:120)
prob_list <- c()
prob2_list <- c()
for (i in 1:21) {
  prob_list <- c(prob_list, sprintf(sprintf("prob[%%d, %d]", 1:21)[i], 1:120))
  prob2_list <- c(prob2_list, sprintf(sprintf("prob2[%%d, %d]", 1:21)[i], 1:120)) 
}

simMCMC <- buildMCMC(CsimModel, monitors = c(c("beta0", "beta02", 
                                               "beta[1]", "beta[2]", "beta[3]",
                                               "beta2[1]", "beta2[2]", "beta2[3]", "tau2_b"),
                                             gamma_list, prob_list, prob2_list))

simMCMC$monitors

CsimMCMC <- compileNimble(simMCMC, project = simModel)



samples <- runMCMC(CsimMCMC, niter = 13000, nburnin = 5000, thin = 50, summary = TRUE)

par(mfrow=c(1,2))
# 
# plot(samples$samples[ , 'beta0'], type = 'l', xlab = 'iteration',  ylab = expression(beta0))
# plot(samples$samples[ , 'beta[1]'], type = 'l', xlab = 'iteration', ylab = expression(beta1))
# plot(samples$samples[ , 'beta[2]'], type = 'l', xlab = 'iteration', ylab = expression(beta2))
# plot(samples$samples[ , 'beta[3]'], type = 'l', xlab = 'iteration', ylab = expression(beta3))


##########
# 
# plot(samples$samples[ , 'beta02'], type = 'l', xlab = 'iteration',  ylab = expression(beta02))
# plot(samples$samples[ , 'beta2[1]'], type = 'l', xlab = 'iteration', ylab = expression(beta21))
# plot(samples$samples[ , 'beta2[2]'], type = 'l', xlab = 'iteration', ylab = expression(beta22))
# plot(samples$samples[ , 'beta2[3]'], type = 'l', xlab = 'iteration', ylab = expression(beta23))


mu_gammaCols_23_chain1 <- samples$samples[, gamma_list]
mu_gammaMn_23_chain1 <- colMeans(mu_gammaCols_23_chain1)
hist(mu_gammaMn_23_chain1[1:30], xlab = 'posterior means 23', main = 'mu_gamma', probability = FALSE, breaks = 20, col='blue', xlim = c(-4,4), ylim = c(0,15), density = 5, angle = 0, border = 'blue')
hist(mu_gammaMn_23_chain1[31:60], col='red', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 20, density = 5, angle = 30, border = 'red')
hist(mu_gammaMn_23_chain1[61:90], col='orange', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 20, density = 5, angle = 60, border = 'orange')
hist(mu_gammaMn_23_chain1[91:120], col='green', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 20, density = 5, angle = 80, border = 'green')


plot(samples$samples[ , 'tau2_b'], type = 'l', xlab = 'iteration', ylab = expression(tau2_b))









prob_list_estim <- samples$samples[, prob_list]
prob_list_estim_mean <- matrix(colMeans(prob_list_estim), nrow = 120, ncol = 21)

prob2_list_estim <- samples$samples[, prob2_list]
prob2_list_estim_mean <- matrix(colMeans(prob2_list_estim), nrow = 120, ncol = 21)

typeI1 <- rep(0,21)
typeII1 <- rep(0,21)

true_positive1 <- rep(0,21)
true_negative1 <- rep(0,21) 

typeI2 <- rep(0,21)
typeII2 <- rep(0,21)

true_positive2 <- rep(0,21)
true_negative2 <- rep(0,21) 

for (i in 1:120) {
  for (j in 1:21) {
    
    if (prob_list_estim_mean[i,j] <= 0.5) {
      if (simulatedY1s[i,j] == 0) {
        true_negative1[j] <- true_negative1[j] + 1
      } else {typeII1[j] <- typeII1[j] + 1}
    }
    if (prob_list_estim_mean[i,j] >= 0.5) {
      if (simulatedY1s[i,j] == 1) {
        true_positive1[j] <- true_positive1[j] + 1
      } else {typeI1[j] <- typeI1[j] + 1}
    }
    
    if (prob2_list_estim_mean[i,j] <= 0.5) {
      if (simulatedY2s[i,j] == 0) {
        true_negative2[j] <- true_negative2[j] + 1
      } else {typeII2[j] <- typeII2[j] + 1}
    }
    if (prob2_list_estim_mean[i,j] >= 0.5) {
      if (simulatedY2s[i,j] == 1) {
        true_positive2[j] <- true_positive2[j] + 1
      } else {typeI2[j] <- typeI2[j] + 1}
    }
  }
}

typeI1 
typeII1 
true_positive1
true_negative1 
typeI2 
typeII2 
true_positive2
true_negative2

sample_precision1 <- true_positive1/(typeII1 + true_positive1)
sample_consistency1 <- true_negative1/(typeI1 + true_negative1)

sample_precision2 <- true_positive2/(typeII2 + true_positive2)
sample_consistency2 <- true_negative2/(typeI2 + true_negative2)

sample_precision1
sample_consistency1
sample_precision2
sample_consistency2

summary(sample_precision1)
summary(sample_consistency1)
summary(sample_precision2)
summary(sample_consistency2)




limit <- (1:100)/100
samples_precision1_ROC <- rep(0, 100)
samples_consistency1_ROC <- rep(0, 100)

samples_precision2_ROC <- rep(0, 100)
samples_consistency2_ROC <- rep(0, 100)

for (k in 1:100) {
  
  typeI1 <- rep(0,21)
  typeII1 <- rep(0,21)
  
  true_positive1 <- rep(0,21)
  true_negative1 <- rep(0,21) 
  
  typeI2 <- rep(0,21)
  typeII2 <- rep(0,21)
  
  true_positive2 <- rep(0,21)
  true_negative2 <- rep(0,21) 
  
  for (i in 1:120) {
    for (j in 1:21) {
      
      if (prob_list_estim_mean[i,j] <= limit[k]) {
        if (simulatedY1s[i,j] == 0) {
          true_negative1[j] <- true_negative1[j] + 1
        } else {typeII1[j] <- typeII1[j] + 1}
      }
      if (prob_list_estim_mean[i,j] >= limit[k]) {
        if (simulatedY1s[i,j] == 1) {
          true_positive1[j] <- true_positive1[j] + 1
        } else {typeI1[j] <- typeI1[j] + 1}
      }
      
      if (prob2_list_estim_mean[i,j] <= limit[k]) {
        if (simulatedY2s[i,j] == 0) {
          true_negative2[j] <- true_negative2[j] + 1
        } else {typeII2[j] <- typeII2[j] + 1}
      }
      if (prob2_list_estim_mean[i,j] >= limit[k]) {
        if (simulatedY2s[i,j] == 1) {
          true_positive2[j] <- true_positive2[j] + 1
        } else {typeI2[j] <- typeI2[j] + 1}
      }
    }
  }
  
  sum_typeI1 <- sum(typeI1)
  sum_true_positive1 <- sum(true_positive1)
  sum_typeII1 <- sum(typeII1)
  sum_true_negative1 <- sum(true_negative1)
  
  sum_typeI2 <- sum(typeI2)
  sum_true_positive2 <- sum(true_positive2)
  sum_typeII2 <- sum(typeII2)
  sum_true_negative2 <- sum(true_negative2)
  
  samples_precision1_ROC[k] <- sum_true_positive1/(sum_typeII1 + sum_true_positive1)
  samples_consistency1_ROC[k] <- sum_typeI1/(sum_typeI1 + sum_true_negative1)
  
  samples_precision2_ROC[k] <- sum_true_positive2/(sum_typeII2 + sum_true_positive2)
  samples_consistency2_ROC[k] <- sum_typeI2/(sum_typeI2 + sum_true_negative2)
}

samples_precision1_ROC
samples_consistency1_ROC

samples_precision2_ROC
samples_consistency2_ROC