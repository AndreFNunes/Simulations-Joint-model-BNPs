
simCode <- nimbleCode({
  for (i in 1:nsite) {
    gamma[i] ~ dnorm(mu[i], var = tau2[i])  # random effects from mixture dist.
    mu[i] <- muTilde[xi[i]]                 # mean for random effect from cluster xi[i]
    tau2[i] <- tau2Tilde[xi[i]]             # var for random effect from cluster xi[i] 
    xi[i] ~ dcat(w[1:L])
    
    for (j in 1:nvisit) {
      logit(prob[i,j]) <- beta0 + inprod(beta[1:p], x[i, 1:p]) + gamma[i]# True occupancy status
      logit(prob2[i,j]) <- beta02 + inprod(beta2[1:p], x[i, 1:p]) + gamma[i]# True occupancy status
      
      y1[i,j] ~ dbern(prob[i,j]) # Observed data
      y2[i,j] ~ dbern(prob2[i,j]) # Observed data
    }
  }
  
  beta0 ~ dnorm(0, sd = 1)
  beta02 ~ dnorm(0, sd = 1)
  for(k in 1:p) {
    beta[k] ~ dnorm(0, sd = 2)
    beta2[k] ~ dnorm(0, sd = 2)
  }
  # mixture component parameters drawn from base measures
  for(i in 1:L) {
    muTilde[i] ~ dnorm(mu0, var = var0)
    tau2Tilde[i] ~ dinvgamma(a0, b0)
  }
  
  for(i in 1:(L-1)){
    v[i] ~ dbeta(1, alpha)
  }
  
  alpha ~ dgamma(1, 1)
  
  w[1:L] <- stick_breaking(v[1:(L-1)])
})

p <- 3
L <- 30

simModel <- nimbleModel(code = simCode,
                        constants = list(nsite = 120, nvisit = 21, p = p, L = L,
                                         a0 = 2, b0 = 0.2, mu0 = 0, var0 = 25,
                                         x = X), 
                        inits = list(beta0 = 1, beta = rnorm(p, 0, 1),
                                     beta02 = -1, beta2 = rnorm(p, 0, 1),
                                     alpha = 1, muTilde = rnorm(L, 0, 4), 
                                     tau2Tilde = rinvgamma(L, 10, 10),
                                     xi = rep(seq(1,30),4),
                                     v = rbeta(L-1, 1, 1),
                                     gamma = rnorm(120, 0, 4))
)




CsimModel <- compileNimble(simModel)

CsimModel$setData(list(y1 = simulatedY1s))

CsimModel$setData(list(y2 = simulatedY2s))




xi_list <- sprintf("xi[%d]", 1:120)
mu_list <- sprintf("mu[%d]", 1:120)
w_list <- sprintf("w[%d]", 1:L)

tau2_list <- sprintf("tau2[%d]", 1:120)
gamma_list <- sprintf("gamma[%d]", 1:120)

prob_list <- c()
prob2_list <- c()
for (i in 1:21) {
  prob_list <- c(prob_list, sprintf(sprintf("prob[%%d, %d]", 1:21)[i], 1:120))
  prob2_list <- c(prob2_list, sprintf(sprintf("prob2[%%d, %d]", 1:21)[i], 1:120)) 
}

simMCMC <- buildMCMC(CsimModel, monitors = c(c("beta0", "beta02", 
                                               "beta[1]", "beta[2]", "beta[3]",
                                               "beta2[1]", "beta2[2]", "beta2[3]",
                                               "alpha"),
                                             xi_list, mu_list, w_list, tau2_list, gamma_list, prob_list, prob2_list))


CsimMCMC <- compileNimble(simMCMC, project = simModel)



samples1 <- runMCMC(CsimMCMC, niter = 13000, nburnin = 5000, thin = 50, summary = TRUE)


par(mfrow=c(1,2))
# 
# plot(samples1$samples[ , 'beta0'], type = 'l', xlab = 'iteration',  ylab = expression(beta0))
# plot(samples1$samples[ , 'beta[1]'], type = 'l', xlab = 'iteration', ylab = expression(beta1))
# plot(samples1$samples[ , 'beta[2]'], type = 'l', xlab = 'iteration', ylab = expression(beta2))
# plot(samples1$samples[ , 'beta[3]'], type = 'l', xlab = 'iteration', ylab = expression(beta3))
# 
# ##########
# 
# plot(samples1$samples[ , 'beta02'], type = 'l', xlab = 'iteration',  ylab = expression(beta02))
# plot(samples1$samples[ , 'beta2[1]'], type = 'l', xlab = 'iteration', ylab = expression(beta21))
# plot(samples1$samples[ , 'beta2[2]'], type = 'l', xlab = 'iteration', ylab = expression(beta22))
# plot(samples1$samples[ , 'beta2[3]'], type = 'l', xlab = 'iteration', ylab = expression(beta23))

n <- 120
########################## GAMMA
{
  gamma_list <- c()
  mu_gamma_list <- c()
  tau2_gamma_list <- c()
  for (i in 1:n) {
    eq <- paste('gamma_list <- c(gamma_list, "gamma[',i,']")', sep = "")
    eval(parse(text = eq))
    eq <- paste('mu_gamma_list <- c(mu_gamma_list, "mu[',i,']")', sep = "")
    eval(parse(text = eq))
    eq <- paste('tau2_gamma_list <- c(tau2_gamma_list, "tau2[',i,']")', sep = "")
    eval(parse(text = eq))
  }
}

#par(mfcol=c(1, 1))

###############################################################################################
{{# gamma
  # gammaCols_23_chain1 <- samples1$samples[, gamma_list]
  # gammaMn_23_chain1 <- colMeans(gammaCols_23_chain1)
  # # gammaCols_23_chain2 <- data_samples_23$samples[, gamma_list]
  # # gammaMn_23_chain2 <- colMeans(gammaCols_23_chain2)
  # # gammaCols_23_chain3 <- data_samples_23$samples[, gamma_list]
  # # gammaMn_23_chain3 <- colMeans(gammaCols_23_chain3)
  # # gammaCols_23_chain4 <- data_samples_23$samples[, gamma_list]
  # # gammaMn_23_chain4 <- colMeans(gammaCols_23_chain4)
  # hist(gammaMn_23_chain1[1:30], xlab = 'posterior means 23', main = 'gamma', probability = FALSE, breaks = 100, col='blue', xlim = c(-6,6), ylim = c(0,30), density = 5, angle = 0, border = 'blue') 
  # hist(gammaMn_23_chain1[31:60], col='red', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 100, density = 5, angle = 30, border = 'red')
  # hist(gammaMn_23_chain1[61:90], col='orange', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 100, density = 5, angle = 60, border = 'orange')
  # hist(gammaMn_23_chain1[91:120], col='green', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 100, density = 5, angle = 60, border = 'green')
  # 
}##############################################################################################
  {# mu_gamma
    mu_gammaCols_23_chain1 <- samples1$samples[, mu_gamma_list]
    mu_gammaMn_23_chain1 <- colMeans(mu_gammaCols_23_chain1)
    # mu_gammaCols_23_chain2 <- data_samples_23$samples[, mu_gamma_list]
    # mu_gammaMn_23_chain2 <- colMeans(mu_gammaCols_23_chain2)
    # mu_gammaCols_23_chain3 <- data_samples_23$samples[, mu_gamma_list]
    # mu_gammaMn_23_chain3 <- colMeans(mu_gammaCols_23_chain3)
    # mu_gammaCols_23_chain4 <- data_samples_23$samples[, mu_gamma_list]
    # mu_gammaMn_23_chain4 <- colMeans(mu_gammaCols_23_chain4)
    hist(mu_gammaMn_23_chain1[1:30], xlab = 'posterior means 23', main = 'mu_gamma', probability = FALSE, breaks = 20, col='blue', xlim = c(-4,4), ylim = c(0,30), density = 5, angle = 0, border = 'blue')
    hist(mu_gammaMn_23_chain1[31:60], col='red', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 20, density = 5, angle = 30, border = 'red')
    hist(mu_gammaMn_23_chain1[61:90], col='orange', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 20, density = 5, angle = 60, border = 'orange')
    hist(mu_gammaMn_23_chain1[91:120], col='green', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 20, density = 5, angle = 80, border = 'green')
    
  }##############################################################################################
  {# tau2_gamma
    tau2_gammaCols_23_chain1 <- samples1$samples[, tau2_gamma_list]
    tau2_gammaMn_23_chain1 <- colMeans(tau2_gammaCols_23_chain1)
    # tau2_gammaCols_23_chain2 <- data_samples_23$samples[, tau2_gamma_list]
    # tau2_gammaMn_23_chain2 <- colMeans(tau2_gammaCols_23_chain2)
    # tau2_gammaCols_23_chain3 <- data_samples_23$samples[, tau2_gamma_list]
    # tau2_gammaMn_23_chain3 <- colMeans(tau2_gammaCols_23_chain3)
    # tau2_gammaCols_23_chain4 <- data_samples_23$samples[, tau2_gamma_list]
    # tau2_gammaMn_23_chain4 <- colMeans(tau2_gammaCols_23_chain4)
    hist(tau2_gammaMn_23_chain1[1:30], xlab = 'posterior means 23', main = 'tau2_gamma', probability = FALSE, breaks = 100, col='blue', xlim = c(0,1), density = 5, angle = 0, border = 'blue')
    hist(tau2_gammaMn_23_chain1[31:60], col='red', xlim = c(0,5), add = TRUE, breaks = 100, density = 5, angle = 30, border = 'red')
    hist(tau2_gammaMn_23_chain1[61:90], col='orange', xlim = c(0,5), add = TRUE, breaks = 100, density = 5, angle = 60, border = 'orange')
    hist(tau2_gammaMn_23_chain1[91:120], col='green', xlim = c(0,5), add = TRUE, breaks = 100, density = 5, angle = 60, border = 'green')

  }##############################################################################################
  # hist(B[1:30], xlab = 'posterior means 23', main = 'B', probability = FALSE, breaks = 100, col='blue', xlim = c(-6,6), ylim = c(0,30), density = 5, angle = 0, border = 'blue') 
  # hist(B[31:60], col='red', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 100, density = 5, angle = 30, border = 'red')
  # hist(B[61:90], col='orange', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 100, density = 5, angle = 60, border = 'orange')
  # hist(B[91:120], col='green', xlim = c(-5,5), ylim = c(0,30), add = TRUE, breaks = 100, density = 5, angle = 60, border = 'green')
  # 
  # hist(samples1$samples[, "alpha"])
}

samples1$summary['beta0',]
samples1$summary['beta[1]',]
samples1$summary['beta[2]',]
samples1$summary['beta[3]',]

samples1$summary['beta02',]
samples1$summary['beta2[1]',]
samples1$summary['beta2[2]',]
samples1$summary['beta2[3]',]








prob_list_estim <- samples1$samples[, prob_list]
prob_list_estim_mean <- matrix(colMeans(prob_list_estim), nrow = 120, ncol = 21)

prob2_list_estim <- samples1$samples[, prob2_list]
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

precision1 <- typeI1/(typeI1 + true_positive1)
consistency1 <- typeII1/(typeII1 + true_negative1)

precision2 <- typeI2/(typeI2 + true_positive2)
consistency2 <- typeII2/(typeII2 + true_negative2)

precision1
consistency1
precision2
consistency2

prob_list_sample <- c()
prob2_list_sample <- c()
sumY1s <- rep(0,120)
sumY2s <- rep(0,120)
for (i in 1:120) {
  for (j in 1:21) {
    if (simulatedY1s[i,j] == 1) {
      sumY1s[i] <- sumY1s[i] + 1
    }
    if (simulatedY2s[i,j] == 1) {
      sumY2s[i] <- sumY2s[i] + 1
    }
  }
}
prob_list_sample <- sumY1s/21
prob2_list_sample <- sumY2s/21

par(mfrow=c(2,4))

plot(samples2$samples[ , 'prob[71]'], type = 'l', xlab = 'iteration',  ylab = expression(beta0))
abline(h = prob_list_sample[71], col = "blue") 
plot(samples2$samples[ , 'prob[72]'], type = 'l', xlab = 'iteration', ylab = expression(beta1))
abline(h = prob_list_sample[72], col = "blue") 
plot(samples2$samples[ , 'prob[73]'], type = 'l', xlab = 'iteration', ylab = expression(beta2))
abline(h = prob_list_sample[73], col = "blue") 
plot(samples2$samples[ , 'prob[74]'], type = 'l', xlab = 'iteration', ylab = expression(beta3))
abline(h = prob_list_sample[74], col = "blue") 

##########

plot(samples2$samples[ , 'prob[75]'], type = 'l', xlab = 'iteration',  ylab = expression(beta02))
abline(h = prob_list_sample[75], col = "blue") 
plot(samples2$samples[ , 'prob[76]'], type = 'l', xlab = 'iteration', ylab = expression(beta21))
abline(h = prob_list_sample[76], col = "blue") 
plot(samples2$samples[ , 'prob[77]'], type = 'l', xlab = 'iteration', ylab = expression(beta22))
abline(h = prob_list_sample[77], col = "blue") 
plot(samples2$samples[ , 'prob[78]'], type = 'l', xlab = 'iteration', ylab = expression(beta23))
abline(h = prob_list_sample[78], col = "blue") 

par(mfrow=c(2,4))

plot(samples2$samples[ , 'prob2[71]'], type = 'l', xlab = 'iteration',  ylab = expression(beta0))
abline(h = prob2_list_sample[71], col = "blue") 
plot(samples2$samples[ , 'prob2[72]'], type = 'l', xlab = 'iteration', ylab = expression(beta1))
abline(h = prob2_list_sample[72], col = "blue") 
plot(samples2$samples[ , 'prob2[73]'], type = 'l', xlab = 'iteration', ylab = expression(beta2))
abline(h = prob2_list_sample[73], col = "blue") 
plot(samples2$samples[ , 'prob2[74]'], type = 'l', xlab = 'iteration', ylab = expression(beta3))
abline(h = prob2_list_sample[74], col = "blue") 

##########

plot(samples2$samples[ , 'prob2[75]'], type = 'l', xlab = 'iteration',  ylab = expression(beta02))
abline(h = prob2_list_sample[75], col = "blue") 
plot(samples2$samples[ , 'prob2[76]'], type = 'l', xlab = 'iteration', ylab = expression(beta21))
abline(h = prob2_list_sample[76], col = "blue") 
plot(samples2$samples[ , 'prob2[77]'], type = 'l', xlab = 'iteration', ylab = expression(beta22))
abline(h = prob2_list_sample[77], col = "blue") 
plot(samples2$samples[ , 'prob2[78]'], type = 'l', xlab = 'iteration', ylab = expression(beta23))
abline(h = prob2_list_sample[78], col = "blue") 


summary(samples2$samples[ , 'prob[71]'])[2]

type1v1y1 <- 0
type1v1y2 <- 0
for (i in 1:120) {
  if (summary(samples2$samples[ , prob_list[i]])[2] <= prob_list_sample[i] && prob_list_sample[i] <= summary(samples2$samples[ , prob_list[i]])[5]) {
    type1v1y1 <- type1v1y1 + 1
  }
  if (summary(samples2$samples[ , prob2_list[i]])[2] <= prob2_list_sample[i] && prob2_list_sample[i]  <= summary(samples2$samples[ , prob2_list[i]])[5]) {
    type1v1y2 <- type1v1y2 + 1
  }
}
type1v1y1
type1v1y2

######

type1v2y1 <- 0
type1v2y2 <- 0
for (i in 1:120) {
  if (samples2$summary[prob_list[i],][4] <= prob_list_sample[i] && prob_list_sample[i] <= samples2$summary[prob_list[i],][5]) {
    type1v2y1 <- type1v2y1 + 1
  }
  if (samples2$summary[prob2_list[i],][4] <= prob2_list_sample[i] && prob2_list_sample[i]  <= samples2$summary[prob2_list[i],][5]) {
    type1v2y2 <- type1v2y2 + 1
  }
}
type1v2y1
type1v2y2















prob_list_estim <- samples1$samples[, prob_list]
prob_list_estim_mean <- matrix(colMeans(prob_list_estim), nrow = 120, ncol = 21)

prob2_list_estim <- samples1$samples[, prob2_list]
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

sample1_precision1 <- true_positive1/(typeII1 + true_positive1)
sample1_consistency1 <- true_negative1/(typeI1 + true_negative1)

sample1_precision2 <- true_positive2/(typeII2 + true_positive2)
sample1_consistency2 <- true_negative2/(typeI2 + true_negative2)

sample1_precision1
sample1_consistency1
sample1_precision2
sample1_consistency2

summary(sample1_precision1)
summary(sample1_consistency1)
summary(sample1_precision2)
summary(sample1_consistency2)




limit <- (1:100)/100
samples1_precision1_ROC <- rep(0, 100)
samples1_consistency1_ROC <- rep(0, 100)

samples1_precision2_ROC <- rep(0, 100)
samples1_consistency2_ROC <- rep(0, 100)

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
  
  samples1_precision1_ROC[k] <- sum_true_positive1/(sum_typeII1 + sum_true_positive1)
  samples1_consistency1_ROC[k] <- sum_typeI1/(sum_typeI1 + sum_true_negative1)
  
  samples1_precision2_ROC[k] <- sum_true_positive2/(sum_typeII2 + sum_true_positive2)
  samples1_consistency2_ROC[k] <- sum_typeI2/(sum_typeI2 + sum_true_negative2)
}

samples1_precision1_ROC
samples1_consistency1_ROC

samples1_precision2_ROC
samples1_consistency2_ROC