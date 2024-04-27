library(nimble, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)

set.seed(20)

simCode <- nimbleCode({
  for (i in 1:nsite) {
    for (j in 1:nvisit) {
      logit(prob[i,j]) <- beta0 + inprod(beta[1:p], x[i, 1:p]) + b[i] # True occupancy status
      logit(prob2[i,j]) <- beta02 + inprod(beta2[1:p], x[i, 1:p]) + b[i] # True occupancy status
      
      y1[i,j] ~ dbern(prob[i,j]) # Observed data
      y2[i,j] ~ dbern(prob2[i,j]) # Observed data
    }
  }
  
  beta0 ~ dnorm(0, sd = 1)
  beta02 ~ dnorm(0, sd = 1)
  for(k in 1:p) {
    beta[k] ~ dnorm(0, sd = 4)
    beta2[k] ~ dnorm(0, sd = 4)
  }
  
  for(i in 1:(nsite%/%2)) {
    b[i] ~ dnorm(2, sd = 1)
  }
  for(i in (nsite%/%2 + 1):nsite) {
    b[i] ~ dnorm(-1, sd = 1)
  }
  
})

p <- 3

X <- matrix(c(c(rep(0,30), rep(1,30), rep(0,60)),
              c(rep(0,60), rep(1,30), rep(0,30)),
              c(rep(0,90), rep(1,30))), nrow = 120, ncol = 3)
B <- c()

B[1:15] <- rnorm(15, 3, sd = 0.1)
B[16:30] <- rnorm(15, -1, sd = 0.1)
B[31:45] <- rnorm(15, 1, sd = 0.1)
B[46:60] <- rnorm(15, 1, sd = 0.1)
B[61:65] <- rnorm(5, 2, sd = 0.1)
B[66:90] <- rnorm(25, 0, sd = 0.1)
B[91:105] <- rnorm(15, -1, sd = 0.1)
B[106:120] <- rnorm(15, 3, sd = 0.1)


simModel <- nimbleModel(code = simCode,
                        constants = list(nsite = 120, nvisit = 21, p = p, alpha = 0.25,
                                         x = X), 
                        inits = list(beta0 = 1, beta = rnorm(p, 0, 1),
                                     beta02 = -1, beta2 = rnorm(p, 0, 1),
                                     b = B))


simModel$calculate('y1')
simModel$calculate('y2')


nodesToSim <- simModel$getDependencies(c("beta0", "beta", "beta02", "beta2"),
                                       self = F, downstream = T)
nodesToSim



simModel$simulate(nodesToSim)
simModel$y1
simModel$y2

simModel$beta0
simModel$beta

simModel$beta02
simModel$beta2



CsimModel <- compileNimble(simModel)



nodesToSim <- CsimModel$getDependencies(c("beta0", "beta", "beta02", "beta2"),
                                        self = F, downstream = T)
CsimModel$simulate(nodesToSim)

CsimModel$y1

CsimModel$y2





simulatedY1s <- CsimModel$y1

simulatedY2s <- CsimModel$y2


####### FIM DA SIMULAÇÃO ########



