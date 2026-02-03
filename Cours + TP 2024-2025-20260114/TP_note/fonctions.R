library(MCMCpack)
library(LearnBayes)
library(mvtnorm)
library(mnormt)



# BayesA ----


BayesA <- function(y,X,a,b,c,d,muinit,nbiter,nburn)
{
  p <- dim(X)[2]
  n <- dim(X)[1]
  
  # resultats a garder
  resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn))
  ressigma2beta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn))
  resmu <- rep(0,nbiter-nburn)
  ressigma2eps <- rep(0,nbiter-nburn)
  # initialisation
  beta <- rep(0,p)
  mu <- muinit
  sigma2beta <- rinvgamma(p,a,b) # initialisation des valeurs 
  sigma2eps <- rinvgamma(1,c,d)
  #iterations
  for (iter in 1 :nbiter)
  {
    print(iter)
    Sigmabeta <- solve(t(X)%*%X/sigma2eps + diag(1/sigma2beta))
    beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%(y-mu*rep(1,n))/sigma2eps, Sigmabeta))
    mu <- rnorm(1,t(rep(1,n))%*%(y-X%*%beta)/n,sqrt(sigma2eps/n))
    for (j in 1 :p)
    {
      sigma2beta[j] <- rinvgamma(1,a+1/2,b+1/2*beta[j]^2)
    }
    sigma2eps <- rinvgamma(1,c+n/2,d+1/2*t(y-mu*rep(1,n)-X%*%beta)%*%(y-mu*rep(1,n)-X%*%beta))
    if (iter > nburn)
    {
      resbeta[,iter-nburn] <- beta
      ressigma2beta[,iter-nburn] <- sigma2beta
      resmu[iter-nburn] <- mu
      ressigma2eps[iter-nburn] <- sigma2eps 
    }
  }
  return(list(resbeta,ressigma2beta,resmu,ressigma2eps)) 
}


# Fonction de prédiction ----

predictions <- function(X_test,mu_hat,beta_hat){
  y_pred <- mu_hat*as.vector(rep(1,dim(X_test)[1]))+ as.matrix(X_test[,]) %*% beta_hat
  return(y_pred)
}






# Fonction de sélection basée sur boxplot ----
subsetSelected <- function(coefs, seuil_bas, seuil_haut) {
  # Sélectionne les coefficients hors des seuils (outliers du boxplot)
  selected_idx <- which(coefs < seuil_bas | coefs > seuil_haut)
  return(list(
    indices = selected_idx,
    valeurs = coefs[selected_idx],
    noms = colnames(X_train)[selected_idx]
  ))
}



# ABC ----


# Prior sur mu
simMU <- function() {
  runif(1, -10, 10)
}

# Prior sur sigma (écart-type, >0)
simSIGMA <- function() {
  runif(1, 0.01, 10)
}

# Simulateur du modèle
simul <- function(n, mu, sigma) {
  rnorm(n, mean = mu, sd = sigma)
}

MeanVar <- function(obs, simu) {
  c(
    abs(mean(obs) - mean(simu)),
    abs(sd(obs)   - sd(simu))
  )
}


acceptMeanVar <- function(obs, simu, seuil) {
  d <- MeanVar(obs, simu)
  as.integer(d[1] < seuil[1] & d[2] < seuil[2])
}


ABC_MeanVar <- function(obs, nIter, seuil) {
  
  n <- length(obs)
  
  decision     <- integer(nIter)
  sampledMU    <- numeric(nIter)
  sampledSIGMA <- numeric(nIter)
  
  for (i in 1:nIter) {
    mu    <- simMU()
    sigma <- simSIGMA()
    
    simuDATA <- simul(n, mu, sigma)
    
    decision[i] <- acceptMeanVar(obs, simuDATA, seuil)
    sampledMU[i] <- mu
    sampledSIGMA[i] <- sigma
  }
  
  data.frame(
    decision = decision,
    mu = sampledMU,
    sigma = sampledSIGMA
  )
}



quant_stat <- function(x) {
  quantile(x, probs = c(0.1, 0.5, 0.9))
}

comparequant <- function(q1, q2) {
  sqrt(sum((q1 - q2)^2))
}

acceptQuant <- function(obs, simu, seuil) {
  d <- comparequant(quant_stat(obs), quant_stat(simu))
  as.integer(d < seuil)
}



ABC_Quant <- function(obs, nIter, seuil) {
  
  n <- length(obs)
  
  decision     <- integer(nIter)
  sampledMU    <- numeric(nIter)
  sampledSIGMA <- numeric(nIter)
  
  for (i in 1:nIter) {
    mu <- simMU()
    sigma <- simSIGMA()
    
    simuDATA <- simul(n, mu, sigma)
    
    decision[i] <- acceptQuant(obs, simuDATA, seuil)
    sampledMU[i] <- mu
    sampledSIGMA[i] <- sigma
  }
  
  data.frame(
    decision = decision,
    mu = sampledMU,
    sigma = sampledSIGMA
  )
}

