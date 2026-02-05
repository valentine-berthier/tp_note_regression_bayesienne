# TP0
# ENSAI
# Library utiles  
library(mvtnorm) # pour simuler des vecteurs gaussiens
library(mnormt)
library(invgamma) # pour simuler une inverse gamma

####################################################### 
# Simulations des donn?es que l'on ?tudiera 
set.seed(123)
nobs=200
simuvarexpl <- matrix(rep(0,300*nobs),ncol=300,nrow=nobs)
for (i in 1 :300) simuvarexpl[,i] <- runif(nobs,-5,5)
#simuvarexpl = as.data.frame(simuvarexpl)
simuvarexpl=as.matrix(simuvarexpl) 
trueind <- c(10,20,30,40,50) #5 variables qui influencent Y 
beta <- c(1,-1,2,-2,3)
ysimu <- simuvarexpl[,trueind]%*% beta + rnorm(nobs,0,2) 


# remarque : on peut normaliser (centrer réduire) les covariables 
# ainsi les beta sont comparables entre eux
simuvarexpl=scale(simuvarexpl)
###############################################################
# Fonction de pr?diction 
predictions <- function(matableTest,muchap,betachap)
{
  ychap <- muchap * rep(1,dim(matableTest)[1])+ as.matrix(matableTest[,]) %*% betachap
  return(ychap)
}

########################################################################
########################################################################
# Question 2
##############Simulation d'un coefficient aléatoire
# on a calculé la loi de beta | y (c'est une gaussienne)
# on voit dans la variance a posteriori que le bayésien fonctionne comme du ridge 
# la matrice devient inversible 
# on va la simuler (ici le Gibbs n'a qu'une seule variable !!!)
# Initialisaton des paramètres 
s2eps=4  # supposé connu ici 
s2beta= 1 # a priori ici ?a influence le résultat 
# un s2beta trop petit va trop contraindre les beta 
# et ils auront du mal à s'éloigner de zéro
# pour aller vite on prend des tailles de burn in et d'itérations petites
niter=1000
nburn=200

###############################################################
# ci dessous la fonction qui simule le beta a posteriori 

BayesBasic  <- function(y,X,sigma2beta,sigma2eps,nbiter,nburn)
{
  p <- dim(X)[2]
  n <- dim(X)[1]
  identite=diag(1,p,p)
  ### resultats a garder
  resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,
                    ncol=(nbiter-nburn))
  ### initialisation
  beta <- rep(0,p)
  ### iterations
  for (iter in 1 :nbiter)
  {
    #print(iter)
    Sigmabeta <- solve(t(X)%*%X/sigma2eps + identite/sigma2beta)
    beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%y/sigma2eps,
                              Sigmabeta))
    if (iter > nburn)
    {
      resbeta[,iter-nburn] <- beta
    }
  }
  return(resbeta)
}
#### Illustrations 
# si on prend s2beta trop grand on a des beta trop variables
# et inversement 
s2beta=1
nburn=10 #ne sert à rien 
niter=210
res=BayesBasic(ysimu,simuvarexpl,s2beta,4,niter,nburn)
plot((rowMeans(res)))
plot(sort(abs(rowMeans(res))))
res=BayesBasic(ysimu,simuvarexpl,200,4,300,10)
# => trop dispers? 
res=BayesBasic(ysimu,simuvarexpl,0.1,4,300,10)
plot(sort(abs(rowMeans(res))))
# Interprétation 
# 5 coefficients ressortent fortement 
# 
plot(res[1,])
# Interprétation 
# beta1 simulé 300 fois 

plot(density(res[1,]))
# éparpillé autour de 0 



#s2 plus grand 
# les points se détachent 
#s2 plus petit 
# empêche les beta d'être trop éparpillés 
# tous les coeff (beta) vont être très ressérés et les vrais seront
# sous estimés en valeur absoolue




plot(density(res[50,]))
mean(res[50,])
plot(density(res[51,])) 
# si on prend s2beta = 0.1 alors tous les coef vont ?tre tr?s resserr?s 
# et les vrais seront sous estim?s en valeur absolue 
# Inversement si on prend s2beta= 1000 alors les coef vont ?tre plus dispers?s 
# et ils seront surestim?s en valeur absolue  
niter=500
nburn=100
# attention on devrait prendre seulement les 100 premi?res observations (training)
resbayes=BayesBasic(ysimu[1:100],simuvarexpl[1:100,],s2beta,s2eps,niter,nburn)
plot(density(resbayes[1,])) # bruit 
plot(density(resbayes[5,])) # bruit 
mean(resbayes[5,])
plot(density(resbayes[40,])) # vrai coef
mean(resbayes[40,]) # vrai coef 
plot(density(resbayes[50,])) # vrai coef
mean(resbayes[50,]) # vrai coef
#rowMeans(resBayes)
#colMeans(t(resbayes))
plot(sort(abs(colMeans(t(resbayes))))) # on voit les coef qui se d?gagent 
# on voit les 5 coef sortir du lot 

########################################################################
########################################################################
# Question 3
#### Algorithme avec variance de Zellner
# Pour extraire une  matrice inversible on peut prendre les premi?res variables
# on peut aussi se r?f?rer ? l'indice de conditionnement
# Par exemple calculer les valeurs propres
X=as.matrix(simuvarexpl[1:100,]) 
valpX=eigen(t(X)%*%X, symmetric =TRUE, only.values = TRUE)
valpX
valmax=max(valpX$values)
valmin=min(valpX$values)
indcond=valmax/valmin
indcond
# on trouve ici quelque chose de n?gatif car ce sont les approximations 
# (normalement ?a vaut z?ro) 
# Prenons une sous matrice de taille 50*50 par exemple 
indalea=sample(c(1:300),50,replace=FALSE)
X2=X[,indalea] 


# c = coefficiet d'?chelle 
# L = coefficient ridge
BayesZellner  <- function(y,X,c ,L, sigma2beta,sigma2eps,nbiter,nburn)
{
  p <- dim(X)[2]
  n <- dim(X)[1]
  identite=diag(1,p,p)
  ### resultats a garder
  resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,
                    ncol=(nbiter-nburn))
  ### initialisation
  beta <- rep(0,p)
  ### iterations
  for (iter in 1 :nbiter)
  {
    print(iter)
    Sigmabeta <- c*solve(t(X)%*%X/sigma2eps + L*identite/sigma2beta)
    beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%y/sigma2eps,
                              Sigmabeta))
    if (iter > nburn)
    {
      resbeta[,iter-nburn] <- beta
    }
  }
  return(resbeta)
}
# Algorithme avec Zelnner (+ Ridge ? => question suppl?mentaire) 
s2eps=4
s2beta=200
# si on prend s2beta = 0.1 alors tous les coef vont ?tre tr?s resserr?s 
# et les vrais seront sous estim?s en valeur absolue 
# Inversement si on prend s2beta= 1000 alors les coef vont ?tre plus dispers?s 
# et ils seront surestim?s en valeur absolue  
niter=300
nburn=100
c=1
L=100
# attention ici on prend seulement les 100 premi?res observations (training)
resbayes=BayesZellner(ysimu[1:100],as.matrix(simuvarexpl[1:100,]),c,L, s2beta,s2eps,niter,nburn)
######################################################################
plot(density(resbayes[1,])) # bruit 
plot(density(resbayes[5,])) # bruit 
mean(resbayes[5,])
plot(density(resbayes[40,])) # vrai coef
mean(resbayes[40,]) # vrai coef 
plot(density(resbayes[50,])) # vrai coef
mean(resbayes[50,]) # vrai coef
#colMeans(t(resbayes))
plot(sort(abs(colMeans(t(resbayes))))) # on voit les coef qui se d?gagent 
# Rermarque : on pourrait ajouter un coefficient ridge pour stabiliser l'inverse.


#######################################################
#######################################################
# Question 4
# Algorithme avec  s2 al?atoire ? 
# prenons sigma2eps al?atoire de loi InvGamma(c,d) 
# On montre que sa loi a posteriori est une InvGamma de param?tre 
# (c+n/2,d+(Y-mu-X beta)'(Y-mu-X beta)/2)
# on programme le Gibbs ainsi : 
library(invgamma)
BayesBasic2  <- function(y,X,sigma2beta,c,d,nbiter,nburn)
{
  p <- dim(X)[2]
  n <- dim(X)[1]
  identite=diag(1,p,p)
  ### resultats a garder
  resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn)) # param?tre vectoriel
  reseps <- rep(0,(nbiter-nburn)) # param?tre univari? 
  ### initialisation
  beta <- rep(0,p) # on aurait pu simuler un vecteur gaussien N(0,s2beta I)
  sigma2eps <- var(y)   # on intialise avec l'empirique
  # on aurait pu aussi simuler une inv gamma (c,d)
  ### iterations
  for (iter in 1 :nbiter)
  {
    print(iter)
    Sigmabeta <- solve(t(X)%*%X/as.numeric(sigma2eps) + identite/sigma2beta)
    # premi?re marge conditionnelle sachant s2eps
    beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%y/as.numeric(sigma2eps),Sigmabeta))
    s=t(y-X%*%beta)%*%(y-X%*%beta)
    # deuxi?me loi conditionnelle sachant beta
    sigma2eps <- rinvgamma(1,c+n/2,d+s/2)
    if (iter > nburn)
    {
      resbeta[,iter-nburn] <- beta
      reseps[iter-nburn] <- sigma2eps
    }
  }
  resu=list(resbeta,reseps)
  return(resu)
}
############ Illustrations 
# on fixe c=1=d pour l'inverse Gamma => vague car esp?rance et variance = infini
c=d=1  # => la variance du prior vaut l'infini !
# moyenne inverse gamma(c,d) = c/(d-1)
# variance inverse gamma(c,d) = c^2/((d-1)^2(d-2)) 
sigma2beta=100 # a priori sur la variance de beta
res2=BayesBasic2(ysimu,simuvarexpl,sigma2beta,c,d,600,300)
resbeta=res2[[1]]
reseps=res2[[2]] 
plot(resbeta[1,])
plot(resbeta[50,])
plot(reseps)
plot(density(resbeta[1,]))
plot(density(resbeta[50,]))
plot(density(resbeta[40,]))
plot(density(reseps))
mean(resbeta[50,]) # moyenne a posteriori qui est l'estimateur bay?sien
mean(resbeta[40,])
mean(reseps) # moyenne a posteriori 
plot(sort(abs(colMeans(t(resbeta)))))
plot(reseps)
plot(resbeta[50,])

########################################################################
########################################################################
# Question 5
# Ici on est tr?s proche du Bayes A qui sera ?tudier au TP1 
# ce qui change ici c'est que mu est suppos?e ?gale ? 0
# Mais on va le faire 
BayesBasic3  <- function(y,X,a,b,c,d,nbiter,nburn)
{
  p <- dim(X)[2]
  n <- dim(X)[1]
  identite=diag(1,p,p)
  ### resultats a garder
  resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn)) # param?tre vectoriel
  reseps <- rep(0,(nbiter-nburn)) # param?tre univari? 
  resvarbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn))
  ### initialisation
  beta <- rep(0,p)
  sigma2beta = rinvgamma(1,a,b)
  sigma2eps <- var(y)   # on intialise avec l'empirique
  ### iterations
  for (iter in 1 :nbiter)
  {
    print(iter)
    Sigmabeta <- solve(t(X)%*%X/as.numeric(sigma2eps) + identite/sigma2beta)
    # premi?re marge conditionnelle sachant s2eps (et Y)
    beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%y/as.numeric(sigma2eps),Sigmabeta))
    s=t(y-X%*%beta)%*%(y-X%*%beta)
    # deuxi?me loi conditionnelle sachant beta (et Y)
    sigma2eps <- rinvgamma(1,c+n/2,d+s/2)
    # troisi?me loi conditionnelle sachant beta (et Y)
    sigma2beta  = rinvgamma(p,a+1/2,b+1/2*beta^2)
    if (iter > nburn)
    {
      resbeta[,iter-nburn] <- beta
      reseps[iter-nburn] <- sigma2eps
      resvarbeta[,iter-nburn] <- sigma2beta
    }
  }
  resu=list(resbeta,reseps,resvarbeta)
  return(resu)
}
a=b=c=d=1  # tr?s vague car exp?rance et variance infinies 
res3=BayesBasic3(ysimu,simuvarexpl,a,b,c,d,800,300)
resbeta=res3[[1]]
reseps=res3[[2]]
resvarbeta=res3[[3]] 
par(mfrow=c(1,1))
plot(resbeta[1,])
plot(resbeta[50,])
plot(reseps)
plot(resvarbeta[2,])
plot(resvarbeta[40,])
plot(resvarbeta[200,])
plot(density(resbeta[1,]))
plot(density(resbeta[50,]))
plot(density(resbeta[40,]))
plot(density(reseps))
plot(density(resvarbeta[1,]))
plot(density(resvarbeta[50,]))
plot(density(resvarbeta[40,]))
mean(resbeta[50,]) # moyenne a posteriori qui est l'estimateur bay?sien
mean(resbeta[40,])
mean(reseps) # moyenne a posteriori 
plot(sort(abs(colMeans(t(resbeta)))))
# on rep?re les 5 coefficients les plus significatifs du mod?le 
# on peut regarder leurs trajectoires, leurs densit?s estim?es 
par(mfrow=c(2,1))
plot(resbeta[10,])
plot(resbeta[20,])
plot(resbeta[30,])
plot(resbeta[40,])
plot(resbeta[50,])
plot(density(resbeta[10,]))
plot(density(resbeta[20,]))
plot(density(resbeta[30,]))
plot(density(resbeta[40,]))
plot(density(resbeta[50,]))
# Critiques : il faudrait prendre 1) un burn-in plus grand (10.000) 
# 2) prendre une silu sur 10 (pour diminuer la d?pdendance) et avoir 
# ainsi une ?chantillon "presque" "i"id pour construire par exemepl des intervalles
# de confiance a posteriori => c'est un moyen de s?lectionner les betas 
# significatifs (ceux dont 0 n'est pas dans l'intervalle). 
# Le mod?le repose sur la normalit? des Y|beta,sigma 
# mais c'est assez robuste ? la non-normalit?  
# une analyse de sensibilit? serait de faire varier les hyperparam?tres : 
# a,b,c,d (les prendre plus grands...) 

########################################################################


########################################################################
# Question 6
########################### # R?gression bay?sienne et ABC 
# Illustration avec un mod?le de scale mixture
# Simulation du mod?le
# on pose 
# Y_i = mu_i + X_i beta_i + e_i 
# avec e ~ N(0,s)
# mu ~uniforme (-a,a)
# beta_i ~ N(m,v) 
#  
#####################################################
# Les param?tres 
n=200 # nbre observations (simulations)
d=1 # taille beta 
m=6
v=1
s=1
m0=0
a=5
b=5
x <- matrix(rep(0,d*n),ncol=d,nrow=n)
for (i in 1:d){
  x[,i] <- runif(n,-5,5)}
x <- as.data.frame(x)
X=as.matrix(x) 
beta <- rnorm(n,m,v) # chaque beta est al?atoire 
m0 <- runif(n,-a,a)
s <- rchisq(n,b)
y <- m0 + X%*%beta + rnorm(n,0,s)  # nos y simul?s 



############################################################################################
# ABC basique 
############################################################################################
ABCbasic  = function(y, X, a, b, V, seuil, K)
{ 
  d = dim(X)[2]
  identite=diag(rep(1,d)) # matrice identit? =1 (dimension =1)
  identn=diag(rep(1,n))    # matrice identit? n x n 
  initbeta=rep(0,(d*K))
  resbeta=matrix(initbeta,ncol=K)
  resmu=rep(0,K)
  ress2=rep(0,K)
  n=length(y)
  compteur=0
  for (j in 1:(K-1))
  {
    mu=runif(1,-a,a) # a priori uniforme 
    s2=rchisq(1,b) # a priori sur la variance du bruit 
    s=sqrt(s2)
    beta=rmnorm(1,0,V*identite ) # a priori sur beta 
    m=mu+X%*%beta  # moyenne des y sachant les param?tres
    res = rmnorm(1,m,s*identn) # simulation du y sachant les param?tres 
    if (sqrt(sum((y-res)^2)/n) < seuil)  # acceptation rejet 
    {
      compteur = compteur+1  
      resbeta[,compteur] =beta 
      ress2[compteur] =s2 
      resmu[compteur] = mu 
    }
  }
  accept=sum(resmu!=0) # nbre d'acceptation 
  betafinal=resbeta[,c(1:accept)]
  resmufinal=resmu[c(1:accept)]
  ress2final=ress2[c(1:accept)]
  return(list(mu=resmufinal,beta=betafinal,s2=ress2final, taux=100*accept/(K-1)))
}
##################################################
a=5
b=5
V=40
seuil=8
K=20000
res = ABCbasic(y, X, a, b, V, seuil, K)
res$taux
mean(res$mu)
mean(res$s2)
mean(res$beta)
plot(density(res$beta))
# le d?s?quilibre peut venir du fait que l'on pr?sente plus souvent des beta 
# ? gauche qu'? droite de 5 
# il s'estompe si on diminue le seuil ou si on augmente V
# il y a un pb d'identifiabilit? avec mu !!! 
# il faut enlever mu ou alors reporter la moyenne de beta sur mu !!! 


# TP1
# ENSAI
library(mvtnorm)
library(mnormt)
# Remarque : il faut normaliser le design X en faisant X <- scale(X)
# Ici les X sont toutes de m?me loi 
# Simulations des données



# ETAPE 1 
simuvarexpl <- matrix(rep(0,300*200),ncol=300,nrow=200)
for (i in 1 :300) simuvarexpl[,i] <- runif(200,-5,5)
simuvarexpl = as.data.frame(simuvarexpl)
trueind <- c(10,20,30,40,50)
beta <- c(1,-1,2,-2,3)
ysimu <- as.matrix(simuvarexpl)[,trueind]%*% beta + rnorm(200,0,2)
###############################################################
# Fonction de pr?diction 
predictions <- function(matableTest,muchap,betachap)
{
  ychap <- muchap * rep(1,dim(matableTest)[1])+ as.matrix(matableTest[,]) %*% betachap
  return(ychap)
}
############################
subsetSelected <- function(resAlgo, varexpl, mini, maxi) 
{
  numselected <- c(which(resAlgo < mini), which(resAlgo > maxi))
  selected <- character()
  valeurs <- numeric()
  for (i in 1 :length(numselected)) 
  {
    selected[i] <- names(varexpl)[numselected[i]]
    valeurs[i] <- resAlgo[numselected[i]] 
  }
  subset <- cbind(selected, valeurs)
  subset <- as.data.frame(subset)
  subset$valeurs <- as.numeric(as.vector(subset$valeurs))
  return(subset) 
}

# Remarque : il faudrait faire de la validation crois?e en 
# changeant le jeu de donn?es TEST 

#####################Mod?le RR 
# Y = 1*mu + X beta + eps # notation du cours 
# avec les notations du package mu = beta 
# beta = u 
# X = Z 
# le 1 = X 
# Y = X*beta + Z u + e  
library(rrBLUP)

######## le mod?le sur le training (les 100 premi?res observations)
resBLUP <- mixed.solve(ysimu[1 :100], Z = as.matrix(simuvarexpl[1 :100, ]), X = as.matrix(rep(1,
                                                                                              100)), method = "REML")


# les notations de mixed.solve sont : 
# notre mu = beta 
# notre beta = u 
# notre epsilon = e
######## les estimations 
muchap <- resBLUP$beta  # moyenne  
muchap
betachap <- resBLUP$u   # les coefficients 
par(mfrow=c(1,1))
plot(sort(abs(betachap)))   # on voit un fort shrinkage 


# Interprétation 
# - fort effet shrinkage 
betachap[50]
resBLUP$Ve  # variance des epsilon
resBLUP$Vu # variance des betas
# on a des variances tr?s instables  
# des beta sous-etim?s !!! 

################### Les pr?dictions sur les donn?es du jeu test (les 100 derni?res valeurs)
predBLUP <- predictions(simuvarexpl[101 :200, ], muchap, betachap)

################ corr?lation en y pr?dits et y observ?s (sur jeu test)
cor(ysimu[101 :200], predBLUP)
# 0.5659269
plot(predBLUP ~ ysimu[101 :200])
abline(0,1)
# on s'aper?oit de l'effet shrinkage 

############ S?lection des coefficients 
plot(betachap,xlab="Indice variable",ylab="Param?tre de r?gression associ?",main="rrBLUP")
plot(sort(abs(betachap)))
boxplot(betachap)
bb <- boxplot(betachap)
varselectedBLUP <- subsetSelected(betachap,simuvarexpl,bb$stats[1,], bb$stats[5,])
varselectedBLUP
#On retient les valeurs suivantes 

#  selected    valeurs
#1      V40 -0.8154402
#2      V61 -0.3654459
#3      V30  0.8491002
#4      V50  1.1367013


################ Bayes A
library(MCMCpack) 
library(LearnBayes)
library(invgamma)
library(rmnorm)
################## le mod?le
# y|beta, sigmaSeps, mu ~N(mu + X beta , sigma2eps)
# u uniforme 
# beta ~N(0, sigma2(1), ..., sigma2(p)) 
# sigma2(i) ~ InvGamma(a,b) iid
# sigma2eps ~InvGamma(c,d)
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

####################### Application sur nos donn?es
# les a priori sont des inverses gamma 
# de param?tres (a,b) et (c,d) 
# leurs moyennes et variances sont 
a=c=1 ; b=d=1
b/(a-1) ; d/(c-1)   # esp?rance d'une inv gamma
b^2/((a-1)^2*(a-2)) ; d^2/((c-1)^2*(c-2)) # variance inv gamma
# avec ces valeurs ont met des a priori tr?s vague, de variance infinie 
priora <- 2
priorb <- 1
priorc <- 2
priord <- 1
resBAYESA <- BayesA(ysimu[1 :100], as.matrix(simuvarexpl[1 :100, ]), priora, priorb, priorc,
                    priord, mean(ysimu[1 :100]), 2000, 1500)
moybeta <- apply(resBAYESA[[1]], 1, mean) # estimateurs des beta
moysigma2beta <- apply(resBAYESA[[2]], 1, mean)
moymu <- mean(resBAYESA[[3]])
moysigma2eps <- mean(resBAYESA[[4]])
moybeta[50]
plot(moybeta)
plot(sort(abs(moybeta)))
moymu
plot(moysigma2beta)
moysigma2eps
# regardons les trajectoires : 

png('test.png')
par(mfrow=c(2,2))
plot(resBAYESA[[1]][50, ])  # pour beta50
plot(resBAYESA[[2]][1, ])   # pour s2beta1
plot(resBAYESA[[3]])        # pour mu
plot(resBAYESA[[4]])        # pour s2eps 
dev.off()
# en conclusion il faut un plus grand burn-in. 50.000 par exemple. 
# ou alors lancer plusieurs cha?nes et comparer les trajectoires (empirique)
# dans tous les cas il peut y avoir avoir une d?pendance des variables g?n?r?es 
# par l'algo. 

############# Comparaison des lois a priori et a posteriori
tracesigmabeta1 <- resBAYESA[[2]][30,]
plot(density(tracesigmabeta1),col="blue",lwd=2, main="sigma_beta_1 : Bayes A")
# attention il faudrait prendre "une variable sur 10" car ce n'est pas iid (c'est juste id)
curve(dinvgamma(x,priora,priorb), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))

tracesigmabeta50 <- resBAYESA[[2]][50,]
plot(density(tracesigmabeta50),col="blue",lwd=2, main="sigma_beta_50 : Bayes A", ylim=c(0,2))
curve(dinvgamma(x,2,1), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))

tracesigmaeps <- resBAYESA[[4]]
plot(density(tracesigmaeps),col="blue",lwd=2, main="sigma_eps : Bayes A")
curve(dinvgamma(x,priorc,priord), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))


# Int?r?t en bay?sien s?quentiel : si on fait le m?me mod?le avec des donn?es diff?rentes 
# on peut utiliser les lois a posteriori comme a priori ensuite... 
# => ?a rend l'algorithme compl?tement automatique 
# => pas d'hyperparam?tre ? r?gler 

# On pourrait faire la m?me chose si l'?chantillon est assez grand 
# on en prend 10% par exemple pour r?cup?rer des lois a posteriori 
# Et ensuite on travaille avec ces lois comme a posteriori sur les 90% des donn?es restantes 

##################### Pr?dictions 
predBAYESA <- predictions(simuvarexpl[101 :200, ], moymu, moybeta)  # sur le jeu test
cor(ysimu[101 :200], predBAYESA)
# 0.6516713
plot(predBAYESA ~ ysimu[101 :200])  # pr?dictions versus observations sur jeu test 
abline(0,1)

#################### S?lection des betas
plot(moybeta,xlab="Indice variable",ylab="Param?tre de r?gression associ?",main="Bayes
A")
plot(sort(abs(moybeta)))
bb <- boxplot(moybeta)
varselectedBAYESA <- subsetSelected(moybeta,simuvarexpl, bb$stats[1,],bb$stats[5,])
varselectedBAYESA
# Valeurs
# 1      V40 -0.5848762
# 2     V275 -0.3434829
# 3     V299 -0.3073950
# 4      V10  0.4964229
# 5      V30  0.7446974
# 6      V50  1.4801090
# 7      V54  0.3512319
# 8      V97  0.3584435

####################### LASSO Bay?sien 
# On utilise le package BLR(Y,XF,XR,XL,Z,prior,nIter, burnIn
# car Y = mu*1 + XF betaF + XR betaR + XL betaL + Z u + epsilon 
# avec betaF = fixe (= mod?le classique, ici non) 
# mu = constante 
# betaR = ridge (= bayes A, avec plusieurs variance, ici non)  
# betaL = Lasso (ici oui) 
# u = effet al?atoire du RR (avec une seule variance, ici non)

############# le mod?le 
library(BLR)
LASSO_BLR <- BLR(y=ysimu[1 :100],XL=as.matrix(simuvarexpl[1 :100,]), prior=list(varE=list(df=2,S=1),
                                                                                lambda=list(shape=10, rate=0.1, type='random',value=2)), nIter=3000,burnIn=2000, saveAt="BLR_",
                 thin=1)#,thin2=1)
# thin = 1 signifie qu'il prend toutes les simus (1/10 par d?faut)
# thin = 5 => il garde une simu sur 5 
# attention il garde m?me le burn in 
# 
# remarque : la variance de lambda^2 ici vaut shape/rate^2 
# plus on prend rate petit et plus la variance va ?tre ?lev?e et lambda pourra ?voluer librement
# VarE est la variance de epsilon qui suit une inverse khi deux 
# lambda de type 'random', value=2, signifie que lambda^2 suit une gamma 
# on peut m?me fixer la valeur de lambda, par exemple issue d'un LASSO classique
LASSO_BLR$mu # muchap = moyenne des trajectoires = estimation de la moyenne a posteriori
LASSO_BLR$varE # = moyenne des trajectoires = estimation de la moyenne a posteriori
LASSO_BLR$bL   # betachap = moyenne des trajectoires = estimation de la moyenne a posteriori
LASSO_BLR$tau2  # = moyenne des trajectoires = estimation de la moyenne a posteriori
LASSO_BLR$lambda # = moyenne des trajectoires = estimation de la moyenne a posteriori

# remarque : il serait int?ressant de faire varier le shape et le rate de la loi (a priori) gamma du 
# param?tre LASSO lambda. On regarderait alors si la corr?lation est meilleure

################## Trajectoires simul?es e

# Attention : 
# si on relance alors les trajecgoires se cumulent (concat?nation des r?sultats)
# il faut soit effac? BLR_, soit changer le nom...
# Par ailleurs, on ne garde qu'une valeur sur 10 (par d?faut, pour garantir une "ind?pendance") 

tracevarE <- scan('BLR_varE.dat')
plot(tracevarE,type='o',xlab="iteration",ylab="varE",main="trace varE")
# on observe une rupture => ? explorer sur une cha?ne plus longue
# ou bien sur plusieurs cha?nes 
plot(density(tracevarE),col="blue",lwd=2, main="varE : LASSO bay?sien")
curve(dinvgamma(x,1,1), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))


tracelambda <- scan('BLR_lambda.dat')
plot(tracelambda,type='o',xlab="iteration",ylab="lambda",main="trace lambda")
plot(density(tracelambda),xlim=c(0,5),col="blue",lwd=2, main="lambda : LASSO bay?sien")
curve(dgamma(x,1,0.1), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))


#####################Pr?diction 
predLASSOBLR <- predictions(simuvarexpl[101 :200, ], LASSO_BLR$mu, LASSO_BLR$bL)
cor(ysimu[101 :200], predLASSOBLR)
#  0.9522485 
plot(predLASSOBLR ~ ysimu[101 :200])
abline(0,1)


###################### S?lection
plot(LASSO_BLR$bL,xlab="Indice variable",ylab="Param?tre de r?gression associ?", main="LASSO
bay?sien")
plot(sort(abs(LASSO_BLR$bL)))
bb <- boxplot(LASSO_BLR$bL)
bb
varselectedLASSO <- subsetSelected(LASSO_BLR$bL,simuvarexpl, bb$stats[1,],bb$stats[5,])
varselectedLASSO
# on va faire du seuillage ? la main pour ?liminer quelques variables 
varselectedLASSO <- subsetSelected(LASSO_BLR$bL,simuvarexpl,-0.3,0.3)
varselectedLASSO
d     valeurs
1       V20 -0.53605402
2       V40 -1.52697107
3      V272 -0.11674226
4      V275 -0.09606266
5      V297 -0.11999264
6       V10  0.75508769
7       V30  1.70948647
8       V50  2.78763941
9      V242  0.11050764
10     V260  0.10767629
11     V273  0.11312140

# Pour s?lectionner moins de variable il faudrait un srinkage plus fort, 
# ce qui correspond ? un param?tre lambda du LASSO plus ?lev?
# Ici lambda^2 suit une Gamma(e,f) avdc (e,f) = (1,0.1)
# E(lambda^2)=e/f  et V(lamdba^2) = e^2/f
# Prenons e + grand (=10) 

################# SSVS 
####### Mod?lisation 
selection_SSVS <- function (vardep, varexpl, nbiter, nburn, lec, nbSelecInit,nbToChange,Pi)
{
  X <- as.matrix(varexpl)
  y <- as.numeric(as.vector(vardep))
  y <- y - mean(y)
  nind <- dim(X)[1]
  nvar <- dim(X)[2]
  sumgamma <- rep(0, nvar)
  nbactu <- 0
  nbselecactu <- numeric()
  indgamma10 <- sample(c(1:nvar), nbSelecInit, replace = FALSE)
  gamma0 <- rep(0, nvar)
  for (i in 1:nbSelecInit) {
    gamma0[indgamma10[i]] <- 1
  }
  indgamma1 <- indgamma10
  gamma <- gamma0
  nbSelec <- nbSelecInit
  
  for (iter in 1:nbiter) {
    print(iter)
    gammaprop <- gamma
    indgamma1prop <- indgamma1
    indToChange <- sample(c(1:nvar), nbToChange, replace = FALSE)
    for (i in 1:nbToChange){
      if (gamma[indToChange[i]]==0){
        gammaprop[indToChange[i]] <- 1
        indgamma1prop <- c(indgamma1prop,indToChange[i])
      }
      else {
        gammaprop[indToChange[i]] <- 0
        indremove <- which(indgamma1prop==indToChange[i])
        indgamma1prop <- indgamma1prop[-indremove]
      }
    }
    nbSelecprop <- length(indgamma1prop)
    if (nbSelecprop==0){ # condition pour empecher gamma avec que des 0
      cond <- 0
      while(cond==0){
        gammaprop <- gamma
        indgamma1prop <- indgamma1
        indToChange <- sample(c(1:nvar), nbToChange, replace = FALSE)
        for (i in 1:nbToChange){
          if (gamma[indToChange[i]]==0){
            gammaprop[indToChange[i]] <- 1
            indgamma1prop <- c(indgamma1prop,indToChange[i])
          }
          else {
            gammaprop[indToChange[i]] <- 0
            indremove <- which(indgamma1prop==indToChange[i])
            indgamma1prop <- indgamma1prop[-indremove]
          }
        }
        nbSelecprop <- length(indgamma1prop)
        if (nbSelecprop>0){cond <- 1}
      }
    }
    indgamma1 <- which(gamma == 1)
    nbSelec <- length(indgamma1)
    Xgamma <- X[, indgamma1]
    Xgammaprop <- X[, indgamma1prop]
    temp <- (t(y)%*%(diag(rep(1,nind))-lec/(1+lec)*Xgammaprop%*%solve(t(Xgammaprop)%*%Xgammaprop)%*%t(Xgammaprop))%*%y)/(t(y)%*%(diag(rep(1,nind))-lec/(1+lec)*Xgamma%*%solve(t(Xgamma)%*%Xgamma)%*%t(Xgamma))%*%y)
    A <- (1+lec)^((nbSelec-nbSelecprop)/2)*(Pi/(1-Pi))^(nbSelecprop-nbSelec)*temp^(-(nind-1)/2)
    probaccept1 <- min(1,A)
    seuil <- runif(1)
    if (seuil < probaccept1){
      gamma <- gammaprop
      indgamma1 <- indgamma1prop
      nbSelec <- nbSelecprop
      nbactu <- nbactu+1
      nbselecactu <- c(nbselecactu,nbSelec)
    }
    if (iter > nburn) {
      sumgamma <- sumgamma + gamma
    }
  }
  return(list(sumgamma,nbactu,nbselecactu))
}
#### Initialisation 
nbinit <- 10 # nombre de beta non nuls au d?part 
nbToChange <- 2 # nombre de beta que l'on propose de changer ? chaque fois 
Pi <- 10/300   # la probabilit? de choisir un beta (en moyenne 10 car il y en a 300)
lec <- 50  # pr?conisation : entre 10 et 100 
resSSVS <- selection_SSVS(ysimu[1 :100], as.matrix(simuvarexpl[1 :100, ]), 3000, 2000, lec,
                          nbinit, nbToChange, Pi)
resSSVS[[2]]  # on voit combien il y a eu d'?changes dans les 1000 derniers run
resSSVS[[3]]  # on voit le nbre de variables s?lectionn?es lors de ces ?changes 
resSSVS[[1]]
####### S?lection 
plot(sort(abs(resSSVS[[1]])))
boxplot(resSSVS[[1]])
plot(resSSVS[[1]],xlab="Indice variable",ylab= "Nombre de s?lections post-burn-in", main="SSVS")
varselectedSSVS <- subsetSelected(resSSVS[[1]],simuvarexpl,0,999)
varselectedSSVS
which(resSSVS[[1]]>400)
# 10  20  30  40  50  90 133 147
# On a s?lectionn? les 5 bonnes variables. 
# Les variables suppl?mentairtes ne sont pas celles du LASSO
# Or on sait que le LASSO a tendance ? ne pas s?lectionner des variables 
# corr?l?es. C'est peut ?tre une explication... 

# on peut changer pi 
# en fait le burn-in permet de chercher les bons betas. Une fois qu'ils 
# sont s?l?ctionn?s ils restent car les autres beta ont une densit?s conditionnelle a 
# posteriori moins forte.  




##################### 
# Mod?le lin?aire simple avec les var retenues LASSO vs SSVS 
indLASSO=c(10,20,30,40,50,78,150, 154, 192, 258, 273) 
indSSVS=c(10,20,30,40,50,90,133, 147) 
resLASSO=lm(ysimu~as.matrix(simuvarexpl[,indLASSO])) 
resSSVS=lm(ysimu~as.matrix(simuvarexpl[,indSSVS])) 
summary(resLASSO)
summary(resSSVS)



