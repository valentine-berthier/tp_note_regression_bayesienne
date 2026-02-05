getwd()
setwd("/home/id2467/tp_note_regression_bayesienne/Cours + TP 2024-2025-20260114/TP_note")
getwd()

telecat <- read.csv(file = "telecat.csv")

set.seed(1234)


telecat <- read.csv("telecat.csv")

ytelecat <- telecat$Y
Xtelecat <- telecat[, colnames(telecat) != "Y"]
dim(Xtelecat)

# Train/test split
train_indices <- sample(1:150, size = 100)

y_train <- ytelecat[train_indices]
X_train <- Xtelecat[train_indices, ]

y_test  <- ytelecat[-train_indices]
X_test  <- Xtelecat[-train_indices, ]


# Fonction de pr?diction 
predictions <- function(matableTest,muchap,betachap)
{
  ychap <- muchap * rep(1,dim(matableTest)[1])+ as.matrix(matableTest[,]) %*% betachap
  return(ychap)
}

# Fonction de sélection
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

#Question 3 LASSO Bayésien

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
LASSO_BLR <- BLR(
  y = y_train,
  XL = as.matrix(X_train),
  prior = list(
    varE = list(df=2, S=1),
    lambda = list(shape=1, rate=1, type='random')
  ),
  nIter = 3000,
  burnIn = 2000,
  saveAt = "BLR_",
  thin = 1
)

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

#Graphiques
png("3_3_coefficients_LASSO.png", width=800, height=600)
plot(LASSO_BLR$bL, main="Coefficients LASSO bayésien", ylab="beta_j")
dev.off()

png("3_3_coefficients_tries.png", width=800, height=600)
plot(sort(abs(LASSO_BLR$bL)), main="|beta_j| triés")
dev.off()


# remarque : il serait int?ressant de faire varier le shape et le rate de la loi (a priori) gamma du 
# param?tre LASSO lambda. On regarderait alors si la corr?lation est meilleure

#3.3 representations
pdf("representations.pdf", width = 8, height = 6)
plot(LASSO_BLR$bL, main="Coefficients LASSO bayésien")
plot(sort(abs(LASSO_BLR$bL)))
dev.off()

##Question 4

################## Trajectoires simul?es 
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

#Représentation graphique
png("3_4_prior_posterior_varE.png", width=800, height=600)
plot(density(tracevarE), lwd=2, main="Posterior vs Prior varE")
curve(dinvgamma(x,1,1), add=TRUE, lty=2, lwd=2)
legend("topright", c("Posterior","Prior"), lty=c(1,2), lwd=2)
dev.off()


## Question 3.4

tracelambda <- scan('BLR_lambda.dat')
plot(tracelambda,type='o',xlab="iteration",ylab="lambda",main="trace lambda")
plot(density(tracelambda),xlim=c(0,5),col="blue",lwd=2, main="lambda : LASSO bay?sien")
curve(dgamma(x,1,0.1), add=TRUE, col="red",lwd=2,lty=4)
lalegende=c("Posterior","Prior")
legend("topright",lalegende,lwd=c(2,2),lty=c(1,4),col=c("blue","red"))

# graphique
png("3_4_prior_posterior_lambda.png", width=800, height=600)
plot(density(tracelambda), lwd=2, main="Posterior vs Prior lambda")
curve(dgamma(x,1,1), add=TRUE, lty=2, lwd=2)
legend("topright", c("Posterior","Prior"), lty=c(1,2), lwd=2)
dev.off()

png("3_4_trace_varE.png", width=800, height=600)
plot(tracevarE, type='l', main="Trace varE")
dev.off()

png("3_4_trace_lambda.png", width=800, height=600)
plot(tracelambda, type='l', main="Trace lambda")
dev.off()



## Question 3.5

#####################Pr?diction 
predLASSOBLR <- predictions(X_test, LASSO_BLR$mu, LASSO_BLR$bL)
cor(y_test, predLASSOBLR)
#  0.9522485 
plot(predLASSOBLR ~ y_test)
abline(0,1)

#graphique
png("3_5_prediction_test.png", width=800, height=600)
plot(predLASSOBLR ~ y_test, main="Prédictions vs Observations (test)")
abline(0,1, lwd=2)
dev.off()



## Question 3.6

###################### S?lection
plot(LASSO_BLR$bL,xlab="Indice variable",ylab="Param?tre de r?gression associ?", main="LASSO
bay?sien")
plot(sort(abs(LASSO_BLR$bL)))

## Graphique
png("3_6_selection_beta.png", width=800, height=600)
plot(LASSO_BLR$bL, main="Sélection via LASSO")
dev.off()


bb <- boxplot(LASSO_BLR$bL)
bb
varselectedLASSO <- subsetSelected(LASSO_BLR$bL,X_train, bb$stats[1,],bb$stats[5,])
varselectedLASSO
# on va faire du seuillage ? la main pour ?liminer quelques variables 
varselectedLASSO <- subsetSelected(LASSO_BLR$bL,X_train,-0.3,0.3)
varselectedLASSO

# Pour s?lectionner moins de variable il faudrait un srinkage plus fort, 
# ce qui correspond ? un param?tre lambda du LASSO plus ?lev?
# Ici lambda^2 suit une Gamma(e,f) avdc (e,f) = (1,0.1)
# E(lambda^2)=e/f  et V(lamdba^2) = e^2/f
# Prenons e + grand (=10) 








##### Question 4

################# SSVS 
####### Mod?lisation 

### Question 4.2
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
resSSVS <- selection_SSVS(y_train, as.matrix(X_train), 3000, 2000, lec,
                          nbinit, nbToChange, Pi)
resSSVS[[2]]  # on voit combien il y a eu d'?changes dans les 1000 derniers run
resSSVS[[3]]  # on voit le nbre de variables s?lectionn?es lors de ces ?changes 
resSSVS[[1]]

# représentation

png("4_2_SSVS_histogramme.png", width=800, height=600)
plot(sort(abs(resSSVS[[1]])), main="Fréquence de sélection (triée)")
dev.off()

png("4_2_SSVS_boxplot.png", width=800, height=600)
boxplot(resSSVS[[1]], main="Boxplot fréquences sélection SSVS")
dev.off()

png("4_2_SSVS_par_variable.png", width=800, height=600)
plot(resSSVS[[1]],
     xlab="Indice variable",
     ylab="Nombre de sélections post-burn-in",
     main="SSVS : fréquence de sélection par variable")
dev.off()


### Question 4.3

####### S?lection 
plot(sort(abs(resSSVS[[1]])))
boxplot(resSSVS[[1]])
plot(resSSVS[[1]],xlab="Indice variable",ylab= "Nombre de s?lections post-burn-in", main="SSVS")
varselectedSSVS <- subsetSelected(resSSVS[[1]], X_test,0,999)
varselectedSSVS
which(resSSVS[[1]]>400)


# grpahique
png("4_3_SSVS_variables_retenues.png", width=800, height=600)
plot(resSSVS[[1]] > 400,
     main="Variables retenues par SSVS (>400 sélections)",
     ylab="TRUE = sélectionnée")
dev.off()

# 10  20  30  40  50  90 133 147
# On a s?lectionn? les 5 bonnes variables. 
# Les variables suppl?mentairtes ne sont pas celles du LASSO
# Or on sait que le LASSO a tendance ? ne pas s?lectionner des variables 
# corr?l?es. C'est peut ?tre une explication... 

# on peut changer pi 
# en fait le burn-in permet de chercher les bons betas. Une fois qu'ils 
# sont s?l?ctionn?s ils restent car les autres beta ont une densit?s conditionnelle a 
# posteriori moins forte.  

### Question 4.4

Pi2 <- 20/ncol(X_train)

resSSVS_Pi2 <- selection_SSVS(y_train, as.matrix(X_train), 3000, 2000, lec,
                              nbinit, nbToChange, Pi2)

png("4_4_SSVS_effet_Pi.png", width=800, height=600)
plot(resSSVS[[1]], type="l", lwd=2, main="Effet de Pi sur la sélection",
     ylab="Nb sélections")
lines(resSSVS_Pi2[[1]], lty=2, lwd=2)
legend("topright", c("Pi initial","Pi augmenté"), lty=c(1,2), lwd=2)
dev.off()

