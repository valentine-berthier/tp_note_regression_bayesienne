getwd()
setwd("/home/id2467/tp_note_regression_bayesienne/Cours + TP 2024-2025-20260114/TP_note")
getwd()

telecat <- read.csv(file = "telecat.csv")

set.seed(1234)

library(MCMCpack)
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
mean_tau2 <- mean(LASSO_BLR$tau2)
mean_tau2


png("3_3_importance_beta.png", width=1000, height=600)

beta <- as.vector(LASSO_BLR$bL)
names(beta) <- colnames(X_train)
beta_sorted <- sort(abs(beta), decreasing = TRUE)

barplot(beta_sorted[1:20],
        main="20 plus grands |beta| - LASSO bayésien",
        las=2)

dev.off()



# remarque : il serait int?ressant de faire varier le shape et le rate de la loi (a priori) gamma du 
# param?tre LASSO lambda. On regarderait alors si la corr?lation est meilleure




################## Trajectoires simul?es 
# Attention : 
# si on relance alors les trajecgoires se cumulent (concat?nation des r?sultats)
# il faut soit effac? BLR_, soit changer le nom...
# Par ailleurs, on ne garde qu'une valeur sur 10 (par d?faut, pour garantir une "ind?pendance") 


## Question 3.4

tracevarE   <- scan("BLR_varE.dat")
tracelambda <- scan("BLR_lambda.dat")

varE = list(df=2, S=1)

lambda = list(shape=1, rate=1)

png("3_4_traceplots.png", width=1000, height=500)

par(mfrow=c(1,2))

plot(tracevarE, type='l',
     main="Trace plot varE",
     ylab="varE", xlab="Iteration")

plot(tracelambda, type='l',
     main="Trace plot lambda",
     ylab="lambda", xlab="Iteration")

dev.off()

png("3_4_prior_vs_posterior.png", width=1000, height=500)

par(mfrow=c(1,2))

# varE
plot(density(tracevarE),
     main="Posterior vs Prior : varE",
     lwd=2)
curve(dinvgamma(x, shape=1, scale=1),
      add=TRUE, lty=2, lwd=2)
legend("topright", c("Posterior","Prior"), lty=c(1,2), lwd=2)

# lambda
plot(density(tracelambda),
     main="Posterior vs Prior : lambda",
     lwd=2)
curve(dgamma(x,1,1),
      add=TRUE, lty=2, lwd=2)
legend("topright", c("Posterior","Prior"), lty=c(1,2), lwd=2)

dev.off()




## Question 3.5

#####################Pr?diction 
predLASSOBLR <- predictions(X_test, LASSO_BLR$mu, LASSO_BLR$bL)
cor(y_test, predLASSOBLR)
#  0.9522485 
png("3_5_prediction_test.png", width=800, height=600)
plot(predLASSOBLR ~ y_test, main="Prédictions vs Observations (test)")
abline(0,1, lwd=2)
dev.off()



## Question 3.6

###################### S?lection

# Récupération des coefficients
beta <- as.vector(LASSO_BLR$bL)
names(beta) <- colnames(X_train)
bb <- boxplot(beta, plot = FALSE)
seuil_bas  <- bb$stats[1]
seuil_haut <- bb$stats[5]

varselectedLASSO <- subsetSelected(beta, X_train, seuil_bas, seuil_haut)
varselectedLASSO
selec <- names(beta) %in% varselectedLASSO$selected
cols <- ifelse(selec, "blue", "black")

# Graphique
png("3_6_selection_LASSO.png", width=900, height=600)

plot(beta,
     col = cols,
     pch = 16,
     main = "Sélection des variables par le LASSO bayésien",
     xlab = "Indice des variables",
     ylab = "Valeur du coefficient beta")

legend("topright",
       legend = c("Variable sélectionnée", "Variable non sélectionnée"),
       col = c("blue", "black"),
       pch = 16)

dev.off()

sort(abs(LASSO_BLR$bL), decreasing = TRUE)[1:5]

beta <- as.vector(LASSO_BLR$bL)
names(beta) <- colnames(X_train)

beta[order(abs(beta), decreasing = TRUE)][1:5]

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

# Seuil de sélection (comme tu faisais : > 400)
seuil <- 400
selec <- resSSVS[[1]] > seuil

cols <- ifelse(selec, "blue", "black")

png("4_2_SSVS_par_variable.png", width=900, height=600)

plot(resSSVS[[1]],
     col = cols,
     pch = 16,
     xlab = "Indice variable",
     ylab = "Nombre de sélections post-burn-in",
     main = "SSVS : fréquence de sélection par variable")

legend("topright",
       legend = c("Variable sélectionnée", "Variable non sélectionnée"),
       col = c("blue", "black"),
       pch = 16)

dev.off()


which(resSSVS[[1]] > 400)


### Question 4.3

####### S?lection 
### Question 4.3 – Sélection SSVS

# 4_3_SSVS_triés.png : fréquences de sélection triées
png("4_3_SSVS_triés.png", width=800, height=600)
plot(sort(abs(resSSVS[[1]])), main="SSVS : Fréquence de sélection (triée)", ylab="Nombre de sélections post-burn-in")
dev.off()

# 4_3_SSVS_boxplot.png : boxplot des fréquences
png("4_3_SSVS_boxplot.png", width=800, height=600)
boxplot(resSSVS[[1]], main="SSVS : Boxplot des fréquences de sélection")
dev.off()

# 4_3_SSVS_par_variable.png : fréquence de sélection par variable
png("4_3_SSVS_par_variable.png", width=800, height=600)
plot(resSSVS[[1]], xlab="Indice variable", ylab="Nombre de sélections post-burn-in",
     main="SSVS : fréquence de sélection par variable")
dev.off()

# Sélection des variables avec seuil > 400
varselectedSSVS <- subsetSelected(resSSVS[[1]], X_test,0,999)
varselectedSSVS
which(resSSVS[[1]]>400)


# grpahique
png("4_3_SSVS_variables_retenues.png", width=800, height=600)

# sélection avec seuil >400
selec <- resSSVS[[1]] > 400

# couleur : bleu si sélectionnée, noir sinon
cols <- ifelse(selec, "blue", "black")

plot(1:length(selec), selec*1, col=cols, pch=16, ylim=c(0,1),
     main="Variables retenues par SSVS (>400 sélections)",
     xlab="Indice variable", ylab="Sélectionnée (1=oui)")
legend("topright", legend=c("Sélectionnée","Non sélectionnée"), col=c("blue","black"), pch=16)

dev.off()

names(X_train)[resSSVS[[1]] > 400]



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

