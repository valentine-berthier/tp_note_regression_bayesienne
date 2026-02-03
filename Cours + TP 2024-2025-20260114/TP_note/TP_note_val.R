
# ---
# title: Code du TP Noté de
# author: "Valentine BERTHIER & Maïlis LANNE"
# date: "2026-01-15"
# ---


rm(list=ls())
# à modifier lors du lancement du code
setwd(paste0(getwd(),'/Cours + TP 2024-2025-20260114/TP_note'))
source("fonctions.r")


library(BGLR)      # Pour BayesA, BL, BRR
library(rrBLUP)    # Pour RR-BLUP
library(mvtnorm)   # Pour simulations multivariées
library(mnormt)    # Pour loi normale multivariée
library(invgamma)  # Pour loi inverse gamma
library(glmnet)    # Pour LASSO classique (question 5.6)


# Mise en place ----

# Chargement des données
data <- read.csv("telecat.csv")


# Analyse descriptive ----
colnames(data)
summary(data)
skim(data)


# On définit une clé pour assurer la reproductibilité du code
set.seed(12345)


# Séparation training/test ----
n <- nrow(data)
train_idx <- sample(1:n, 100)
test_idx <- setdiff(1:n, train_idx)

y_train <- data$Y[train_idx]
y_test <- data$Y[test_idx]
X_train <- as.matrix(data[train_idx, -which(names(data) == "Y")])
X_test <- as.matrix(data[test_idx, -which(names(data) == "Y")])




# Partie 1 ----
## 1.1 Estimations ----
# Le modèle: Y = mu*1 + X*beta + epsilon
# beta ~ N(0, sigma2_u * I) - variance commune pour tous les beta
# epsilon ~ N(0, sigma2_e)
fit_RR <- mixed.solve(y_train, 
                      Z = X_train,  # Matrice design
                      X = matrix(1, nrow=100, ncol=1),  # Intercept
                      method = "REML")

# Résultats
mu_RR <- fit_RR$beta        # Estimation de mu (intercept)
beta_RR <- fit_RR$u         # Estimations des coefficients
sigma2_e_RR <- fit_RR$Ve    # Variance des epsilon
sigma2_u_RR <- fit_RR$Vu    # Variance des beta

cat("\nmu estimé:", mu_RR, "\n")
cat("sigma2_epsilon:", sigma2_e_RR, "\n")
cat("sigma2_beta:", sigma2_u_RR, "\n")

# Visualisation des coefficients
par(mfrow=c(1,2))
plot(beta_RR, main="Coefficients RR-BLUP", 
     xlab="Indice variable", ylab="Coefficient")
abline(h=0, col="red", lty=2)
plot(sort(abs(beta_RR)), main="Coefficients triés (valeur absolue)", 
     xlab="Rang", ylab="|Coefficient|")

# Observation: FORT EFFET SHRINKAGE - tous les coefficients sont "rétrécis" vers 0
# C'est l'effet du prior gaussien avec variance commune

# à commenter

## 1.2 Prédiction ----

y_pred_RR <- predictions(X_test, mu_RR, beta_RR)
cor_RR <- cor(y_test, y_pred_RR)
cat("\n1.2 Corrélation RR-BLUP (test):", round(cor_RR, 4), "\n")

# Interprétation: Cette corrélation mesure la qualité prédictive du modèle
# Plus elle est élevée, meilleur est le modèle

par(mfrow=c(1,1))
plot(y_test, y_pred_RR, 
     main=paste("RR-BLUP: Corrélation =", round(cor_RR, 3)),
     xlab="Y observé (test)", ylab="Y prédit")
abline(0, 1, col="red", lwd=2)



# à commenter

## 1.3 Sected_RR ----
threshold <- quantile(abs(fit_RR$ETA[[1]]$b), 0.95)
# Le seuil correspond à la valeur à 95%, c'est à dire la "moustache" supérieure
selected_vars_RR <- which(abs(fit_RR$ETA[[1]]$b) > threshold)

rownames(data.frame(selected_vars_RR))

library(dplyr)
library(tibble)


# install.packages("ggplot2")
library(ggplot2)

coeff <- data.frame(fit_RR$ETA[[1]]$b)
colnames(coeff)<-'y'

coeff$variable <- rownames(coeff)
# Créer la colonne selected_var
vars <- rownames(data.frame(selected_vars_RR))
coeff$selected_var <- ifelse(coeff$variable %in% vars, "oui", "non")



bp_RR <- boxplot(beta_RR, plot=FALSE)
selected_RR <- subsetSelected(beta_RR, bp_RR$stats[1,], bp_RR$stats[5,])

cat("\n1.3 Variables sélectionnées avec RR-BLUP:", length(selected_RR$indices), "\n")
print(data.frame(
  Variable = selected_RR$noms,
  Indice = selected_RR$indices,
  Coefficient = round(selected_RR$valeurs, 4)
))




# Basic box plot
ggplot(coeff, aes(x = "", y = y)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    aes(color = selected_var),
    width = 0.15,
    size = 2
  ) +
  scale_color_manual(
    values = c("oui" = "red", "non" = "black")
  )




# Partie 2 ----
 
## 2.1 Question ----  
library(MCMCpack) 
library(LearnBayes)
library(invgamma)
library(rmnorm)

# La différence : Bayes A permet une variance spécifique pour chaque coefficient 
# (hétéroscédasticité), tandis que RR-BLUP utilise une variance commune.


## 2.2 Estimations ---- 

# Choix des hyper-paramètres
# Pour des priors vagues (non informatifs):
priora <- 1
priorb <- 1
priorc <- 1
priord <- 1

cat("\nHyper-paramètres choisis:\n")
cat("- a =", priora, ", b =", priorb, "(pour sigma2_beta_j ~ InvGamma(a,b))\n")
cat("- c =", priorc, ", d =", priord, "(pour sigma2_epsilon ~ InvGamma(c,d))\n")
cat("\nAvec a=b=c=d=1:\n")
cat("- E[sigma2] = b/(a-1) = INFINI (prior très vague)\n")
cat("- Var[sigma2] = INFINI\n")
cat("=> On laisse les DONNÉES décider des valeurs\n")

# Paramètres MCMC
nbiter_BA <- 5000
nburn_BA <- 3000  # 60% de burn-in

cat("\nParamètres MCMC:\n")
cat("- Nombre d'itérations:", nbiter_BA, "\n")
cat("- Burn-in:", nburn_BA, "\n")
cat("- Itérations utilisées:", nbiter_BA - nburn_BA, "\n")

# Application sur les données d'entraînement
cat("\nLancement de l'algorithme de Gibbs...\n")
res_BayesA <- BayesA(
  y = y_train, 
  X = X_train, 
  a = priora, 
  b = priorb, 
  c = priorc, 
  d = priord, 
  muinit = mean(y_train),  # Initialisation de mu avec la moyenne empirique
  nbiter = nbiter_BA, 
  nburn = nburn_BA
)

# Extraction des résultats
resbeta_BA <- res_BayesA[[1]]
ressigma2beta_BA <- res_BayesA[[2]]
resmu_BA <- res_BayesA[[3]]
ressigma2eps_BA <- res_BayesA[[4]]

# Estimations = moyennes a posteriori
mu_BA <- mean(resmu_BA)
beta_BA <- rowMeans(resbeta_BA)
sigma2_e_BA <- mean(ressigma2eps_BA)
sigma2_beta_BA <- rowMeans(ressigma2beta_BA)


cat("\n========== RÉSULTATS BAYES A ==========\n")
cat("mu estimé:", round(mu_BA, 4), "\n")
cat("sigma2_epsilon estimé:", round(sigma2_e_BA, 4), "\n")
cat("Statistiques sur sigma2_beta_j:\n")
cat("  - Médiane:", round(median(sigma2_beta_BA), 4), "\n")
cat("  - Min:", round(min(sigma2_beta_BA), 4), "\n")
cat("  - Max:", round(max(sigma2_beta_BA), 4), "\n")

# Visualisation des estimations
par(mfrow=c(2,2))
plot(beta_BA, main="Coefficients beta estimés (Bayes A)",
     xlab="Indice variable", ylab="beta_j", pch=20)
abline(h=0, col="red", lty=2)

plot(sort(abs(beta_BA)), main="Coefficients triés (valeur absolue)",
     xlab="Rang", ylab="|beta_j|", type="l", lwd=2, col="darkgreen")

plot(sigma2_beta_BA, main="Variances sigma2_beta_j estimées",
     xlab="Indice variable", ylab="sigma2_beta_j", pch=20, col="blue")

boxplot(beta_BA, main="Boxplot des coefficients", ylab="beta_j")




moybeta <- apply(res_BayesA[[1]], 1, mean) # estimateurs des beta
moysigma2beta <- apply(res_BayesA[[2]], 1, mean)
moymu <- mean(res_BayesA[[3]])
moysigma2eps <- mean(res_BayesA[[4]])
moybeta[50]
plot(moybeta)
plot(sort(abs(moybeta)))
moymu
plot(moysigma2beta)
moysigma2eps

# regardons les trajectoires : 

png('test.png')
par(mfrow=c(2,2))
plot(res_BayesA[[1]][50, ])  # pour beta50
plot(res_BayesA[[2]][1, ])   # pour s2beta1
plot(res_BayesA[[3]])        # pour mu
plot(res_BayesA[[4]])        # pour s2eps 
dev.off()


## 2.3 Burn-in plus grande ----

cat("\n========== 2.3 POURQUOI UN GRAND BURN-IN? ==========\n")
cat("Le burn-in est CRUCIAL pour la validité des résultats:\n")
cat("\n1. PROBLÈME D'INITIALISATION:\n")
cat("   - Au départ: paramètres initialisés ARBITRAIREMENT\n")
cat("   - Ex: beta = 0, sigma2_beta ~ InvGamma(a,b) aléatoire\n")
cat("   - Ces valeurs NE SUIVENT PAS la distribution a posteriori\n")

cat("\n2. PHASE DE CONVERGENCE:\n")
cat("   - Les premières itérations servent à 'explorer' l'espace\n")
cat("   - La chaîne se déplace vers les zones de HAUTE DENSITÉ a posteriori\n")
cat("   - Il faut du temps pour 'oublier' les valeurs initiales\n")

cat("\n3. CONSÉQUENCES D'UN BURN-IN TROP PETIT:\n")
cat("   - BIAIS dans les estimations (influence de l'initialisation)\n")
cat("   - Les moyennes a posteriori sont FAUSSÉES\n")
cat("   - Les intervalles de crédibilité sont INCORRECTS\n")

cat("\n4. RÈGLE EMPIRIQUE:\n")
cat("   - Burn-in = 50-80% du nombre total d'itérations\n")
cat("   - Ici:", nburn_BA, "/", nbiter_BA, "=", 
    round(100*nburn_BA/nbiter_BA, 1), "%\n")
cat("   - Pour être sûr: regarder les TRACES (question 2.4)\n")

cat("\n5. ALTERNATIVE: CRITÈRES DE CONVERGENCE:\n")
cat("   - Diagnostic de Gelman-Rubin (comparer plusieurs chaînes)\n")
cat("   - Critère de Geweke (comparer début et fin de chaîne)\n")
cat("   - Inspection visuelle des traces\n")


## 2.4 Trace après burn-in ----


cat("\n========== 2.4 EXAMEN DES TRACES ==========\n")
cat("Analysons les trajectoires de quelques paramètres...\n")

# On prend les 500 premières itérations après burn-in
n_trace <- min(500, ncol(resbeta_BA))

par(mfrow=c(2,3))

# Traces de quelques beta
plot(resbeta_BA[10, 1:n_trace], type="l", 
     main="Trace de beta_10", xlab="Itération", ylab="beta_10")
abline(h=mean(resbeta_BA[10,]), col="red", lwd=2)

plot(resbeta_BA[50, 1:n_trace], type="l",
     main="Trace de beta_50", xlab="Itération", ylab="beta_50")
abline(h=mean(resbeta_BA[50,]), col="red", lwd=2)

# Trace de mu
plot(resmu_BA[1:n_trace], type="l",
     main="Trace de mu", xlab="Itération", ylab="mu", col="blue")
abline(h=mean(resmu_BA), col="red", lwd=2)

# Trace de sigma2_eps
plot(ressigma2eps_BA[1:n_trace], type="l",
     main="Trace de sigma2_epsilon", xlab="Itération", ylab="sigma2_eps", col="darkgreen")
abline(h=mean(ressigma2eps_BA), col="red", lwd=2)

# Autocorrélogramme de beta_50
acf(resbeta_BA[50,], main="ACF de beta_50", lag.max=50)

# Densité a posteriori de beta_50
plot(density(resbeta_BA[50,]), main="Densité a posteriori de beta_50",
     xlab="beta_50", col="purple", lwd=2)
abline(v=mean(resbeta_BA[50,]), col="red", lwd=2, lty=2)

cat("\n========== INTERPRÉTATION DES TRAJECTOIRES ==========\n")
cat("\nBONNE CONVERGENCE (ce qu'on veut voir):\n")
cat("- Trajectoire STABLE, oscillant autour d'une valeur moyenne\n")
cat("- PAS de tendance montante ou descendante\n")
cat("- PAS de rupture brutale ou de changement de régime\n")
cat("- Ressemble à du 'BRUIT BLANC' autour d'une moyenne\n")

cat("\nMAUVAISE CONVERGENCE (problèmes):\n")
cat("- Tendance claire (montée/descente) => burn-in insuffisant\n")
cat("- Ruptures, sauts => problème d'exploration\n")
cat("- Variance qui change => non-stationnarité\n")

cat("\nAUTOCORRÉLOGRAMME (ACF):\n")
cat("- Mesure la CORRÉLATION entre valeurs successives\n")
cat("- Idéal: décroissance RAPIDE vers 0\n")
cat("- Si décroissance lente: forte dépendance temporelle\n")
cat("  => Solution: prendre 1 valeur sur 10 (thinning)\n")

cat("\nDENSITÉ A POSTERIORI:\n")
cat("- Forme régulière (souvent proche de normale)\n")
cat("- Symétrie ou légère asymétrie acceptable\n")
cat("- Multi-modalité = problème potentiel\n")



## 2.5 Choix Hypp ----
# 2.5 Choix d'hyper-paramètres ayant peu d'influence
cat("\n========== 2.5 HYPER-PARAMÈTRES À FAIBLE INFLUENCE ==========\n")
cat("OBJECTIF: Choisir des priors VAGUES (non informatifs)\n")

cat("\nCHOIX PROPOSÉ: a = b = c = d = 1\n")
cat("\nPOUR LES VARIANCES sigma2_beta_j ~ InvGamma(a, b):\n")
cat("- Espérance: E[sigma2] = b/(a-1) = 1/0 = INFINI\n")
cat("- Variance: Var[sigma2] = b²/((a-1)²(a-2)) = INFINI\n")
cat("- Interprétation: AUCUNE contrainte a priori sur la variance\n")

cat("\nPOUR LA VARIANCE sigma2_epsilon ~ InvGamma(c, d):\n")
cat("- Même raisonnement: espérance et variance infinies\n")
cat("- Prior très DIFFUS, très étalé\n")

cat("\nPOURQUOI CES VALEURS?\n")
cat("1. On veut que les DONNÉES 'parlent' (soient informatives)\n")
cat("2. On n'impose PAS de contraintes a priori fortes\n")
cat("3. Le prior ne doit pas dominer la vraisemblance\n")
cat("4. Maximum de FLEXIBILITÉ pour le modèle\n")

cat("\nALTERNATIVE (si on veut un prior légèrement informatif):\n")
cat("- a = c = 2, b = d = 1\n")
cat("- E[sigma2] = 1/(2-1) = 1 (prior centré sur 1)\n")
cat("- Var[sigma2] = 1 (variance finie mais grande)\n")
cat("- Utile si on a une idée approximative de l'échelle\n")


## 2.6 Corrélations des prédictions----
# A l'aide des estimations obtenues, donner la corrélation des prédictions avec les
# observations du jeu de validation et comparer cette corrélation à celle précédente du
# Random Regression.

y_pred_BA <- predictions(X_test, mu_BA, beta_BA)
cor_BA <- cor(y_test, y_pred_BA)

cat("\n2.6 COMPARAISON DES CORRÉLATIONS:\n")
cat("RR-BLUP:", round(cor_RR, 4), "\n")
cat("Bayes A:", round(cor_BA, 4), "\n")
cat("Différence:", round(cor_BA - cor_RR, 4), "\n")

if(cor_BA > cor_RR) {
  cat("=> Bayes A MEILLEUR: les variances spécifiques améliorent la prédiction\n")
} else {
  cat("=> RR-BLUP comparable ou meilleur (peut dépendre des données)\n")
}

par(mfrow=c(1,2))
plot(y_test, y_pred_RR, main=paste("RR-BLUP: r =", round(cor_RR, 3)),
     xlab="Y observé", ylab="Y prédit")
abline(0, 1, col="red", lwd=2)
plot(y_test, y_pred_BA, main=paste("Bayes A: r =", round(cor_BA, 3)),
     xlab="Y observé", ylab="Y prédit")
abline(0, 1, col="red", lwd=2)



## 2.7 selected_BA ----
# A l'aide des estimations obtenues, sélectionner les variables explicatives qui paraissent
# les plus pertinentes. On pourra comparer les meilleures variables retenues ici par rapport 
# aux meilleures retenues avec RR-BLUP.


bp_BA <- boxplot(beta_BA, plot=FALSE)
selected_BA <- subsetSelected(beta_BA, bp_BA$stats[1,], bp_BA$stats[5,])

cat("\n2.7 Variables sélectionnées avec Bayes A:", length(selected_BA$indices), "\n")
print(data.frame(
  Variable = selected_BA$noms,
  Indice = selected_BA$indices,
  Coefficient = round(selected_BA$valeurs, 4)
))

cat("\nCOMPARAISON RR-BLUP vs BAYES A:\n")
common_vars <- intersect(selected_vars_RR$indices, selected_BA$indices)
cat("Variables communes:", length(common_vars), "sur", 
    max(length(selected_vars_RR$indices), length(selected_BA$indices)), "\n")




# Partie 3 ----






















