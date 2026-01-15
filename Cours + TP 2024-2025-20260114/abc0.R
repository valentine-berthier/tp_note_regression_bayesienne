# ABC pour une loi normale
n = 100
data = rnorm(n, mean = 2, sd = 2)

#On veut estimer la moyenne µ et la variance σ
#2 par la méthode ABC.
#Nous n'avons pas d'a priori précis sur µ et σ et nous mettons des lois uniformes (on peut
#faire varier les hyperparamètres) :
simMU <- function () 
{
return (runif(1, min=-10, max=10)) 
}

simSIGMA = function () 
{
return (runif(1, min=0, max=10)) 
}
#On simule des Y candidats à ressembler aux observations :
simul = function (n, mu, sigma) 
{
return(rnorm(n, mean = mu, sd = sigma)) 
}
#Il faut maintenant comparer les simulations aux observations et accepter ou non l'égalité.
#Pour cela il faut résumer les données et comparer les "summary statistics". On peut par
#exemple
# Faire un test d'égalité des lois des données.
# Comparer des fractiles des données.
# Compare les moyennes et les variances des données.

# I. Distance basée sur la moyenne et la variance :
MeanVar = function (obs, simu) 
{
diffmean <- abs(mean(obs) - mean(simu))
diffsd <- abs(sd(obs) - sd(simu))
return(c(diffmean, diffsd)) 
}
#La fonction d'acceptation-rejet va être de la forme (1 si accepté, 0 sinon) :
acceptMeanVar = function (obs, simu, seuil) 
{
differences = MeanVar(obs, simu)
if((differences[1] < seuil[1]) & (differences[2] < seuil[2])) return(1) else return(0)
}

#Test : 
obs0=rnorm(100)
simu0=rnorm(100)
MeanVar(obs0,simu0)
acceptMeanVar(obs0,simu0,c(0.,0.1))

# Programme :
sampledMeanVar = function (obs, nIter, seuil) 
{
n = length(obs)
decision = vector(length = nIter) # rep(0,nIter)
sampledMU = vector(length = nIter, mode = "numeric")
sampledSIGMA = vector (length = nIter, mode = "numeric")
for (i in 1 :nIter)
{
mu <- simMU()
sigma <- simSIGMA()
#parameters = list("mu"=mu, "sigma"=sigma) # A QUOI SERT CETTE LIGNE ?
simuDATA = simul(n, mu, sigma)
decision[i] = acceptMeanVar(obs, simuDATA, seuil)
sampledMU[i] = mu
sampledSIGMA[i] = sigma 
}
return(data.frame(cbind('DECISION' = decision, 'sampledMU' = sampledMU, 'sampledSIGMA' = sampledSIGMA))) 
}

#Mise en pratique :
resuMeanVar = sampledMeanVar(data, 500000, c(0.1,0.1))
sum(resuMeanVar$DECISION)
ind=which(resuMeanVar$DECISION==1)
ind
resumu=resuMeanVar$sampledMU[ind]
resuSIGMA=resuMeanVar$sampledSIGMA[ind]
plot(resumu)
plot(density(resumu))
mean(resumu)
plot(resuSIGMA)
plot(density(resuSIGMA))
mean(resuSIGMA)


# II. Distance basée sur les quantiles :
quant = function(data) 
{
return (quantile(data, probs=c(0.1, 0.5, 0.9))) 
}
comparequant = function (obs, simu) 
{
compare = sqrt(sum(mapply(function(x,y) (x-y)**2, obs, simu)))
return(compare) 
}

#La fonction d'acceptation-rejet va être de la forme (1 si accepté, 0 sinon) :
acceptQuant = function (obs, simu, seuil) 
{
distance = comparequant(quant(obs), quant(simu))
if (distance < seuil) return(1) else return(0)
}

#Test : 
obs0=rnorm(100)
simu0=rnorm(100)
quant(obs0)
quant(simu0)
comparequant(quant(obs0),quant(simu0))
acceptQuant(obs0,simu0,1)

# Programme : 
sampledQuant = function (obs, nIter, seuil) 
{
n = length(obs)
decision = vector(length = nIter)
sampledMU = vector(length = nIter, mode = "numeric")
sampledSIGMA = vector (length = nIter, mode = "numeric")
for (i in 1 :nIter)
{
mu <- simMU()
sigma <- simSIGMA()
# parameters = list("mu"=mu, "sigma"=sigma) # A QUOI SERT CETTE LIGNE ?
simuDATA = simul(n, mu, sigma)
decision[i] = acceptQuant(obs, simuDATA, seuil)
sampledMU[i] = mu
sampledSIGMA[i] = sigma 
}
return(data.frame(cbind('DECISION' = decision, 'sampledMU' = sampledMU, 'sampledSIGMA' = sampledSIGMA))) 
}

#Mise en pratique :
resuQuant = sampledQuant(data, 100000, 0.5)
sum(resuQuant$DECISION)
ind=which(resuQuant$DECISION==1)
ind
resumu=resuQuant$sampledMU[ind]
resuSIGMA=resuQuant$sampledSIGMA[ind]
plot(resumu)
plot(density(resumu))
mean(resumu)
plot(resuSIGMA)
plot(density(resuSIGMA))
mean(resuSIGMA)

# III. Programme finalglobal :
sampled = function (obs, nIter, seuil, acceptfunction) 
{
n = length(obs)
decision = vector(length = nIter)
sampledMU = vector(length = nIter, mode = "numeric")
sampledSIGMA = vector (length = nIter, mode = "numeric")
for (i in 1 :nIter)
{
mu <- simMU()
sigma <- simSIGMA()
parameters = list("mu"=mu, "sigma"=sigma) # A QUOI SERT CETTE LIGNE ?
simuDATA = simul(n, mu, sigma)
decision[i] = acceptfunction(obs, simuDATA, seuil)
sampledMU[i] = mu
sampledSIGMA[i] = sigma 
}
return(data.frame(cbind('DECISION' = decision, 'sampledMU' = sampledMU, 'sampledSIGMA' = sampledSIGMA))) 
}

#Mise en pratique :
resu = sampled(data, 100000, c(2,2), acceptMeanVar)
#Pourcentage d'acceptation :
sum(resu$decision)/n
#Trace et densités des paramètres :
resMU= res$sampledMU
resSIGMA2=res$sampledSIGM**2
plot(resMU)
plot(density(resMU))
plot(resSIGMA2)
plot(density(resSIGMA2))

Conclusion : 
L'utilisation de la moyenne et de la variance comme critère marche bien ici, 
puisque l'on sait que la loi normale est caractérisée par ces deux paramètres. 
Mais en général, il vaut mieux avoir le plus de critère possible, 
quitte à faire un test d'adéquation. 

On peut ajouter des covariables (fixed design) et se placer dans un cadre de régression. 

On pourrait s'intéresser à des mélanges gaussiens. Refaire l'ABC mais avec k mélanges? 
On simulerait alors les proportions et les paramètres des lois (moyennes et variance). 
Ce serait évidemment un plus long (en fonction de k) à converger... 

Pour les GMM on peut aussi faire du bayésien avec Gibbs, ou Metropolis Hasting.  Et on peut 
ajouter des covariables. 



