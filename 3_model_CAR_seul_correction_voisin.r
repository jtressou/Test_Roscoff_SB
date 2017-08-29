# j'ajoute un commentaire pour essayer
### je chage vraiment
#=================================================================================================
# Modele avec structure CAR propre seul 
# sans variable explicative
#------------------------------------
# parametres du CAR.propre
mynblist = nb2WB(voisinage0) # Convert the neighbourhood object in an object usable by WinBUGS
#
# -------------------
# Parametres du CAR:
# -------------------
# 1. Nb voisin de chaque PARCELLE
nb_voisin <- unlist(mynblist[[3]] )
length(nb_voisin)
head(nb_voisin)
#
# 2. Liste des voisins de chaque parcelle
liste_voisins <- unlist(mynblist[[1]])
#
# 3. Critère de sélection du lien: on supp les liens qui passent sur la N6
cote_parcelle <- parcelle@data$N6
cote_parcelle
# 
# 4. Comparaison du critere de selection "cote_parcelle"
indice <- 0
voisin_a_suppr <- rep(0, length(liste_voisins))
nb_voisin_a_suppr <- rep(0, length(nb_voisin))
#
# a. Sélection de la parcelle  "p"
for (id_parc in 1:(length(nb_voisin))){
  cote_parcelle_parcsel <-  cote_parcelle[id_parc]     #as.character(type$type[p])
  #
  # b. analyse des voisins de la parcelle
  for (v in 1: nb_voisin[id_parc])
  {
  id_v <-  liste_voisins[indice + v]
  #
  cote_parcelle_id_v <- cote_parcelle[id_v]
  
  if (cote_parcelle_parcsel==cote_parcelle_id_v){voisin_a_suppr[indice + v] <- "ok"}else{voisin_a_suppr[indice + v] <- "type_different"
      nb_voisin_a_suppr[id_parc] <-nb_voisin_a_suppr[id_parc]+1} # fin du test de comparaison
  # voisin_a_suppr[indice + v] <- ifelse(cote_parcelle_parcsel==cote_parcelle_id_v, "ok", "type_different")
  # nb_voisin_a_suppr[p] <-ifelse(type_parc==type_id_v, nb_voisin_a_suppr[p], nb_voisin_a_suppr[p]+1) 
 } # fin test sur voisins 
indice <- indice + nb_voisin[id_parc]
} # fin parcelles
#
voisin_a_suppr
nb_voisin_a_suppr
#
# verif
liste_voisins <- liste_voisins[voisin_a_suppr=="ok"]
length(liste_voisins)
nb_voisin <- nb_voisin - nb_voisin_a_suppr
sum(nb_voisin)
tt <- data.frame("nom_parcelle"=parcelle@data$iidt_prf, nb_voisin_a_suppr)
tt
sum(nb_voisin[1:41])
nb_voisin[42]
liste_voisins[(sum(nb_voisin[1:41])):(sum(nb_voisin[1:41])+nb_voisin[42])]
voisin_a_suppr[(sum(nb_voisin[1:41])):(sum(nb_voisin[1:41])+nb_voisin[42])]
tt[210:214,]
#------------


#-------------------------------------
# PARAMETRES DU CAR maj
adj <- liste_voisins
num <- nb_voisin
#-------------------------------------



# 3. Poids des voisins de chaque unité i = 1/n_i
Ci <- NULL # creation du vecteur C prenant 1/n_i si voisin, 0 sinon
for(i in 1:length(num)){tt <- rep((1/ num[i]), num[i])
Ci <- c(Ci, tt)}
Ci  
#=======================

# diagonale de la matrice Mii of the conditional variance matrix 
M <- rep(sigma2, N)
#=================================================================================================









voir comment mettre le sigma2 utilisé pour construire M[]

 

datalist=c(Y=Y, N=N, C[]=Ci, adj[], num[], M[], prec, rho, sigma2)
#
init1=list( u=rnorm(N,0,4), mu=runif(1,-50,50),sigma=0.1, tau=0.1)
init2=list( u=rnorm(N,0,4), mu=runif(1,-50,50),sigma=0.1, tau=0.1)
inits = list(init1,init2)

genDataFile(datalist,"cancerdata.txt")
genInitsFile(2, inits, "cancer.init")

out = bugs(datalist, inits=list(init1,init2), parameters.to.save=c('mu', 'sigma2', 'tau2', 'u'), model.file="model0.txt",
           n.chains=2, n.iter=10000, n.burnin=1000, debug=FALSE, codaPkg=FALSE,n.thin=10,
           #bugs.directory='C://Program Files (x86)//WinBUGS14',   #'C:/MATHS/WINBUGS14',program="WinBUGS",
           bugs.directory='C://Program Files (x86)//OpenBUGS//OpenBUGS323',
           program="OpenBUGS",DIC=FALSE,working.directory=getwd())





# Likelihood
for (p in 1:nbP) {
  y[p] ~ dpois(lambda[p])
  log(lambda[p]) <-  mu + P[p]
}
#
# SPATIAL EFFECT AT PARCEL LEVEL
P[1:nbP] ~ car.proper(0, C[], adj[], num[], M[], prec, rho)
variance_spP <-1/prec
#
# PRIORS
mu ~ dnorm(0,0.001)
prec  ~ dgamma(0.5, 0.0005)
# standard deviation
rho.min <- min.bound(C[], adj[], num[], M[])
rho.max <- max.bound(C[], adj[], num[], M[])
rho ~ dunif(rho.min, rho.max)





DATA:
  data <- list ("y", "mu_S", "C", "adj", "num", "M")

y<- obs
mu_S[] <- A vector giving the mean for each area
C[]   <- A vector the same length as adj[] giving normalised weights associated with each pair of  areas 
adj[] <- ID numbers of the adjacent areas
num[] <- A vector of length N (the total number of areas) giving the number of neighbours  ni for each area.
M[] <-  A vector of length N giving the diagonal elements Mii of the conditional variance matrix 



tau <-  A scalar parameter representing the overall precision (in verse variance) parameter.
rho <-  A scalar parameter representing the overall degree of spatial dependence. This parameter is
constrained to lie between bounds given by the inverse of the minimum and maximum eigenvalues of the matrix


spatial.exp: 
  mu[]   A vector giving the mean for each area 
(this can either be entered as data, assigned a prior distribution, or specified deterministically within the model code).

x[] and y[] : Vectors of length N giving the x and y coordinates of the location of each point, or the centroid of
each area

tau :  A scalar parameter representing the overall precision (inverse variance) parameter.

phi : A scalar parameter representing the rate of decline of correlation with distance between points.
Note that the magnitude of this parameter will depend on the units in wh ich the x and y coordinates of each location are
measured (e.g. metres, km etc.).
kappa: A scalar parameter controlling the amount of spatial smoothing. This is constrained to lie in the interval
[0, 2).




data <- list ("y", "y", "sigma.y")

inits <- function(){
  list(theta = rnorm(J, 0, 100), mu.theta = rnorm(1, 0, 100),
       sigma.theta = runif(1, 0, 100))
}



