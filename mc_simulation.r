##############################################
###############MCMC###########################
##############################################
#install.packages("devtools")
library(devtools)
#install.packages("moments")
library(moments)

#source_url("https://github.com/wendtvi/dissertation_code/blob/main/continous_nonlinear.r?raw=TRUE")
N=1000
matriz_resultados=mc_function(N=1000)
vetor_mc_resultados_vies=vector()
vetor_mc_resultados_sd=vector()
vetor_mc_resultados_media_pop=vector()
vetor_mc_resultados_media_est=vector()
vetor_mc_resultados_sd_vies=vector()

xvetor_mc_resultados_vies=vector()
xvetor_mc_resultados_sd=vector()
xvetor_mc_resultados_media_pop=vector()
xvetor_mc_resultados_media_est=vector()
xvetor_mc_resultados_sd_vies=vector()

yvetor_mc_resultados_vies=vector()
yvetor_mc_resultados_sd=vector()
yvetor_mc_resultados_media_pop=vector()
yvetor_mc_resultados_media_est=vector()
yvetor_mc_resultados_sd_vies=vector()
c=0
for (k in seq(1,ncol(matriz_resultados)/3,2)){
  c=c+1
  vetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
  vetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
  vetor_mc_resultados_sd[c]=var(matriz_resultados[,k+1])*sqrt((kurtosis(matriz_resultados[,k+1])-1)/N)
  vetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
  vetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
}

c=0
for (k in seq(ncol(matriz_resultados)/3+1,ncol(matriz_resultados)*2/3,2)){
  c=c+1
  xvetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
  xvetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
  xvetor_mc_resultados_sd[c]=var(matriz_resultados[,k+1])*sqrt((kurtosis(matriz_resultados[,k+1])-1)/N)
  xvetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
  xvetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
}


c=0
for (k in seq(ncol(matriz_resultados)*2/3+1,ncol(matriz_resultados),2)){
  c=c+1
  yvetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
  yvetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
  yvetor_mc_resultados_sd[c]=var(matriz_resultados[,k+1])*sqrt((kurtosis(matriz_resultados[,k+1])-1)/N)
  yvetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
  yvetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
}


par(mfrow=c(1,2))
hist(matriz_resultados[,8],main="Distribuição estimador Wooldridge t=4 e sd 2 no modelo DiD")
hist(matriz_resultados[,8],main="Distribuição estimador CIC t=4 e sd 2 no modelo DiD")
