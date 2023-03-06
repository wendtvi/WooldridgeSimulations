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
WOOLvetor_mc_resultados_vies=vector()
WOOLvetor_mc_resultados_sd=vector()
WOOLvetor_mc_resultados_media_pop=vector()
WOOLvetor_mc_resultados_media_est=vector()
WOOLvetor_mc_resultados_sd_vies=vector()
WOOLvetor_mc_resultados_mcsd=vector()

WOOLxvetor_mc_resultados_vies=vector()
WOOLxvetor_mc_resultados_sd=vector()
WOOLxvetor_mc_resultados_media_pop=vector()
WOOLxvetor_mc_resultados_media_est=vector()
WOOLxvetor_mc_resultados_sd_vies=vector()
WOOLxvetor_mc_resultados_mcsd=vector()

WOOLyvetor_mc_resultados_vies=vector()
WOOLyvetor_mc_resultados_sd=vector()
WOOLyvetor_mc_resultados_media_pop=vector()
WOOLyvetor_mc_resultados_media_est=vector()
WOOLyvetor_mc_resultados_sd_vies=vector()
WOOLyvetor_mc_resultados_mcsd=vector()

c=0
for (k in seq(1,ncol(matriz_resultados)/3,2)){
  c=c+1
  WOOLvetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
  WOOLvetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
  WOOLvetor_mc_resultados_mcsd[c]=var(matriz_resultados[,k+1])*sqrt((kurtosis(matriz_resultados[,k+1])-1)/N)
  WOOLvetor_mc_resultados_sd[c]=sd(matriz_resultados[,k+1])
  WOOLvetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
  WOOLvetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
}

c=0
for (k in seq(ncol(matriz_resultados)/3+1,ncol(matriz_resultados)*2/3,2)){
  c=c+1
  WOOLxvetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
  WOOLxvetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
  WOOLxvetor_mc_resultados_mcsd[c]=var(matriz_resultados[,k+1])*sqrt((kurtosis(matriz_resultados[,k+1])-1)/N)
  WOOLxvetor_mc_resultados_sd[c]=sd(matriz_resultados[,k+1])
  WOOLxvetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
  WOOLxvetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
}


c=0
for (k in seq(ncol(matriz_resultados)*2/3+1,ncol(matriz_resultados),2)){
  c=c+1
  WOOLyvetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
  WOOLyvetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
  WOOLyvetor_mc_resultados_mcsd[c]=var(matriz_resultados[,k+1])*sqrt((kurtosis(matriz_resultados[,k+1])-1)/N)
  WOOLyvetor_mc_resultados_sd[c]=sd(matriz_resultados[,k+1])
  WOOLyvetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
  WOOLyvetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
}


par(mfrow=c(1,2))
hist(matriz_resultados[,8],main="Distribuição estimador Wooldridge t=4 e sd 2 no modelo DiD")
hist(matriz_resultados[,8],main="Distribuição estimador CIC t=4 e sd 2 no modelo DiD")
