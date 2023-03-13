##############################################
###############MCMC###########################
##############################################
#install.packages("devtools")
library(devtools)
#install.packages("moments")
library(moments)

#source_url("https://github.com/wendtvi/dissertation_code/blob/main/continous_nonlinear.r?raw=TRUE")
N=1000
full_return=mc_function(N=1000)
WOOLmatriz_resultados=full_return[[1]]
matriz_resultados=full_return[[2]]

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

simulation_residuals=matrix(0,nrow = N,ncol = 2)

c=0
for (k in seq(1,ncol(WOOLmatriz_resultados)/3,2)){
  c=c+1
  WOOLvetor_mc_resultados_vies[c]=mean(WOOLmatriz_resultados[,k])-mean(WOOLmatriz_resultados[,k+1])
  WOOLvetor_mc_resultados_sd_vies[c]=mean(abs(WOOLmatriz_resultados[,k]-(WOOLmatriz_resultados[,k+1])))
  WOOLvetor_mc_resultados_mcsd[c]=var(WOOLmatriz_resultados[,k+1])*sqrt((kurtosis(WOOLmatriz_resultados[,k+1])-1)/N)
  WOOLvetor_mc_resultados_sd[c]=sd(WOOLmatriz_resultados[,k+1])
  WOOLvetor_mc_resultados_media_pop[c]=mean(WOOLmatriz_resultados[,k])
  WOOLvetor_mc_resultados_media_est[c]=mean(WOOLmatriz_resultados[,k+1])
}

c=0
for (k in seq(ncol(WOOLmatriz_resultados)/3+1,ncol(WOOLmatriz_resultados)*2/3,2)){
  c=c+1
  WOOLxvetor_mc_resultados_vies[c]=mean(WOOLmatriz_resultados[,k])-mean(WOOLmatriz_resultados[,k+1])
  WOOLxvetor_mc_resultados_sd_vies[c]=mean(abs(WOOLmatriz_resultados[,k]-(WOOLmatriz_resultados[,k+1])))
  WOOLxvetor_mc_resultados_mcsd[c]=var(WOOLmatriz_resultados[,k+1])*sqrt((kurtosis(WOOLmatriz_resultados[,k+1])-1)/N)
  WOOLxvetor_mc_resultados_sd[c]=sd(WOOLmatriz_resultados[,k+1])
  WOOLxvetor_mc_resultados_media_pop[c]=mean(WOOLmatriz_resultados[,k])
  WOOLxvetor_mc_resultados_media_est[c]=mean(WOOLmatriz_resultados[,k+1])
  if(simulation_residuals[1]!=0)simulation_residuals[1]=abs(WOOLxvetor_mc_resultados_vies[c]-simulation_residuals[1])
  if(simulation_residuals[1]==0)simulation_residuals[1]=WOOLxvetor_mc_resultados_vies[c]
}


c=0
for (k in seq(ncol(WOOLmatriz_resultados)*2/3+1,ncol(WOOLmatriz_resultados),2)){
  c=c+1
  WOOLyvetor_mc_resultados_vies[c]=mean(WOOLmatriz_resultados[,k])-mean(WOOLmatriz_resultados[,k+1])
  WOOLyvetor_mc_resultados_sd_vies[c]=mean(abs(WOOLmatriz_resultados[,k]-(WOOLmatriz_resultados[,k+1])))
  WOOLyvetor_mc_resultados_mcsd[c]=var(WOOLmatriz_resultados[,k+1])*sqrt((kurtosis(WOOLmatriz_resultados[,k+1])-1)/N)
  WOOLyvetor_mc_resultados_sd[c]=sd(WOOLmatriz_resultados[,k+1])
  WOOLyvetor_mc_resultados_media_pop[c]=mean(WOOLmatriz_resultados[,k])
  WOOLyvetor_mc_resultados_media_est[c]=mean(WOOLmatriz_resultados[,k+1])
}


vetor_mc_resultados_vies=vector()
vetor_mc_resultados_sd=vector()
vetor_mc_resultados_media_pop=vector()
vetor_mc_resultados_media_est=vector()
vetor_mc_resultados_sd_vies=vector()
vetor_mc_resultados_mcsd=vector()

xvetor_mc_resultados_vies=vector()
xvetor_mc_resultados_sd=vector()
xvetor_mc_resultados_media_pop=vector()
xvetor_mc_resultados_media_est=vector()
xvetor_mc_resultados_sd_vies=vector()
xvetor_mc_resultados_mcsd=vector()

yvetor_mc_resultados_vies=vector()
yvetor_mc_resultados_sd=vector()
yvetor_mc_resultados_media_pop=vector()
yvetor_mc_resultados_media_est=vector()
yvetor_mc_resultados_sd_vies=vector()
yvetor_mc_resultados_mcsd=vector()

c=0
for (k in seq(1,ncol(matriz_resultados)/3,2)){
  c=c+1
  vetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
  vetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
  vetor_mc_resultados_mcsd[c]=var(matriz_resultados[,k+1])*sqrt((kurtosis(matriz_resultados[,k+1])-1)/N)
  vetor_mc_resultados_sd[c]=sd(matriz_resultados[,k+1])
  vetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
  vetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
}

c=0
for (k in seq(ncol(matriz_resultados)/3+1,ncol(matriz_resultados)*2/3,2)){
  c=c+1
  xvetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
  xvetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
  xvetor_mc_resultados_mcsd[c]=var(matriz_resultados[,k+1])*sqrt((kurtosis(matriz_resultados[,k+1])-1)/N)
  xvetor_mc_resultados_sd[c]=sd(matriz_resultados[,k+1])
  xvetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
  xvetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
  if(simulation_residuals[2]!=0)simulation_residuals[2]=abs(WOOLxvetor_mc_resultados_vies[c]-simulation_residuals[2])
  if(simulation_residuals[2]==0)simulation_residuals[2]=WOOLxvetor_mc_resultados_vies[c]
}


c=0
for (k in seq(ncol(matriz_resultados)*2/3+1,ncol(matriz_resultados),2)){
  c=c+1
  yvetor_mc_resultados_vies[c]=mean(matriz_resultados[,k])-mean(matriz_resultados[,k+1])
  yvetor_mc_resultados_sd_vies[c]=mean(abs(matriz_resultados[,k]-(matriz_resultados[,k+1])))
  yvetor_mc_resultados_mcsd[c]=var(matriz_resultados[,k+1])*sqrt((kurtosis(matriz_resultados[,k+1])-1)/N)
  yvetor_mc_resultados_sd[c]=sd(matriz_resultados[,k+1])
  yvetor_mc_resultados_media_pop[c]=mean(matriz_resultados[,k])
  yvetor_mc_resultados_media_est[c]=mean(matriz_resultados[,k+1])
}



###################################################################
######################MEAN SQUARED ERROR###########################
###################################################################
MSE_woold=WOOLxvetor_mc_resultados_vies^2+WOOLxvetor_mc_resultados_sd^2
MSE_cic=xvetor_mc_resultados_vies^2+xvetor_mc_resultados_sd^2
MSE_woold
MSE_cic
MSE_woold/MSE_cic
