mc_function=function(N){
for(p in 1:N){
##########################################################
#####################CENÁRIO SIMULAÇÃO####################
##########################################################
n1=250
n0=250
TT=6
q=4

#Gerando covariável
Z_cov_mean=vector()
Z_cov=matrix(NA,ncol=TT,nrow=n1)
for (k in 1:nrow(Z_cov)){
  Z_cov[k,]=rexp(TT,1)
  Z_cov_mean[k]=mean(Z_cov[k,])
}

#Variável de tratamento
V=rlogis(n1,0,1)
D=vector()
D=rep(-0.5,n1)+(Z_cov_mean-1)+V>0
sum(D)/length(D) #Incidência de tratamento na população

#Termos de erro U
U0=rlogis(n1,0,1)
U1=rlogis(n1,0,1)

#Modelo de diferenças em diferença para cada cohort (01 primero é grupo e segundo é tempo: neste caso é grupo de controle no período pós trat)
f5=rep(1,n1)
f6=rep(1,n1)
X_estrela_0=(Z_cov_mean-1)/2-2*D+(Z_cov_mean-1)*D/4+U0 #para t=1,...,TT
X_estrela_1_ntratados=rep(0.5,n0)+(Z_cov_mean-1)+U1 #nao tratados após tratamento
X_estrela_1_t4=rep(0.5,n1)+(Z_cov_mean-1)-2*D+0.2*f5+0.3*f6+U1 #para t=4
X_estrela_1_t5=rep(0.5,n1)+(Z_cov_mean-1)-2*D+0.2*f5+U1 #para t=5
X_estrela_1_t6=rep(0.5,n1)+(Z_cov_mean-1)-2*D+0.2*f5+0.3*f6+U1 #para t=6
  
#Gerando matriz de dados latentes.
matriz_X_estrela=matrix(NA, ncol = TT+1,nrow = n1+n0)
matriz_X_estrela[,TT+1]=t(c(rep(1,n1),rep(0,n0)))
for (k in 1:ncol(matriz_X_estrela)-1){
  matriz_X_estrela[,k]=t(c(X_estrela_0,X_estrela_0))
  if (k==4)matriz_X_estrela[,k]=t(c(X_estrela_1_t4,X_estrela_1_ntratados))
  if (k==5)matriz_X_estrela[,k]=t(c(X_estrela_1_t5,X_estrela_1_ntratados))
  if (k==6)matriz_X_estrela[,k]=t(c(X_estrela_1_t6,X_estrela_1_ntratados))
}

#Transformação W() apenas replica X estrela
#Tranformação não linear S() exp()
matriz_Y=matriz_X_estrela
for (k in 1:ncol(matriz_Y)-1){
  for (i in 1:nrow(matriz_Y)){
    matriz_Y[i,k]=exp(matriz_X_estrela[i,k])
  }
}
hist(matriz_Y[matriz_Y[,TT+1]==0,q-1],breaks = 200,main="Histograma grupo de controle em t=3 (período pré tratamento)")

##########################################################
#####################ESTIMADORES##########################
##########################################################
Y10= matriz_Y[matriz_Y[,TT+1]==1,(q-1)] #variável resposta observada para grupo dos tratados no período pré tratamento t=3
#Suponho que sei que variável latente segue distribuição logistica com parâmetros 0,1
F_Y10=plogis(Y10,location=mean((Z_cov_mean-1)/2),scale = 1)
hist(F_Y10)
F_inver_F_Y10=qlogis(F_Y10,location = mean(rep(0.5,n1)+(Z_cov_mean-1)),scale = 1 )
plot(F_inver_F_Y10)



##########################################################
#####################RESULTADOS###########################
##########################################################
tau_4=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q])-mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q]) #resultados populacionais
tau_5=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+1])-mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q+1])
tau_6=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+2])-mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,q+2])

tau_4_hat=mean(log(matriz_Y[matriz_Y[,TT+1]==1,q]))-mean(log(F_inver_F_Y10[F_inver_F_Y10!=Inf]))
tau_5_hat=mean(log(matriz_Y[matriz_Y[,TT+1]==1,q+1]))-mean(log(F_inver_F_Y10[F_inver_F_Y10!=Inf]))
tau_6_hat=mean(log(matriz_Y[matriz_Y[,TT+1]==1,q+2]))-mean(log(F_inver_F_Y10[F_inver_F_Y10!=Inf]))

matriz_resultados=matrix(NA, ncol = (TT-q+1)*2,nrow = N)
matriz_resultados[p,1]=tau_4
matriz_resultados[p,2]=tau_4_hat
matriz_resultados[p,3]=tau_5
matriz_resultados[p,4]=tau_5_hat
matriz_resultados[p,5]=tau_6
matriz_resultados[p,6]=tau_6_hat

}
  return(matriz_resultados)
}
