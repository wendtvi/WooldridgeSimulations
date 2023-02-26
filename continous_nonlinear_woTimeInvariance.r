mc_function=function(N){
  n1=500
  TT=6
  q=4
  matriz_resultados=matrix(NA, ncol = 3*((TT-q+1)*2),nrow = N)
  
  for(p in 1:N){
    ##########################################################
    #####################CENÁRIO SIMULAÇÃO####################
    ##########################################################
    
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
    U0=rnorm(n1,0,1.81)
    U1=rnorm(n1,0,1.81)
    epsolon=rnorm(n1,0,1)
    
    #Modelo de diferenças em diferença para cada cohort (01 primero é grupo e segundo é tempo: neste caso é grupo de controle no período pós trat)
    f5=rep(1,n1)
    f6=rep(1,n1)
    X_estrela_0=(Z_cov_mean-1)/2-2*D+(Z_cov_mean-1)*D/4+(epsolon+rnorm(n1,0,1.81)) #para t=1,...,TT caso não houvesse tratamento
    X_estrela_1_ntratados=X_estrela_0 #se houvesse tratamento, mas var observada antes do tratamento
    X_estrela_1_t4=(Z_cov_mean-1)*3/4-2*D+(epsolon+rnorm(n1,0,1.81)) #para t=4
    X_estrela_1_t5=(Z_cov_mean-1)*3/4-2*D+0.2*f5*D+(epsolon+rnorm(n1,0,1.81)) #para t=5
    X_estrela_1_t6=(Z_cov_mean-1)*3/4-2*D+0.2*f5*D+0.3*f6*D+(epsolon+rnorm(n1,0,1.81)) #para t=6
    
    matriz_estado_naotratamento=matrix(NA,nrow = n1,ncol = TT+1)
    matriz_estado_naotratamento[,TT+1]=t(D)
    for (k in 1:TT){
      if (k<q) X_estrela_0=(Z_cov_mean-1)/2-2*D+(Z_cov_mean-1)*D/4+(epsolon+rnorm(n1,0,1.81))
      if (k>=q) X_estrela_0=(Z_cov_mean-1)/2-2*D+(Z_cov_mean-1)*D/4+(epsolon+rnorm(n1,0,1.81))
      matriz_estado_naotratamento[,k]=t(X_estrela_0)
    }
    
    #Gerando matriz de dados latentes.
    matriz_X_estrela=matrix(NA, ncol = TT+1,nrow = n1)
    matriz_X_estrela[,TT+1]=t(D)
    for (k in 1:ncol(matriz_X_estrela)-1){
      matriz_X_estrela[,k]=matriz_estado_naotratamento[,k]     
      for (i in 1:nrow(matriz_X_estrela)){
        if (k==4 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=X_estrela_1_t4[i]
        if (k==5 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=X_estrela_1_t5[i]
        if (k==6 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=X_estrela_1_t6[i]
      }
    }
    
    #Transformação W() apenas replica X estrela
    #Tranformação não linear S() exp()
    matriz_Y=matriz_X_estrela
    for (k in 1:ncol(matriz_Y)-1){
      for (i in 1:nrow(matriz_Y)){
        matriz_Y[i,k]=exp(matriz_X_estrela[i,k])
      }
    }
    #hist(matriz_Y[matriz_Y[,TT+1]==0,q-1],breaks = 200,main="Histograma grupo de controle em t=3 (período pré tratamento)")
    
    ##########################################################
    #####################ESTIMADORES##########################
    ##########################################################
    Y10= matriz_Y[matriz_Y[,TT+1]==1,(q-1)] #variável resposta observada para grupo dos tratados no período pré tratamento t=3
    #Suponho que sei que variável latente segue distribuição logistica com parâmetros 0,1
    F_Y10=plnorm((Y10),meanlog =mean(log(matriz_Y[matriz_Y[,TT+1]==0,(q-1)])),
                 sdlog =sqrt(var(log(matriz_Y[matriz_Y[,TT+1]==0,(q-1)]))))
    F_Y10[F_Y10==1]=0.9999999999
    #hist(F_Y10)
    F_inver_F_Y10_t4=qlnorm(F_Y10,meanlog =mean(log(matriz_Y[matriz_Y[,TT+1]==0,(q)])),
                         sdlog = sqrt(var(log(matriz_Y[matriz_Y[,TT+1]==0,(q)]))))
    F_inver_F_Y10_t5=qlnorm(F_Y10,meanlog =mean(log(matriz_Y[matriz_Y[,TT+1]==0,(q+1)])),
                         sdlog = sqrt(var(log(matriz_Y[matriz_Y[,TT+1]==0,(q+1)]))))
    F_inver_F_Y10_t6=qlnorm(F_Y10,meanlog =mean(log(matriz_Y[matriz_Y[,TT+1]==0,(q+2)])),
                         sdlog = sqrt(var(log(matriz_Y[matriz_Y[,TT+1]==0,(q+2)]))))
    
    
    ##########################################################
    #####################RESULTADOS###########################
    ##########################################################
    tau_4=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q])-mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q]) #resultados populacionais
    tau_5=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+1])-mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q+1])
    tau_6=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+2])-mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q+2])
    
    tau_4_hat=log(exp(mean(log(matriz_Y[matriz_Y[,TT+1]==1,q])))/exp(mean(log(F_inver_F_Y10_t4))))
    tau_5_hat=log(exp(mean(log(matriz_Y[matriz_Y[,TT+1]==1,q+1])))/exp(mean(log(F_inver_F_Y10_t5))))
    tau_6_hat=log(exp(mean(log(matriz_Y[matriz_Y[,TT+1]==1,q+2])))/exp(mean(log(F_inver_F_Y10_t6))))
    
    gamma_4=exp(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q]+var(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q])/2))-exp(mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q-1]+var(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q-1])/2))
    gamma_5=exp(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+1]+var(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+1])/2))-exp(mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q-1]+var(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q-1])/2))
    gamma_6=exp(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+2]+var(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+2])/2))-exp(mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q-1]+var(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q-1])/2))
    
    
    gamma4_hat=exp(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q]+var(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q])/2))-mean(F_inver_F_Y10_t4)
    gamma5_hat=exp(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+1]+var(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+1])/2))-mean(F_inver_F_Y10_t5)
    gamma6_hat=exp(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+2]+var(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+2])/2))-mean(F_inver_F_Y10_t6)
    
    kappa_4=mean(exp(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q]))-mean(exp(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q-1]))
    kappa_5=mean(exp(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+1]))-mean(exp(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q-1]))
    kappa_6=mean(exp(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+2]))-mean(exp(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q-1]))
    
    
    kappa4_hat=mean(exp(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q]))-mean(F_inver_F_Y10_t4)
    kappa5_hat=mean(exp(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+1]))-mean(F_inver_F_Y10_t5)
    kappa6_hat=mean(exp(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+2]))-mean(F_inver_F_Y10_t6)
    
    matriz_resultados[p,1]=tau_4
    matriz_resultados[p,2]=tau_4_hat
    matriz_resultados[p,3]=tau_5
    matriz_resultados[p,4]=tau_5_hat
    matriz_resultados[p,5]=tau_6
    matriz_resultados[p,6]=tau_6_hat
    matriz_resultados[p,7]=gamma_4
    matriz_resultados[p,8]=gamma4_hat
    matriz_resultados[p,9]=gamma_5
    matriz_resultados[p,10]=gamma5_hat
    matriz_resultados[p,11]=gamma_6
    matriz_resultados[p,12]=gamma6_hat
    matriz_resultados[p,13]=kappa_4
    matriz_resultados[p,14]=kappa4_hat
    matriz_resultados[p,15]=kappa_5
    matriz_resultados[p,16]=kappa5_hat
    matriz_resultados[p,17]=kappa_6
    matriz_resultados[p,18]=kappa6_hat
    
  }
  return(matriz_resultados)
}
