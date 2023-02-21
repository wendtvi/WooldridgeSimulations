mc_function=function(N){
  n1=500
  TT=6
  q=4
  matriz_resultados=matrix(NA, ncol = 2*((TT-q+1)*2),nrow = N)
  
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
    

    
    #Modelo de diferenças em diferença para cada cohort (01 primero é grupo e segundo é tempo: neste caso é grupo de controle no período pós trat)
    f5=rep(1,n1)
    f6=rep(1,n1)
    X_estrela_0=-(Z_cov_mean-1)*D/4-0.02*100+rnorm(n1,0,1) #para t=1,...,TT caso não houvesse tratamento
    X_estrela_1_ntratados=X_estrela_0 #se houvesse tratamento, mas var observada antes do tratamento
    X_estrela_1_t4=-(Z_cov_mean-1)*D/4-0.02*100+rnorm(1,0,1) #para t=4
    X_estrela_1_t5=-(Z_cov_mean-1)*D/4-0.02*100+0.2*f5*D+rnorm(1,0,1) #para t=5
    X_estrela_1_t6=-(Z_cov_mean-1)*D/4-0.02*100+0.2*f5*D+0.3*f6*D+rnorm(1,0,1) #para t=6
    

    
    matriz_estado_naotratamento=matrix(NA,nrow = n1,ncol = TT+1)
    matriz_estado_naotratamento[,TT+1]=t(D)
    for (k in 1:TT){
      if (k<q) matriz_estado_naotratamento[,k]=-(Z_cov_mean-1)*D/4-0.02*100+rnorm(n1,0,1) 
      if (k>=q) matriz_estado_naotratamento[,k]=-(Z_cov_mean-1)*D/4-0.02*100+rnorm(n1,0,1) 
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
    
    #Tranformação não linear W() exp()
    matriz_X=matriz_X_estrela
    for (k in 1:ncol(matriz_X)-1){
      for (i in 1:nrow(matriz_X)){
        matriz_X[i,k]=exp(matriz_X_estrela[i,k])
      }
    }
    
    sample_00=vector()
    sample_01=vector()
    sample_10=vector()
    sample11_4=vector()
    sample11_5=vector()
    sample11_6=vector()
    for (k in 1:n1){
      sample_00[k]=rbinom(1,1,exp(mean(-(Z_cov_mean-1)*D/4-0.02*100)))
      sample_01[k]=rbinom(1,1,exp(mean(-(Z_cov_mean-1)*D/4-0.02*100)))
      sample_10[k]=rbinom(1,1,exp(mean(-(Z_cov_mean-1)*D/4-0.02*100)))
      sample11_4[k]=rbinom(1,1,exp(mean(-(Z_cov_mean-1)*D/4-0.02*100+0.00)))
      sample11_5[k]=rbinom(1,1,exp(mean(-(Z_cov_mean-1)*D/4-0.02*100+0.20)))
      sample11_6[k]=rbinom(1,1,exp(mean(-(Z_cov_mean-1)*D/4-0.02*100+0.30)))
    }    
    
    
    #Tranformação não linear S() F()
    matriz_Y=matriz_X
    for (k in 1:ncol(matriz_Y)-1){
      for (i in 1:nrow(matriz_Y)){
        if(matriz_Y[i,TT+1]==0 && k<q)matriz_Y[i,k]=pbinom(sample_00[matriz_Y[,TT+1]==0],size=length(sample_00[matriz_Y[,TT+1]==0]),prob = mean(sample_00[matriz_Y[,TT+1]==0]))[i]
        if(matriz_Y[i,TT+1]==1 && k<q)matriz_Y[i,k]=pbinom(sample_10[matriz_Y[,TT+1]==1],size=length(sample_10[matriz_Y[,TT+1]==1]),prob = mean(sample_10[matriz_Y[,TT+1]==1]))[i]
        if(matriz_Y[i,TT+1]==0 && k>=q)matriz_Y[i,k]=pbinom(sample_01[matriz_Y[,TT+1]==0],size=length(sample_01[matriz_Y[,TT+1]==0]),prob = mean(sample_01[matriz_Y[,TT+1]==0]))[i]
        if(matriz_Y[i,TT+1]==1 && k==q)matriz_Y[i,k]=pbinom(sample11_4[matriz_Y[,TT+1]==1],size=length(sample11_4[matriz_Y[,TT+1]==1]),prob = mean(sample11_4[matriz_Y[,TT+1]==1]))[i]
        if(matriz_Y[i,TT+1]==1 && k==q+1)matriz_Y[i,k]=pbinom(sample11_5[matriz_Y[,TT+1]==1],size=length(sample11_5[matriz_Y[,TT+1]==1]),prob = mean(sample11_5[matriz_Y[,TT+1]==1]))[i]
        if(matriz_Y[i,TT+1]==1 && k==q+2)matriz_Y[i,k]=pbinom(sample11_6[matriz_Y[,TT+1]==1],size=length(sample11_6[matriz_Y[,TT+1]==1]),prob = mean(sample11_6[matriz_Y[,TT+1]==1]))[i]
      }
    }
    
    #hist(matriz_Y[matriz_Y[,TT+1]==0,q-1],breaks = 200,main="Histograma grupo de controle em t=3 (período pré tratamento)")
  
    
    ##########################################################
    #####################RESULTADOS###########################
    ##########################################################
    tau_4=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q])-mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q]) #resultados populacionais
    tau_5=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+1])-mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q+1])
    tau_6=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+2])-mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q+2])
    
    tau_4_hat=log(mean(sample11_4[matriz_Y[,TT+1]==1])/(mean(sample_01[matriz_Y[,TT+1]==0])/mean(sample_00[matriz_Y[,TT+1]==0])*mean(sample_10[matriz_Y[,TT+1]==1])))
    tau_5_hat=log(mean(sample11_5[matriz_Y[,TT+1]==1])/(mean(sample_01[matriz_Y[,TT+1]==0])/mean(sample_00[matriz_Y[,TT+1]==0])*mean(sample_10[matriz_Y[,TT+1]==1])))
    tau_6_hat=log(mean(sample11_6[matriz_Y[,TT+1]==1])/(mean(sample_01[matriz_Y[,TT+1]==0])/mean(sample_00[matriz_Y[,TT+1]==0])*mean(sample_10[matriz_Y[,TT+1]==1])))

    
    gamma_4=exp(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q]))-exp(mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q]))
    gamma_5=exp(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+1]))-exp(mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q+1]))
    gamma_6=exp(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q+2]))-exp(mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,q+2]))
    
    gamma4_hat=exp(log(mean(sample_01[matriz_Y[,TT+1]==0])/mean(sample_00[matriz_Y[,TT+1]==0])*mean(sample_10[matriz_Y[,TT+1]==1]))+tau_4_hat)-(mean(sample_01[matriz_Y[,TT+1]==0])/mean(sample_00[matriz_Y[,TT+1]==0])*mean(sample_10[matriz_Y[,TT+1]==1]))
    gamma5_hat=exp(log(mean(sample_01[matriz_Y[,TT+1]==0])/mean(sample_00[matriz_Y[,TT+1]==0])*mean(sample_10[matriz_Y[,TT+1]==1]))+tau_5_hat)-(mean(sample_01[matriz_Y[,TT+1]==0])/mean(sample_00[matriz_Y[,TT+1]==0])*mean(sample_10[matriz_Y[,TT+1]==1]))
    gamma6_hat=exp(log(mean(sample_01[matriz_Y[,TT+1]==0])/mean(sample_00[matriz_Y[,TT+1]==0])*mean(sample_10[matriz_Y[,TT+1]==1]))+tau_6_hat)-(mean(sample_01[matriz_Y[,TT+1]==0])/mean(sample_00[matriz_Y[,TT+1]==0])*mean(sample_10[matriz_Y[,TT+1]==1]))
    
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
    
  }
  return(matriz_resultados)
}
