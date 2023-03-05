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
    pD=sum(D)/length(D) #Incidência de tratamento na população
    
    
    ct_0=rlogis(n1,0,1) #erro idiossincrático estado zero
    ct_inf=rlogis(n1,0,1)#erro idiossincrático estado inf
    
    #Modelo de diferenças em diferença para cada cohort (01 primero é grupo e segundo é tempo: neste caso é grupo de controle no período pós trat)
    f5=rep(1,n1)*D
    f6=rep(1,n1)*D
    
    matriz_estado_naotratamento=matrix(NA,nrow = n1,ncol = TT+1)
    matriz_estado_naotratamento[,TT+1]=t(D)
    for (k in 1:TT){
      if (k<q) matriz_estado_naotratamento[,k]=as.numeric(-(Z_cov_mean-1)/2-2*D+(Z_cov_mean-1)*D/4+rnorm(n1,1,3)+ct_0>0)
      if (k>=q) matriz_estado_naotratamento[,k]=as.numeric(-(Z_cov_mean-1)/2-2*D+(Z_cov_mean-1)*D/4+rnorm(n1,1,3)+ct_0>0)
    }
    
    #Gerando matriz de dados latentes.
    matriz_X_estrela=matrix(NA, ncol = TT+1,nrow = n1)
    matriz_X_estrela[,TT+1]=t(D)
    for (k in 1:ncol(matriz_X_estrela)-1){
      matriz_X_estrela[,k]=matriz_estado_naotratamento[,k]
      if (k<q) matriz_X_estrela[,k]=as.numeric(-(Z_cov_mean-1)/2-2*D+(Z_cov_mean-1)*D/4+rnorm(n1,1,3)+ct_inf>0)
      for (i in 1:nrow(matriz_X_estrela)){
        if (k==4 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=as.numeric(0.5+(Z_cov_mean[i]-1)-2*D[i]+rnorm(1,0,1)+ct_inf[i]>0)
        if (k==5 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=as.numeric(0.5+(Z_cov_mean[i]-1)-2*D[i]+0.2*f5[i]+rnorm(1,0,1)+ct_inf[i]>0)
        if (k==6 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=as.numeric(0.5+(Z_cov_mean[i]-1)-2*D[i]+0.2*f5[i]+0.3*f6[i]+rnorm(1,0,1)+ct_inf[i]>0)
      }
    }
    
    
    ####################################################
    ##GERANDO VAR DO ESTADO NAO TRATADO#################
    ####################################################
    #GERANDO VARIAVEL X^*
    X_estrela_vetor_0_0=vector()
    X_estrela_vetor_1_0=vector()
    for (k in 1:TT){
      X_estrela_vetor_0_0[k]=mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==0,k])/(1-mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==0,k]))
      X_estrela_vetor_1_0[k]=mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,k])/(1-mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,k]))
    }
    
    #GERANDO VARIAVEL X
    X_vetor_0_0=vector()
    X_vetor_1_0=vector()
    for (k in 1:TT){
      X_vetor_0_0[k]= X_estrela_vetor_0_0[k]/(1+ X_estrela_vetor_0_0[k])
      X_vetor_1_0[k]=X_estrela_vetor_1_0[k]/(1+X_estrela_vetor_1_0[k])
    }
    
    
    #GERANDO VARIAVEL Y
    Y_vetor_0_0=vector()
    Y_vetor_1_0=vector()
    for (k in 1:TT){
      Y_vetor_0_0[k]=sum(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==0,k])
      Y_vetor_1_0[k]=sum(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,k])
    }
    
    
    
    ####################################################
    ##GERANDO VAR DO ESTADO TRATADO#####################
    ####################################################
    #GERANDO VARIAVEL X^*
    X_estrela_vetor_0_inf=vector()
    X_estrela_vetor_1_inf=vector()
    for (k in 1:TT){
      X_estrela_vetor_0_inf[k]=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,k])/(1-mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,k]))
      X_estrela_vetor_1_inf[k]=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,k])/(1-mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,k]))
    }
    
    
    #GERANDO VARIAVEL X
    X_vetor_0_inf=vector()
    X_vetor_1_inf=vector()
    for (k in 1:TT){
      X_vetor_0_inf[k]=X_estrela_vetor_0_inf[k]/(1+X_estrela_vetor_0_inf[k])
      X_vetor_1_inf[k]=X_estrela_vetor_1_inf[k]/(1+X_estrela_vetor_1_inf[k])
    }
    

    #GERANDO VARIAVEL Y
    Y_vetor_0_inf=vector()
    Y_vetor_1_inf=vector()
    for (k in 1:TT){
      Y_vetor_0_inf[k]=sum(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,k])
      Y_vetor_1_inf[k]=sum(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,k])
    }
    
    ##########################################################
    #####################ESTIMADORES##########################
    ##########################################################
    Ybar_11_t4=X_vetor_1_inf[q]
    Ybar_11_t5=X_vetor_1_inf[q+1]
    Ybar_11_t6=X_vetor_1_inf[q+2]
    
    Ybar_10_t3=X_vetor_1_inf[q-1]
    Ybar_00_t3=X_vetor_0_inf[q-1]
    Ybar_01_t4=X_vetor_0_inf[q]
    Ybar_01_t5=X_vetor_0_inf[q+1]
    Ybar_01_t6=X_vetor_0_inf[q+2]
    
    G_Ybar_11_t4=X_vetor_1_inf[q]/(1+X_vetor_1_inf[q])
    G_Ybar_11_t5=X_vetor_1_inf[q+1]/(1+X_vetor_1_inf[q+1])
    G_Ybar_11_t6=X_vetor_1_inf[q+2]/(1+X_vetor_1_inf[q+2])
    
    G_Ybar_10_t3=X_vetor_1_inf[q-1]/(1+X_vetor_1_inf[q-1])
    G_Ybar_00_t3=X_vetor_0_inf[q-1]/(1+X_vetor_0_inf[q-1])
    G_Ybar_01_t4=X_vetor_0_inf[q]/(1+X_vetor_0_inf[q])
    G_Ybar_01_t5=X_vetor_0_inf[q+1]/(1+X_vetor_0_inf[q+1])
    G_Ybar_01_t6=X_vetor_0_inf[q+2]/(1+X_vetor_0_inf[q+2])
    
    
    ##########################################################
    #####################RESULTADOS###########################
    ##########################################################
    tau_4=X_estrela_vetor_1_inf[q]-X_estrela_vetor_1_0[q]
    tau_5=X_estrela_vetor_1_inf[q+1]-X_estrela_vetor_1_0[q+1]
    tau_6=X_estrela_vetor_1_inf[q+2]-X_estrela_vetor_1_0[q+2]
    
    tau_4_hat=X_estrela_vetor_1_inf[q]-((G_Ybar_01_t4-G_Ybar_00_t3+G_Ybar_10_t3))
    tau_5_hat=X_estrela_vetor_1_inf[q+1]-((G_Ybar_01_t5-G_Ybar_00_t3+G_Ybar_10_t3))
    tau_6_hat=X_estrela_vetor_1_inf[q+2]-((G_Ybar_01_t6-G_Ybar_00_t3+G_Ybar_10_t3))
    
    gamma_4=X_vetor_1_inf[q]-X_vetor_1_0[q]
    gamma_5=X_vetor_1_inf[q+1]-X_vetor_1_0[q+1]
    gamma_6=X_vetor_1_inf[q+2]-X_vetor_1_0[q+2]
    
    gamma4_hat=X_vetor_1_inf[q]-((G_Ybar_01_t4-G_Ybar_00_t3+G_Ybar_10_t3)/(1-(G_Ybar_01_t4-G_Ybar_00_t3+G_Ybar_10_t3)))
    gamma5_hat=X_vetor_1_inf[q+1]-((G_Ybar_01_t5-G_Ybar_00_t3+G_Ybar_10_t3)/(1-(G_Ybar_01_t5-G_Ybar_00_t3+G_Ybar_10_t3)))
    gamma6_hat=X_vetor_1_inf[q+2]-((G_Ybar_01_t6-G_Ybar_00_t3+G_Ybar_10_t3)/(1-(G_Ybar_01_t6-G_Ybar_00_t3+G_Ybar_10_t3)))
    
    
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
