mc_function=function(N){
  n1=50
  TT=6
  q=4
  matriz_resultados=matrix(NA, ncol = (2*((TT-q+1)*2)),nrow = N)
  WOOLmatriz_resultados=matrix(NA, ncol = (2*((TT-q+1)*2)),nrow = N)
  
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
    
    
    #ct_0_controle_pre=+rbeta(n1,0.00001,.5)#erro idiossincrático estado zero
    #ct_0_controle_pos=+rbeta(n1,0.00001,1)#erro idiossincrático estado zero
    #ct_inf_controle_pre=+rbeta(n1,0.00001,.5)#erro idiossincrático estado inf
    #ct_inf_controle_pos=+rbeta(n1,.00001,1)#erro idiossincrático estado inf
    #ct_0_trat_pre=+rbeta(n1,0.00001,.5)#erro idiossincrático estado zero
    #ct_0_trat_pos=rbeta(n1,0.00001,1)#erro idiossincrático estado zero
    #ct_inf_trat_pre=+rbeta(n1,0.00001,.5)#erro idiossincrático estado inf
    #ct_inf_trat_pos=+rbeta(n1,0.00001,1)#erro idiossincrático estado inf
    
    ct_0_controle_pre=rep(0,n1)
    ct_0_controle_pos=rep(0,n1)
    ct_inf_controle_pre=rep(0,n1)
    ct_inf_controle_pos=rep(0,n1)
    ct_0_trat_pre=rep(0,n1)
    ct_0_trat_pos=rep(0,n1)
    ct_inf_trat_pre=rep(0,n1)
    ct_inf_trat_pos=rep(0,n1)
    
    #Modelo de diferenças em diferença para cada cohort (01 primero é grupo e segundo é tempo: neste caso é grupo de controle no período pós trat)
    f5=rep(1,n1)*D
    f6=rep(1,n1)*D
    matriz_estado_naotratamento=matrix(NA,nrow = n1,ncol = TT+1)
    matriz_estado_naotratamento[,TT+1]=t(D)
    for (k in 1:TT){
      for(i in 1:nrow(matriz_estado_naotratamento)){
        if (k<q && matriz_estado_naotratamento[i,TT+1]==0) matriz_estado_naotratamento[i,k]=as.numeric(-(Z_cov_mean[i]-1)*5/2-2*D[i]+(Z_cov_mean[i]-1)*D[i]/4+rnorm(1,0,1)+ct_0_controle_pre[i])
        if (k>=q&& matriz_estado_naotratamento[i,TT+1]==0) matriz_estado_naotratamento[i,k]=as.numeric(-(Z_cov_mean[i]-1)*5/2-2*D[i]+(Z_cov_mean[i]-1)*D[i]/4+rnorm(1,0,1)+ct_0_controle_pos[i])
        if (k<q && matriz_estado_naotratamento[i,TT+1]==1) matriz_estado_naotratamento[i,k]=as.numeric(-(Z_cov_mean[i]-1)*5/2-2*D[i]+(Z_cov_mean[i]-1)*D[i]/4+rnorm(1,0,1)+ct_0_trat_pre[i])
        if (k>=q&& matriz_estado_naotratamento[i,TT+1]==1) matriz_estado_naotratamento[i,k]=as.numeric(-(Z_cov_mean[i]-1)*5/2-2*D[i]+(Z_cov_mean[i]-1)*D[i]/4+rnorm(1,0,1)+ct_0_trat_pos[i])
        
      }
    }
    
    #Gerando matriz de dados latentes.
    matriz_X_estrela=matrix(NA, ncol = TT+1,nrow = n1)
    matriz_X_estrela[,TT+1]=t(D)
    for (k in 1:ncol(matriz_X_estrela)-1){
      for (i in 1:nrow(matriz_X_estrela)){
        if (k<q && matriz_X_estrela[i,TT+1]==0) matriz_X_estrela[i,k]=as.numeric(-(Z_cov_mean[i]-1)*5/2-2*D[i]+(Z_cov_mean[i]-1)*D[i]/4+rnorm(1,0,1)+ct_0_controle_pre[i])
        if (k>=q&& matriz_X_estrela[i,TT+1]==0) matriz_X_estrela[i,k]=as.numeric(-(Z_cov_mean[i]-1)*5/2-2*D[i]+(Z_cov_mean[i]-1)*D[i]/4+rnorm(1,0,1)+ct_0_controle_pos[i])
        if (k<q && matriz_X_estrela[i,TT+1]==1) matriz_X_estrela[i,k]=as.numeric(-(Z_cov_mean[i]-1)*5/2-2*D[i]+(Z_cov_mean[i]-1)*D[i]/4+rnorm(1,0,1)+ct_0_trat_pre[i])
        if (k>=q&& matriz_X_estrela[i,TT+1]==1) matriz_X_estrela[i,k]=as.numeric(-(Z_cov_mean[i]-1)*5/2-2*D[i]+(Z_cov_mean[i]-1)*D[i]/4+rnorm(1,0,1)+ct_0_trat_pos[i])
        if (k==4 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=as.numeric(0.5+(Z_cov_mean[i]-1)-2*D[i]+rnorm(1,0,1)+ct_inf_trat_pos[i])
        if (k==5 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=as.numeric(0.5+(Z_cov_mean[i]-1)-2*D[i]+0.2*f5[i]+rnorm(1,0,1)+ct_inf_trat_pos[i])
        if (k==6 && matriz_X_estrela[i,TT+1]==1)matriz_X_estrela[i,k]=as.numeric(0.5+(Z_cov_mean[i]-1)-2*D[i]+0.2*f5[i]+0.3*f6[i]+rnorm(1,0,1)+ct_inf_trat_pos[i])
      }
    }
    
    ####################################################
    ##GERANDO VAR DO ESTADO NAO TRATADO#################
    ####################################################
    #GERANDO VARIAVEL X^*
    X_estrela_vetor_0_0=vector()
    X_estrela_vetor_1_0=vector()
    for (k in 1:TT){
      X_estrela_vetor_0_0[k]=exp(mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==0,k]))
      X_estrela_vetor_1_0[k]=exp(mean(matriz_estado_naotratamento[matriz_estado_naotratamento[,TT+1]==1,k]))
    }
    
    #GERANDO VARIAVEL X
    X_vetor_0_0=vector()
    X_vetor_1_0=vector()
    for (k in 1:TT){
      X_vetor_0_0[k]=log(X_estrela_vetor_0_0[k])
      X_vetor_1_0[k]=log(X_estrela_vetor_1_0[k])
    }
    
    
    
    
    ####################################################
    ##GERANDO VAR DO ESTADO TRATADO#####################
    ####################################################
    #GERANDO VARIAVEL X^*
    X_estrela_vetor_0_inf=vector()
    X_estrela_vetor_1_inf=vector()
    for (k in 1:TT){
      X_estrela_vetor_0_inf[k]=exp(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==0,k]))
      X_estrela_vetor_1_inf[k]=exp(mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,k]))
    }
    
    
    #GERANDO VARIAVEL X
    X_vetor_0_inf=vector()
    X_vetor_1_inf=vector()
    for (k in 1:TT){
      X_vetor_0_inf[k]=log(X_estrela_vetor_0_inf[k])
      X_vetor_1_inf[k]=log(X_estrela_vetor_1_inf[k])
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
    
    G_Ybar_11_t4=log(X_estrela_vetor_1_inf[q])
    G_Ybar_11_t5=log(X_estrela_vetor_1_inf[q+1])
    G_Ybar_11_t6=log(X_estrela_vetor_1_inf[q+2])
    
    G_Ybar_10_t3=log(X_estrela_vetor_1_inf[q-1])
    G_Ybar_00_t3=log(X_estrela_vetor_0_inf[q-1])
    G_Ybar_01_t4=log(X_estrela_vetor_0_inf[q])
    G_Ybar_01_t5=log(X_estrela_vetor_0_inf[q+1])
    G_Ybar_01_t6=log(X_estrela_vetor_0_inf[q+2])
    
    
    ##########################################################
    #####################RESULTADOS###########################
    ##########################################################
    tau_4=X_estrela_vetor_1_inf[q]-X_estrela_vetor_1_0[q]
    tau_5=X_estrela_vetor_1_inf[q+1]-X_estrela_vetor_1_0[q+1]
    tau_6=X_estrela_vetor_1_inf[q+2]-X_estrela_vetor_1_0[q+2]
    
    tau_4_hat=X_estrela_vetor_1_inf[q]-exp(G_Ybar_01_t4-G_Ybar_00_t3+G_Ybar_10_t3)
    tau_5_hat=X_estrela_vetor_1_inf[q+1]-exp((G_Ybar_01_t5-G_Ybar_00_t3+G_Ybar_10_t3))
    tau_6_hat=X_estrela_vetor_1_inf[q+2]-exp((G_Ybar_01_t6-G_Ybar_00_t3+G_Ybar_10_t3))
    
    
    gamma_4=X_vetor_1_inf[q]-X_vetor_1_0[q]
    gamma_5=X_vetor_1_inf[q+1]-X_vetor_1_0[q+1]
    gamma_6=X_vetor_1_inf[q+2]-X_vetor_1_0[q+2]
    
    gamma4_hat=X_vetor_1_inf[q]-((G_Ybar_01_t4-G_Ybar_00_t3+G_Ybar_10_t3))
    gamma5_hat=X_vetor_1_inf[q+1]-((G_Ybar_01_t5-G_Ybar_00_t3+G_Ybar_10_t3))
    gamma6_hat=X_vetor_1_inf[q+2]-((G_Ybar_01_t6-G_Ybar_00_t3+G_Ybar_10_t3))
    
    
    
    WOOLmatriz_resultados[p,1]=tau_4
    WOOLmatriz_resultados[p,2]=tau_4_hat
    WOOLmatriz_resultados[p,3]=tau_5
    WOOLmatriz_resultados[p,4]=tau_5_hat
    WOOLmatriz_resultados[p,5]=tau_6
    WOOLmatriz_resultados[p,6]=tau_6_hat
    WOOLmatriz_resultados[p,7]=gamma_4
    WOOLmatriz_resultados[p,8]=gamma4_hat
    WOOLmatriz_resultados[p,9]=gamma_5
    WOOLmatriz_resultados[p,10]=gamma5_hat
    WOOLmatriz_resultados[p,11]=gamma_6
    WOOLmatriz_resultados[p,12]=gamma6_hat
    
    
    
    ##########################################################
    #####################ESTIMADORES##########################
    ##########################################################
    Y10=mean(matriz_X_estrela[matriz_X_estrela[,TT+1]==1,q-1]) #variável resposta observada para grupo dos tratados no período pré tratamento t=3
    #Suponho que sei que variável latente segue distribuição logistica com parâmetros 0,1
    F_Y10=pnorm((Y10), mean =5 ,sd = 1)
    
    F_inver_F_Y10_t4=qnorm(F_Y10, mean =5 ,sd = 1)
    F_inver_F_Y10_t5=qnorm(F_Y10,mean =5 ,sd = 1)
    F_inver_F_Y10_t6=qnorm(F_Y10,mean =5 ,sd = 1)
    
    
    ##########################################################
    #####################RESULTADOS###########################
    ##########################################################
    tau_4=X_estrela_vetor_1_inf[q]-X_estrela_vetor_1_0[q]
    tau_5=X_estrela_vetor_1_inf[q+1]-X_estrela_vetor_1_0[q+1]
    tau_6=X_estrela_vetor_1_inf[q+2]-X_estrela_vetor_1_0[q+2]
    
    tau_4_hat=X_estrela_vetor_1_inf[q]-exp((F_inver_F_Y10_t4))
    tau_5_hat=X_estrela_vetor_1_inf[q+1]-exp((F_inver_F_Y10_t5))
    tau_6_hat=X_estrela_vetor_1_inf[q+2]-exp((F_inver_F_Y10_t6))
    
    gamma_4=X_vetor_1_inf[q]-X_vetor_1_0[q]
    gamma_5=X_vetor_1_inf[q+1]-X_vetor_1_0[q+1]
    gamma_6=X_vetor_1_inf[q+2]-X_vetor_1_0[q+2]
    
    gamma4_hat=X_vetor_1_inf[q]-(F_inver_F_Y10_t4)
    gamma5_hat=X_vetor_1_inf[q+1]-(F_inver_F_Y10_t5)
    gamma6_hat=X_vetor_1_inf[q+2]-(F_inver_F_Y10_t6)
    
    
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
  return(list(WOOLmatriz_resultados,matriz_resultados))
}
