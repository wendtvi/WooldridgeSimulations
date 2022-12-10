#Wooldridge scenario for simulation 1
#set.seed(1234)
TT=6
q=4
n_simul=1000
N=500
X=matrix(NA,nrow = TT,ncol = N)
V=matrix(NA,nrow = TT,ncol = N)
D=matrix(NA,nrow = TT,ncol = N)
prob_D=vector()
prob_Y_0t=vector()
prob_Y_1t=vector()
C=0
U_0t=matrix(NA,nrow = TT,ncol = N)
U_1t=matrix(NA,nrow = TT,ncol = N)
Y_0t=matrix(NA,nrow = TT,ncol = N)
Y_1t=matrix(NA,nrow = TT,ncol = N)
Y_1t_temp=matrix(NA,nrow = TT,ncol = N)
f5_t=rep(0,TT)
f5_t[5]=1
f6_t=rep(0,TT)
f6_t[6]=1
N=500
list_simul_Y_0t=list()
list_simul_Y_1t=list()
for (i in 1:n_simul){
  for (k in 1:N){
    X[,k]=(mean(rexp(TT,1)))
    V[,k]=(rlogis(TT,0,1))
    D[,k]=-0.5+(X[,k]-1)+V[,k]>0
    prob_D[k]=sum(D[,k]==TRUE)/TT
    U_0t[,k]=(rep(C,TT)+sqrt(2)*rlogis(TT,0,1))/2
    Y_0t[,k]=((X[,k]-1)/2-2*D[,k]+(X[,k]-1)*(D[,k]/4)+U_0t[,k])>0
    U_1t[,k]=(rep(C,TT)+(sqrt(2))*rlogis(TT,0,1))/2
    Y_1t[1:(q-1),k]=Y_0t[1:(q-1),k]
    Y_1t_temp[,k]=(0.5+(X[,k]-1)-2*D[,k]+0.2*f5_t+0.3*f6_t+U_1t[,k])>0
    Y_1t[q:TT,k]=Y_1t_temp[q:TT,k]
    prob_Y_1t[k]=sum(Y_1t[,k]==TRUE)/TT
    prob_Y_0t[k]=sum(Y_0t[,k]==TRUE)/TT
    
  
  }
  list_simul_Y_0t[[i]]=Y_0t
  list_simul_Y_1t[[i]]=Y_1t
}

#check if it is equal to Wooldridge simulation
mean(prob_D) #ok
mean(prob_Y_1t)#ok
mean(prob_Y_0t)#ok

#sample ATT
Y_04=vector()
Y_05=vector()
Y_06=vector()
Y_14=vector()
Y_15=vector()
Y_16=vector()
for (i in 1:n_simul){
  Y_04[i]=sum(list_simul_Y_0t[[i]][q,]==TRUE)
  Y_05[i]=sum(list_simul_Y_0t[[i]][q+1,]==TRUE)
  Y_06[i]=sum(list_simul_Y_0t[[i]][q+2,]==TRUE)
  Y_14[i]=sum(list_simul_Y_1t[[i]][q,]==TRUE)
  Y_15[i]=sum(list_simul_Y_1t[[i]][q+1,]==TRUE)
  Y_16[i]=sum(list_simul_Y_1t[[i]][q+2,]==TRUE)
}
real_ATT_4=(sum(Y_14)-sum(Y_04))/(1000*500)
real_ATT_5=(sum(Y_15)-sum(Y_05))/(1000*500)
real_ATT_6=(sum(Y_16)-sum(Y_06))/(1000*500)

#CIC estimator
#if the distribution is known
i=k=0
y_10=matrix(NA,nrow = N,ncol = n_simul)
F_Y_01=matrix(NA,nrow = N,ncol = n_simul)
for(k in 1:n_simul){
  for (i in 1:N){
    y_10[i,k]=list_simul_Y_0t[[k]][q-1,i]
    F_Y_01[i,k]=dlogis(y_10[i,k],0,1)
  }
}

y[,k]=Y_0t[q-1,k]
