#install.packages("NTS")
library(NTS)

G_suporte=c(1,0)
T_suporte=c(1,0)
n1=100
n0=100
TT=20
epsolon_tratados=rnorm(n1,0,1)
epsolon_controle=rnorm(n2,1,1)
beta0=rep(n1,0)
beta1=rep(n1,0.3)
gama0=rep(n1,0.4)
gama1=rep(n1,0.1)

#Modelo de diferenças em diferença para cada cohort
X_estrela_11=beta0+beta1*g+gama0*t+gama1*g*t+epsolon_tratados
X_estrela_10=beta0+beta1*g+epsolon_tratados
X_estrela_01=beta0+gama0*t+epsolon_controle
X_estrela_00=beta0+epsolon_controle

#Gerando séries latentes para grupo de controle e tratados
for (i in 1:TT){
  
}

#Tranformação não linear exp()
X_11=exp(X_estrela_11)
X_10=exp(X_estrela_10)
X_01=exp(X_estrela_01)
X_00=exp(X_estrela_00)

nnet()
