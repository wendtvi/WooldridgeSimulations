vector_treat_effect=vector()
q_parameter=vector()
for(k in 1:1000){
TT=as.numeric(runif(500,0,1)>0.5)
GG=as.numeric(runif(500,0,1)>0.5)
X_star=1.2+0.3*TT+0.7*GG+3.2*GG*TT+rnorm(500,0,1)
X=mean(X_star)
Y=exp(X_star)

data=cbind(TT,GG,X_star,Y,as.numeric(TT+GG>1),as.numeric(TT+GG==0))
Y_00=data[data[,6]==1,4]
Y_11=data[data[,5]==1,4]
Y_10=data[data[,2]==1&data[,5]==0,4]
Y_01=data[data[,1]==1&data[,5]==0,4]

X_00=data[data[,6]==1,3]
X_11=data[data[,5]==1,3]
X_10=data[data[,2]==1&data[,5]==0,3]
X_01=data[data[,1]==1&data[,5]==0,3]

treat_effect=log(qlnorm(plnorm(mean(Y_10),meanlog = mean(X_10), sdlog = sd(X_10)),meanlog = mean(X_11), sdlog = sd(X_11)))-
                 (log(qlnorm(plnorm(mean(Y_10),meanlog = mean(X_10), sdlog = sd(X_10)),meanlog = mean(X_10), sdlog = sd(X_10)))+
                    log(qlnorm(plnorm(mean(Y_10),meanlog = mean(X_10), sdlog = sd(X_10)),meanlog = mean(X_01), sdlog = sd(X_01)))-
                    log(qlnorm(plnorm(mean(Y_10),meanlog = mean(X_10), sdlog = sd(X_10)),meanlog = mean(X_00), sdlog = sd(X_00))))

vector_treat_effect[k]=treat_effect
q_parameter[k]=log(qlnorm(plnorm(exp(treat_effect+1.2+0.3+0.7),meanlog = mean(X_11-1.2-0.3-0.7), sdlog = sd(X_11)), meanlog = mean(X_11-1.2-0.3-0.7),sdlog = sd(X_11)))-1.2-0.3-0.7
}

par(mfrow=c(1,2))
hist(vector_treat_effect,main = "Delta_1 generated as in (1)",xlab = "Delta_1")
hist(q_parameter,main = "Delta_1 generated as the inverse of (2)",xlab = "Delta_1")


