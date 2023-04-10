vet_mean=vector()  
mu_F_11_0=vector()

for (k in 1:1000){
  F_01=1/2*rlnorm(1000,3,1)+1/2*rlnorm(1000,3,1)
  F_00=1/2*rlnorm(1000,2,1)+1/2*rlnorm(1000,3,1)
  F_10=1/2*rlnorm(1000,2,1)+1/2*rlnorm(1000,4,1)
  F_11_0=1/2*rlnorm(1000,3,1)+1/2*rlnorm(1000,4,1)
  mu_F01=log(median(F_01))
  mu_F00=log(median(F_00))
  mu_F10=log(median(F_10))
  mu_F_11_0[k]=log(median(F_11_0))
  
 # mu_F01-mu_F00+mu_F10+cov(log(F_01),log(F_00))+cov(log(F_10),log(F_00))-cov(log(F_10),log(F_01))
  
  #mean(log(F_10)*log(F_00))-mean(log(F_10))*mean(log(F_00))
  
  model_F_11_0=1/2*rlnorm(1000,mu_F01-mu_F00-cov(log(F_01),log(F_00)),sqrt(var(log(F_01))+var(log(F_00))-2*cov(log(F_01),log(F_00))))+
    1/2*rlnorm(1000,mu_F10-mu_F00-cov(log(F_10),log(F_00)),sqrt(var(log(F_00))+var(log(F_10))-2*cov(log(F_10),log(F_00))))
  vet_mean[k]=mean(model_F_11_0)
  
  
  model_F_11_0=rlnorm(1000,mu_F01-mu_F00+mu_F10+cov(log(F_01),log(F_00))+cov(log(F_10),log(F_00))-cov(log(F_10),log(F_01)),
                      sqrt(var(log(F_01))+var(log(F_00))+var(log(F_10))-2*cov(log(F_01),log(F_00))+2*cov(log(F_10),log(F_01))-2*cov(log(F_10),log(F_00))))
  vet_mean[k]= log(median(model_F_11_0))
  
}
