y_11=rnorm(100,0,1)
y_10=rnorm(100,0,1)
y_01=rnorm(100,2,1)
y_00=rnorm(100,2,1)

par(mfrow=c(1,4))
hist(y_11,prob = TRUE)
lines(density(y_11), col = 2, lwd = 2)
hist(y_10,prob = TRUE)
lines(density(y_10), col = 2, lwd = 2)
hist(y_01,prob = TRUE)
lines(density(y_01), col = 2, lwd = 2)
hist(y_00,prob = TRUE)
lines(density(y_00), col = 2, lwd = 2)

exp_y_11=exp(y_11)
exp_y_10=exp(y_10)
exp_y_01=exp(y_01)
exp_y_00=exp(y_00)

par(mfrow=c(1,4))
hist(exp_y_11,prob = TRUE,xlim = c(0,30))
lines(density(exp_y_11), col = 2, lwd = 2)
hist(exp_y_10,prob = TRUE,xlim = c(0,30))
lines(density(exp_y_10), col = 2, lwd = 2)
hist(exp_y_01,prob = TRUE,xlim = c(0,30))
lines(density(exp_y_01), col = 2, lwd = 2)
hist(exp_y_00,prob = TRUE,xlim = c(0,30))
lines(density(exp_y_00), col = 2, lwd = 2)



