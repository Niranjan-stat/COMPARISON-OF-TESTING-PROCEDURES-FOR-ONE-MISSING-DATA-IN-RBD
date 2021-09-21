set.seed(50)

# Ordinary Approximation
F_sim = NULL
sim = 10^5
for(i in 1:sim)
{
  vec = numeric(0)
  v = 6 ; b = 5
  mu = 0 ;beta = seq(-(b-1)/2,(b-1)/2,length=5); tau = rep(0,v) 
  data_sim = mu + outer(beta,tau,"+") + matrix(rnorm(b*v),nrow = b)
  data_sim[1,1] <- 0
  B.1 <- apply(data_sim,MARGIN = 1,sum)
  T.1 <- apply(data_sim,MARGIN = 2,sum)
  G.1 <- sum(data_sim)
  x_est <- (b*B.1[1]+v*T.1[1]-G.1)/((b-1)*(v-1))
  data_sim[1,1] <- x_est
  B <- apply(data_sim,MARGIN = 1,sum)
  T <- apply(data_sim,MARGIN = 2,sum)
  G <- sum(data_sim)
  SSBl = sum(B^2)/v - (G^2/(b*v)) ; MSBl = SSBl/(b-1)
  SSTr = sum(T^2)/b - (G^2/(b*v)) ; MSTr = SSTr/(v-1)
  TSS = sum(data_sim^2) - (G^2/(b*v))
  SSE = TSS - SSBl - SSTr ; MSE = SSE/((b-1)*(v-1) - 1)
  F_obs = MSTr/MSE 
  F_sim[i] = F_obs
}
hist(F_sim,probability = TRUE,ylim = c(0,1),breaks = 60,xlim = c(0,6))
sup = seq(0,10,0.1)
lines(sup,df(sup,df1 = (v-1),df2 = (b-1)*(v-1) - 1))

F_ecdf = ecdf(F_sim)
plot(sup,F_ecdf(sup),type="l")
lines(sup,pf(sup,df1 = (v-1),df2 = (b-1)*(v-1) - 1),col = "red")

app_F=pf(sup,df1 = (v-1),df2 = (b-1)*(v-1) - 1)
sum(abs(F_ecdf(sup)-app_F))

#output
#> sum(abs(F_ecdf(sup)-app_F))
#[1] 0.5327011


library(ggplot2)
data=data.frame(hist(F_sim,probability = TRUE,ylim = c(0,1),breaks = 60,xlim = c(0,6)))
ggplot(data,aes(x=F_sim))+
geom_histogram()+
geom_density()

# Power Curve Codes

Power_Curve_1 <- function(tau,v = 6,b = 5,mu = 0,beta = seq(-(b-1)/2,(b-1)/2,length=5),sim = 10^5,alpha = 0.05)
{
  tau = tau-mean(tau)
  power.vec = NULL
  cut_pt = qf(p = alpha,df1 = (v-1),df2 = (b-1)*(v-1) - 1,lower.tail = FALSE)
  for(i in 1:sim)
  {
    data_sim = mu + outer(beta,tau,"+") + matrix(rnorm(b*v),nrow = b)
    data_sim[1,1] <- 0
    B.1 <- apply(data_sim,MARGIN = 1,sum)
    T.1 <- apply(data_sim,MARGIN = 2,sum)
    G.1 <- sum(data_sim)
    x_est <- (b*B.1[1]+v*T.1[1]-G.1)/((b-1)*(v-1))
    data_sim[1,1] <- x_est
    B <- apply(data_sim,MARGIN = 1,sum)
    T <- apply(data_sim,MARGIN = 2,sum)
    G <- sum(data_sim)
    SSBl = sum(B^2)/v - (G^2/(b*v)) ; MSBl = SSBl/(b-1)
    SSTr = sum(T^2)/b - (G^2/(b*v)) ; MSTr = SSTr/(v-1)
    TSS = sum(data_sim^2) - (G^2/(b*v))
    SSE = TSS - SSBl - SSTr ; MSE = SSE/((b-1)*(v-1) - 1)
    F_obs = MSTr/MSE 
    power.vec[i] = (F_obs > cut_pt)
  }
  return(mean(power.vec))
}

#plotting sum(tau_sq) and Power.val_1

v = 6
n = 20
sample.tau = NULL
tau_mat = NULL

for(i in 1:n)
{
  sample.tau = seq(-i/10,i/10,length = v)
  tau_mat = rbind(tau_mat,sample.tau -  mean(sample.tau))
}


tau_mat

sum_tau_sq = apply(tau_mat,MARGIN = 1,FUN = function(x){return(sum(x^2))})
Power.val = NULL

tau_mat[1,]

for(i in 1:n)
{
  Power.val[i] = Power_Curve_1(tau = tau_mat[i,],sim = 10^4)
}

Power.val
plot(sum_tau_sq,Power.val,type = "l")

size_0 = Power_Curve_1(tau = rep(0,v),sim = 10^4)

# OLPO Accurate Test 

F_sim_1 = NULL
sim = 10^5
for(i in 1:sim)
{
  v = 6 ; b = 5
  mu = 0 ; beta = seq(-(b-1)/2,(b-1)/2,length=5) ; tau = rep(0,v) 
  data_sim = mu + outer(beta,tau,"+") + matrix(rnorm(b*v),nrow = b)
  data_sim[1,1] <- 0
  B.1 <- apply(data_sim,MARGIN = 1,sum)
  T.1 <- apply(data_sim,MARGIN = 2,sum)
  G.1 <- sum(data_sim)
  x_est <- (b*B.1[1]+v*T.1[1]-G.1)/((b-1)*(v-1))
  data_sim[1,1] <- x_est
  B <- apply(data_sim,MARGIN = 1,sum)
  T <- apply(data_sim,MARGIN = 2,sum)
  G <- sum(data_sim)
  #calculated in augmented ANOVA
  SSBl = sum(B^2)/v - (G^2/(b*v)) ; MSBl = SSBl/(b-1)
  SSTr = sum(T^2)/b - (G^2/(b*v)) ; MSTr = SSTr/(v-1)
  TSS = sum(data_sim^2) - (G^2/(b*v))
  SSE = TSS - SSBl - SSTr ; MSE = SSE/((b-1)*(v-1) - 1)
  #calculated in unaugmented ANOVA
  SSBl.1 = (B.1[1]^2)/(v-1) + (sum(B[-1]^2)/v) - (G.1^2/(b*v - 1)) ; MSBl.1 = SSBl.1/(b-1)
  TSS.1 = sum(data_sim^2) - sum(data_sim[1,1]^2) - (G.1^2/(b*v - 1))
  SSTr.1 = TSS.1 - SSBl.1 - SSE ;  MSTr.1 = SSTr.1/(v-1)
  F_obs_1 = MSTr.1/MSE 
  F_sim_1[i] = F_obs_1
}
hist(F_sim_1,probability = TRUE,ylim = c(0,1),breaks = 60,xlim = c(0,6))
sup = seq(0,10,0.1)
lines(sup,df(sup,df1 = (v-1),df2 = (b-1)*(v-1) - 1))

F_ecdf_1 = ecdf(F_sim_1)
plot(sup,F_ecdf_1(sup),type="l")
lines(sup,pf(sup,df1 = (v-1),df2 = (b-1)*(v-1) - 1),col = "red")

app_F=pf(sup,df1 = (v-1),df2 = (b-1)*(v-1) - 1)
sum(abs(F_ecdf_1(sup)-app_F))


# Power Curve Codes (2)

Power_Curve_2 <- function(tau,v = 6,b = 5,mu = 0,beta = seq(-(b-1)/2,(b-1)/2,length=5),sim = 10^5,alpha = 0.05)
{
  tau = tau-mean(tau)
  power.vec = NULL
  cut_pt_1 = qf(p = alpha,df1 = (v-1),df2 = (b-1)*(v-1) - 1,lower.tail = FALSE)
  for(i in 1:sim)
  {
    data_sim = mu + outer(beta,tau,"+") + matrix(rnorm(b*v),nrow = b)
    data_sim[1,1] <- 0
    B.1 <- apply(data_sim,MARGIN = 1,sum)
    T.1 <- apply(data_sim,MARGIN = 2,sum)
    G.1 <- sum(data_sim)
    x_est <- (b*B.1[1]+v*T.1[1]-G.1)/((b-1)*(v-1))
    data_sim[1,1] <- x_est
    B <- apply(data_sim,MARGIN = 1,sum)
    T <- apply(data_sim,MARGIN = 2,sum)
    G <- sum(data_sim)
    #calculated in augmented ANOVA
    SSBl = sum(B^2)/v - (G^2/(b*v)) ; MSBl = SSBl/(b-1)
    SSTr = sum(T^2)/b - (G^2/(b*v)) ; MSTr = SSTr/(v-1)
    TSS = sum(data_sim^2) - (G^2/(b*v))
    SSE = TSS - SSBl - SSTr ; MSE = SSE/((b-1)*(v-1) - 1)
    #calculated in unaugmented ANOVA
    SSBl.1 = (B.1[1]^2)/(v-1) + (sum(B[-1]^2)/v) - (G.1^2/(b*v - 1)) ; MSBl.1 = SSBl.1/(b-1)
    TSS.1 = sum(data_sim^2) - sum(data_sim[1,1]^2) - (G.1^2/(b*v - 1))
    SSTr.1 = TSS.1 - SSBl.1 - SSE ;  MSTr.1 = SSTr.1/(v-1)
    F_obs_1 = MSTr.1/MSE 
    power.vec[i] = (F_obs_1 > cut_pt_1)
  }
  return(mean(power.vec))
}

Power_Curve_2(tau = sample(1:2,size = v,replace = TRUE),sim = 10^4)

#plotting sum(tau_sq) and Power.val_1

v = 6
n = 20
sample.tau = NULL
tau_mat = NULL

for(i in 1:n)
{
  sample.tau = seq(-i/10,i/10,length = v)
  tau_mat = rbind(tau_mat,sample.tau -  mean(sample.tau))
}


tau_mat
sum_tau_sq = apply(tau_mat,MARGIN = 1,FUN = function(x){return(sum(x^2))})
Power.val_1 = NULL
tau_mat[1,]


for(i in 1:n)
{
  Power.val_1[i] = Power_Curve_2(tau = tau_mat[i,],sim = 10^4)
}

Power.val_1
plot(sum_tau_sq,Power.val_1,type = "l")
lines(sum_tau_sq,Power.val,type = "l",lty = 2)

size_1 = Power_Curve_2(tau = rep(0,v),sim = 10^4)

# OTI Accurate Test 

F_sim_2 = NULL
sim = 10^5
for(i in 1:sim)
{
  vec = numeric(0)
  v = 6 ; b = 5
  mu = 0 ; beta = seq(-(b-1)/2,(b-1)/2,length=5) ; tau = rep(0,v) 
  data_sim = mu + outer(beta,tau,"+") + matrix(rnorm(b*v),nrow = b)
  data_sim[1,1] <- 0
  B.1 <- apply(data_sim,MARGIN = 1,sum)
  T.1 <- apply(data_sim,MARGIN = 2,sum)
  G.1 <- sum(data_sim)
  x_est_1 <- (b*B.1[1]+v*T.1[1]-G.1)/((b-1)*(v-1))
  data_sim[1,1] <- x_est_1
  B <- apply(data_sim,MARGIN = 1,sum)
  T <- apply(data_sim,MARGIN = 2,sum)
  G <- sum(data_sim)
  SSBl_1 = sum(B^2)/v - (G^2/(b*v)) 
  SSTr_1 = sum(T^2)/b - (G^2/(b*v)) 
  TSS_1 = sum(data_sim^2) - (G^2/(b*v))
  SSE_1 = TSS_1 - SSBl_1 - SSTr_1
  x_est_0 <-B.1[1]/(v-1)
  data_sim[1,1] <- x_est_0
  B <- apply(data_sim,MARGIN = 1,sum)
  T <- apply(data_sim,MARGIN = 2,sum)
  G <- sum(data_sim)
  SSBl_0 = sum(B^2)/v - (G^2/(b*v))
  TSS_0 = sum(data_sim^2) - (G^2/(b*v))
  SSE_0 = TSS_0 - SSBl_0
  F_obs = ((SSE_0 - SSE_1)/(v-1))/(SSE_1/(b*v-b-v))
  F_sim_2[i] = F_obs
}
hist(F_sim_2,probability = TRUE,ylim = c(0,1),breaks = 60,xlim = c(0,6))
sup = seq(0,10,0.1)
lines(sup,df(sup,df1 = (v-1),df2 = b*v-b-v))

F_ecdf_2 = ecdf(F_sim_2)
plot(sup,F_ecdf_2(sup),type = "l")
lines(sup,pf(sup,df1 = (v-1),df2 = b*v-b-v),col = "red")

app_F=pf(sup,df1 = (v-1),df2 = b*v-b-v)
sum(abs(F_ecdf_2(sup)-app_F))


# Power Curve Codes (2)

data.frame(Power.val_2,Power.val_1)


Power_Curve_3 <- function(tau,v = 6,b = 5,mu = 0,beta = seq(-(b-1)/2,(b-1)/2,length=5),sim = 10^4,alpha = 0.05)
{
  tau = tau-mean(tau)
  power.vec = NULL
  cut_pt_2 = qf(p = alpha,df1 = (v-1),df2 = (b*v-b-v),lower.tail = FALSE)
  for(i in 1:sim)
  {
    data_sim = mu + outer(beta,tau,"+") + matrix(rnorm(b*v),nrow = b)
    data_sim[1,1] <- 0
    B.1 <- apply(data_sim,MARGIN = 1,sum)
    T.1 <- apply(data_sim,MARGIN = 2,sum)
    G.1 <- sum(data_sim)
    x_est_1 <- (b*B.1[1]+v*T.1[1]-G.1)/((b-1)*(v-1))
    data_sim[1,1] <- x_est_1
    B <- apply(data_sim,MARGIN = 1,sum)
    T <- apply(data_sim,MARGIN = 2,sum)
    G <- sum(data_sim)
    SSBl_1 = sum(B^2)/v - (G^2/(b*v)) 
    SSTr_1 = sum(T^2)/b - (G^2/(b*v)) 
    TSS_1 = sum(data_sim^2) - (G^2/(b*v))
    SSE_1 = TSS_1 - SSBl_1 - SSTr_1
    x_est_0 <-B.1[1]/(v-1)
    data_sim[1,1] <- x_est_0
    B <- apply(data_sim,MARGIN = 1,sum)
    T <- apply(data_sim,MARGIN = 2,sum)
    G <- sum(data_sim)
    SSBl_0 = sum(B^2)/v - (G^2/(b*v))
    TSS_0 = sum(data_sim^2) - (G^2/(b*v))
    SSE_0 = TSS_0 - SSBl_0
    F_obs_2 = ((SSE_0 - SSE_1)/(v-1))/(SSE_1/(b*v-b-v))
    power.vec[i] = (F_obs_2 > cut_pt_2)
  }
  return(mean(power.vec))
}


v = 6
n = 50
sample.tau = NULL
tau_mat = NULL

for(i in 1:n)
{
  sample.tau = seq(-i/10,i/10,length = v)
  tau_mat = rbind(tau_mat,sample.tau -  mean(sample.tau))
}



sum_tau_sq = apply(tau_mat,MARGIN = 1,FUN = function(x){return(sum(x^2))})
Power.val_2 = NULL

Power_Curve_3(tau = sample(1:10,size = v,replace = TRUE),sim = 10^4)
qf(p = 0.05,df1 = (v-1),df2 = (b*v-b-v),lower.tail = FALSE)

for(i in 1:n)
{
  Power.val_2[i] = Power_Curve_3(tau = tau_mat[i,],sim = 10^3)
}

sim=10^4

Exact_power_2 <- function(s_t_s,v = 6,b = 5,mu = 0,beta = seq(-(b-1)/2,(b-1)/2,length=5),sim = 10^4,alpha = 0.05)
{
  qtile = seq(0,1,0.01)
  cut_pt_2 = qf(p = 0.05,df1 = (v-1),df2 = (b*v-b-v),lower.tail = FALSE)
  power.exact=NULL
  for(i in 1:length(qtile))
  {
    
    F_ext_2 = qf(qtile[i],df1=(v-1),df2=(b*v-b-v),ncp=s_t_s,lower.tail = FALSE)
    power.exact[i] = (F_ext_2 > cut_pt_2)
  }
  return(mean(power.exact))
}


power.ext_2=NULL
for(i in 1:n)
{
  power.ext_2[i] = Exact_power_2(s_t_s  = sum_tau_sq[i],sim = 10^3)
}

sup=seq(0,10,0.1)
plot(sum_tau_sq,power.ext_2,type="l")
plot(sup,df(sup,df1=(v-1),df2=(b*v-b-v),ncp=20),type="l")
lines(sup,df(sup,df1=(v-1),df2=(b*v-b-v),ncp=10),type="l")
lines(sup,df(sup,df1=(v-1),df2=(b*v-b-v),ncp=5),type="l")

qf(0.05,df1=(v-1),df2=(b*v-b-v),ncp=20,lower.tail = FALSE)


size_2 = Power_Curve_3(tau = rep(0,v),sim = 10^4)

#exact F calculation

plot(sum_tau_sq,Power.val_2,type = "l")
lines(sum_tau_sq,power.ext_2,type = "l")

#plot all power curves together
plot(sum_tau_sq,Power.val_2,type = "l")
lines(sum_tau_sq,Power.val_1,type = "l",lty = 2)
lines(sum_tau_sq,Power.val,type = "l",col="blue")

sum(Power.val>Power.val_1)
sum(Power.val_1>Power.val_2)

