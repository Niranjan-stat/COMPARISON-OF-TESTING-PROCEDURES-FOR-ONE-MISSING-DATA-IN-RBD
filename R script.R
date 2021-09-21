set.seed(50)
v = 6 ; b = 5 
mu = 0 ; beta = seq(-(b-1)/2,(b-1)/2,length=5); tau = rep(0,v)
sim = 10^4

F_sim = NULL; F_sim_1 = NULL; F_sim_2 = NULL 

for(i in 1:sim)
{
  data_sim = mu + outer(beta,tau,"+") + matrix(rnorm(b*v),nrow = b)
  data_sim[1,1] <- 0
  
  B.1_dash <- sum(data_sim[1,])          # B.1_dash is vector of B1'
  T.1_dash <- sum(data_sim[,1])          # T.1_dash is vector of T1'
  G_dash <- sum(data_sim)                # G_dash is vector of G'
  
  x_hat <- (b*B.1_dash+v*T.1_dash-G_dash)/((b-1)*(v-1))
  x_curl <- B.1_dash/(v-1)
  
                   #### approximate test ####
  
  data_sim[1,1] <- x_hat
  B <- rowSums(data_sim)
  T <- colSums(data_sim)
  G <- sum(data_sim)  
  SSBl = sum(B^2)/v - (G^2/(b*v))
  SSTr = sum(T^2)/b - (G^2/(b*v)) ; MSTr = SSTr/(v-1)
  TSS = sum(data_sim^2) - (G^2/(b*v)) 
  SSE = TSS - SSBl - SSTr ; MSE = SSE/(b*v-b-v)
  F_sim[i] <- MSTr/MSE
  
  sse_hat=SSE

            #### More accurate testing procedure-1  ####
  
  SSBl = (B.1_dash^2)/(v-1) + (sum(B[-1]^2)/v) - (G_dash^2/(b*v - 1)) 
  TSS = sum(data_sim^2) -sum(data_sim[1,1]^2) - (G_dash^2/(b*v - 1))
  SSTr = TSS - SSBl - sse_hat ;  MSTr = SSTr/(v-1)
  F_sim_1[i] <- MSTr/(sse_hat/(b*v-b-v))
   

            #### More accurate testing procedure-2  ####

  data_sim[1,1] <- x_curl
  B <- rowSums(data_sim)
  T <- colSums(data_sim)
  G <- sum(data_sim)
  
  SSBl = sum(B^2)/v - (G^2/(b*v))
  SSTr = sum(T^2)/b - (G^2/(b*v))
  TSS = sum(data_sim^2) - (G^2/(b*v))
  SSE0 = TSS - SSBl
  F_sim_2[i] <- ((SSE0 - SSE)/(v-1))/(sse_hat/(b*v-b-v))
}


sup = seq(0,6,0.1)
F_ecdf = ecdf(F_sim)
F_ecdf_1 = ecdf(F_sim_1)
F_ecdf_2 = ecdf(F_sim_2)
plot(sup,F_ecdf(sup),type = "l",lwd=0.5,lty=2,col="blue")
lines(sup,F_ecdf_1(sup),type = "l",lwd=0.5,col="red")
lines(sup,F_ecdf_2(sup),type = "l",lwd=0.5,lty=2,col="green")


