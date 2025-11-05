
# The function expected_power outputs the expected power of the reliability test 
# of class probability estimates in multi-class classification.

# p is a vector containing class probability estimates, for which the user wants to test the reliability.
# The minimum dimension (i.e. length) of p is 3.

# k is a positive integer used for the reliability test.

# m is the sample size from which the expected power is computed.

# epsilon is the L1 distance between p and a sampled probability vector q. 
# epsilon must be small enough such that p-epsilon > 0 and p+epsilon < 1.

# alp is the significance level.

expected_power <- function(p,k,m=2000,epsilon=0.1,alp=0.05)
{
  exp_power <- numeric(m) 
  nclass <- length(p)
  ntick <- nclass-2 
  if(ntick==1) tick <- runif(m,min=-epsilon/2,max=epsilon/2) # when nclass = 3
  if(ntick>1)  tick <- matrix(runif(ntick*m,min=-epsilon/2,max=epsilon/2),nrow=m,ncol=ntick) # when nclass > 3
  B <- cbind(rep(-epsilon/2,m),tick,rep(0,m),rep(epsilon/2,m)) 
  u_mat <- t(apply(B,1,sort)) # each row gives u_(0),...,u_(nclass)
  
  v_mat <- matrix(0,nrow=m,ncol=nclass)
  for(j in 1:nclass)  # each row of v_mat gives v_0,...,v_nclass
  {
    v <- u_mat[,(j+1)] - u_mat[,j]
    idx <- which(u_mat[,j]<0)
    v[idx] <- -v[idx]  
    v_mat[,j] <- v
  }
  
  p_mat <- matrix(p,nrow=m,ncol=nclass,byrow=T)
  v_mat2 <- matrix(0,nrow=m,ncol=nclass)
  for(j in 1:nrow(v_mat))
  {
    idx <- sample(nclass,nclass)
    v_mat2[j,] <- v_mat[j,idx]
  }
  q_mat0 <- p_mat + v_mat2
  temp <- (q_mat0<0) + (q_mat0>1)
  idx <- which(apply(temp,1,sum)==0)
  q_mat <- q_mat0[idx,]  # each row of q_mat contains a random q such that |q-p| = epsilon
  
  crit <- qchisq(1-alp,df=nclass-1)
  for(j in 1:nrow(q_mat))
  {
    q <- q_mat[j,]
    U_1 <- sum(q * (q/p))
    U_2 <- sum(q * (q/p)^2)
    sigma2 <- U_2 - U_1^2
    W_q <- sum(k*(q-p)^2/p)
    out <- (crit-W_q)/(2*sqrt(k*sigma2))
    1-(pnorm(out))
    exp_power[j] <- 1-(pnorm(out)) # each element of exp_power gives the power of the test when q = true class probabilities
  }
  
  list(mean(exp_power[1:nrow(q_mat)]),sd(exp_power[1:nrow(q_mat)]))
}







