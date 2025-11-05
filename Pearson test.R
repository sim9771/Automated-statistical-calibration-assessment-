pearson_test <- function(pred_all,subs,alp=0.05)
{ 
  O <- apply(subs,2,sum)
  E <- apply(pred_all,2,sum)
  d2 <- sum((O-E)^2/E)
  pval <- 1-pchisq(d2,df=(ncol(subs)-1))
  pval
}