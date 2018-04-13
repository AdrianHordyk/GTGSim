# Estimate VB parameters and bootstrap
VB <- function(Linf,k,t0,age){
  Linf*(1- exp(-k*(age-t0)))
}
