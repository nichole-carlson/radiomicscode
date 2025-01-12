f_clean <- function(candidates,oldkeeps,newkeep,full_cor,toler){
  howcorrelated <- abs(full_cor[candidates,newkeep]) # compute absolute correlation of candidates with newkeep
  candidates <- candidates[which(howcorrelated<toler)] # update candidates to be undercorrelated with newkeep
  # compute for each candidate the number of candidates (including itself) which are overcorrelated with it
  num_over_toler <- apply(abs(full_cor[candidates,candidates]),2,function(x){sum(x>=toler)})
  kees <- c(oldkeeps,newkeep,candidates[which(num_over_toler==1)]) # keep old, new, and anything which is only overcorrelated with itself
  cans <- candidates[which(num_over_toler!=1)] # update candidates to be undercorrelated with kees
  return(list(kees,cans))
}