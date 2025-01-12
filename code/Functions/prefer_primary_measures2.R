prefer_primary_measures <- function(candidates,keeps=NULL,full_cor,toler){
  # compute for each candidate the number of candidates (including itself) which are overcorrelated with it
  num_over_toler <- apply(abs(full_cor[candidates,candidates]),2,function(x){sum(x>=toler)}) 
  indices <- which(num_over_toler==max(num_over_toler)) # identify maximally overcorrelated
  first_order <- str_which(colnames(rad),"first_order") #
  first_order_indices <- intersect(indices,first_order)
  if(length(first_order_indices > 0)){
    newkeep_index <- sample(rep(first_order_indices,2),1)
  }
  else {
  newkeep_index <- sample(rep(indices,2),1) # pick something maximally overcorrelated to keep
  }
  return(candidates[newkeep_index])
}