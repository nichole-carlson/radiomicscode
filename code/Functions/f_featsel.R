f_featsel <- function(fullcor,tol,definitely_keep=c(),
                      rank_f=rank_f_greedy_then_random){
  keep_these <- which(1==apply(abs(fullcor),2,function(x){sum(x>=tol)})) # necessarily included
  keep_these <- unique(c(definitely_keep,keep_these)) # add desired inclusions
  if(length(keep_these)>0){ # define starting candidates based on starting keeps
    starting_cor <- apply(abs(fullcor[,keep_these,drop=FALSE]),1,max)
    candidate <- which(starting_cor<tol)
  }else{candidate <- 1:ncol(fullcor)}
  while(length(candidate)>0){
    keep_this_one <- rank_f(candidate,keep_these,fullcor,tol) # get new keeper using rank_f
    cl <- f_clean(candidate,keep_these,keep_this_one,fullcor,tol) # update keepers and candidates
    keep_these <- cl[[1]]
    candidate <- cl[[2]]}
  return(sort(keep_these))}


f_featsel_firstorder <- function(fullcor,tol,definitely_keep=c(),
                                 rank_f=greedy_then_random){
  firstorder <- str_which(colnames(fullcor),"first_order")
  fullcor_first <- fullcor[firstorder,firstorder]
  definitely_keep_first <- intersect(definitely_keep,firstorder)
  keep_these_first <- f_featsel(fullcor_first,tol,definitely_keep_first,rank_f)
  keep_these <- f_featsel(fullcor,tol,union(definitely_keep,firstorder[keep_these_first]),rank_f)
  return(sort(keep_these))}