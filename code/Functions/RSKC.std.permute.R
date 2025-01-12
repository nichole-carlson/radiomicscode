RSKC.std.permute <- function (x, 
                          K = NULL, # Number of clusters, can be a range to test
                          nstart = 25, # Number of random starts for the clustering method
                          iter.max = 10, # Number of max iterations for kmeans
                          nperms = 25, # Number of permutations to find the best Gap statistic
                          alpha = NULL,
                          wbounds = NULL, # This is only for certain clustering methods
                          nvals = 10, # Not sure what this is for
                          silent = TRUE
)
{
  x <- as.matrix(x)
  
  # Calculate the objective function for the data
  if(silent == FALSE){
    message('Calculating the WBSS from the data')
  }
  tots <- NULL 
  for(wi in 1:length(wbounds)){
    out <- RSKC(d = x, ncl = K, alpha = alpha, L1 = wbounds[wi], nstart = nstart, silent = silent)
    bcss <- max(out$WBSS)
    tots <- c(tots, bcss)
  }
  
  # Permute the data
  if(silent == FALSE){
    message('Calculating the WBSS from the permuted data')
  }
  permtots <- matrix(NA, nrow = length(wbounds), ncol = nperms)
  for (i in 1:nperms) {
    if(silent == FALSE){
      message(paste0('Permutation #', i))
    }
    permx <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (j in 1:ncol(x)) permx[, j] <- sample(x[, j])
    for(wi in 1:length(wbounds)){
      out <- RSKC(d = permx, ncl = K, alpha = alpha, L1 = wbounds[wi], nstart = nstart, silent = silent)
      bcss <- max(out$WBSS)
      permtots[wi, i] <- bcss
    }
  }
  
  # Calculate the Gap statistic
  gaps <- (log(tots) - apply(log(permtots), 1, mean))/sd(log(permtots))
  out <- list(tots = tots, 
              permtots = permtots, 
              gaps = gaps,  
              bestw = wbounds[which.max(gaps)])
  return(out)
}