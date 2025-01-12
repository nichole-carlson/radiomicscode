myGap <- function (x, 
                   method = 'kmeans', # Can also take options of 'KMeansSparseCluster' and 'RSKC'
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
  for(ki in 1:length(K)){
    
    if(method == 'kmeans'){  
     out <- kmeans(x = x, centers = K[ki], nstart = nstart, iter.max = iter.max) 
     bcss <- out$betweenss
     tots <- c(tots, bcss)
    }
    if(method == 'KMeansSparseCluster'){
      out <- KMeansSparseCluster(x = x, K = K[ki], wbounds = wbounds[ki], nstart = nstart, silent = silent)
      bcss <- GetWCSS(x, out[[1]]$Cs)$bcss.perfeature
      tots <- c(tots, sum(out[[1]]$ws * bcss))
    }
    if(method == 'RSKC'){
      out <- RSKC(d = x, ncl = K[ki], alpha = alpha, L1 = wbounds[ki], nstart = nstart, silent = silent)
      bcss <- max(out$WBSS)
      tots <- c(tots, bcss)
    }
  }
  
  # Permute the data
  if(silent == FALSE){
    message('Calculating the WBSS from the permuted data')
  }
  permtots <- matrix(NA, nrow = length(K), ncol = nperms)
  for (i in 1:nperms) {
    if(silent == FALSE){
      message(paste0('Permutation #', i))
    }
    permx <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (j in 1:ncol(x)) permx[, j] <- sample(x[, j])
    for(ki in 1:length(K)){
      
      if(method == 'kmeans'){  
        out <- kmeans(x = permx, centers = K[ki], nstart = nstart, iter.max = iter.max) 
        bcss <- out$betweenss
        permtots[ki, i] <- bcss
      }
      if(method == 'KMeansSparseCluster'){
        out <- KMeansSparseCluster(x = permx, K = K[ki], wbounds = wbounds[ki], nstart = nstart, silent = silent)
        bcss <- GetWCSS(permx, out[[1]]$Cs)$bcss.perfeature
        permtots[ki, i] <- sum(out[[1]]$ws * bcss)
      }
      if(method == 'RSKC'){
        out <- RSKC(d = permx, ncl = K[ki], alpha = alpha, L1 = wbounds[ki], nstart = nstart, silent = silent)
        bcss <- max(out$WBSS)
        permtots[ki, i] <- bcss
      }
    }
  }
  
  # Calculate the Gap statistic
  gaps <- (log(tots) - apply(log(permtots), 1, mean))
  out <- list(tots = tots, 
              permtots = permtots, 
              gaps = gaps,  
              bestk = K[which.max(gaps)])
  return(out)
}


# this file contains the hidden functions for the sparcl package (used for the KMeansSparseCluster)
GetWCSS <- function(x, Cs, ws=NULL){
  wcss.perfeature <- numeric(ncol(x))
  for(k in unique(Cs)){
    whichers <- (Cs==k)
    if(sum(whichers)>1) wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers,],center=TRUE, scale=FALSE)^2, 2, sum)
  }
  bcss.perfeature <- apply(scale(x, center=TRUE, scale=FALSE)^2, 2, sum)-wcss.perfeature
  if(!is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), wcss.ws=sum(wcss.perfeature*ws),
                               bcss.perfeature=bcss.perfeature))
  if(is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), bcss.perfeature=bcss.perfeature))
}

