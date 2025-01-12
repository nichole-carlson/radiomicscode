my_fviz_gap_stat <- function (gap_stat, linecolor = "steelblue", maxSE = list(method = "firstSEmax", 
                                                          SE.factor = 1)) 
{
  if (!inherits(gap_stat, "clusGap")) 
    stop("Only an object of class clusGap is allowed. (cluster package)")
  if (is.list(maxSE)) {
    if (is.null(maxSE$method)) 
      maxSE$method = "firstmax"
    if (is.null(maxSE$SE.factor)) 
      maxSE$SE.factor = 1
  }
  else stop("The argument maxSE must be a list containing the parameters method and SE.factor")
  gap <- gap_stat$Tab[, "gap"]
  se <- gap_stat$Tab[, "SE.sim"]
  decr <- diff(gap) <= 0
  k <- .maxSE(gap, se, method = maxSE$method, SE.factor = maxSE$SE.factor)
  return(k)
  # df <- as.data.frame(gap_stat$Tab)
  # df$clusters <- as.factor(1:nrow(df))
  # df$ymin <- gap - se
  # df$ymax <- gap + se
  # p <- ggpubr::ggline(df, x = "clusters", y = "gap", group = 1, 
  #                     color = linecolor) + ggplot2::geom_errorbar(aes_string(ymin = "ymin", 
  #                                                                            ymax = "ymax"), width = 0.2, color = linecolor) + geom_vline(xintercept = k, 
  #                                                                                                                                         linetype = 2, color = linecolor) + labs(y = "Gap statistic (k)", 
  #                                                                                                                                                                                 x = "Number of clusters k", title = "Optimal number of clusters")
  # p
}


#  Determines the location of the maximum of f see ?cluster::maxSE
# +++++++++++++++++++++++++++++++++++++++++++
# f: numeric vector containing the gap statistic
# SE.f : standard error of the gap statistic
# method : character string indicating how the "optimal" number of clusters, k, 
# is computed from the gap statistics (and their standard deviations), 
# or more generally how the location k^ of the maximum of f[k] should be determined.
# SE.factor:  Determining the optimal number of clusters, Tibshirani et al. proposed the "1 S.E."-rule.
.maxSE <- function (f, SE.f, method = c("firstSEmax", "Tibs2001SEmax", 
                                        "globalSEmax", "firstmax", "globalmax"), SE.factor = 1) 
{
  method <- match.arg(method)
  stopifnot((K <- length(f)) >= 1, K == length(SE.f), SE.f >= 
              0, SE.factor >= 0)
  fSE <- SE.factor * SE.f
  switch(method, firstmax = {
    decr <- diff(f) <= 0
    if (any(decr)) which.max(decr) else K
  }, globalmax = {
    which.max(f)
  }, Tibs2001SEmax = {
    g.s <- f - fSE
    if (any(mp <- f[-K] >= g.s[-1])) which.max(mp) else K
  }, firstSEmax = {
    decr <- diff(f) <= 0
    nc <- if (any(decr)) which.max(decr) else K
    if (any(mp <- f[seq_len(nc - 1)] >= f[nc] - fSE[nc])) which(mp)[1] else nc
  }, globalSEmax = {
    nc <- which.max(f)
    if (any(mp <- f[seq_len(nc - 1)] >= f[nc] - fSE[nc])) which(mp)[1] else nc
  })
}
