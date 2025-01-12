myRSKC  <- function (d, ncl, alpha, L1 = 12, nstart = 200, silent = TRUE, 
          scaling = FALSE, correlation = FALSE) 
{
  if (alpha > 1 | alpha < 0) 
    stop("alpha must be between 0 and 1")
  if (!is.null(L1)) 
    if (L1 < 1) 
      stop("L1 value must be greater or equal to 1 or NULL!")
  if (is.data.frame(d)) 
    d <- as.matrix(d)
  r.ncl <- round(ncl)
  if (ncl != r.ncl) 
    ncl <- r.ncl
  if (ncl <= 1) 
    stop("ncl must be positive integer > 1! but ncl=", ncl)
  if (scaling) 
    d = scale(d)
  if (correlation) 
    d = t(scale(t(d)))
  if (is.null(L1)) 
    sparse <- FALSE
  else {
    sparse <- TRUE
  }
  n <- nrow(d)
  Nout <- floor(n * alpha)
  f <- ncol(d)
  g <- f + 1
  W <- rep(1, f)
  sumW <- f
  if (sum(is.na(d)) == 0) {
    miss <- FALSE
    if (sparse) {
      Result <- RSKC.a1.a2.b(d, L1, ncl, nstart, alpha, 
                             n, f, g, Nout, silent)
    }
    else {
      Result <- RSKC.trimkmeans(d, ncl, trim = alpha, runs = nstart, 
                                maxit = 10000)
    }
  }
  else {
    d[is.nan(d)] <- NA
    miss <- TRUE
    if (sparse) {
      Result <- RSKC.a1.a2.b.missing(d, L1, ncl, nstart, 
                                     alpha, n, f, g, Nout, silent)
    }
    else {
      Result <- RSKC.trimkmeans.missing(d = d, ncl = ncl, 
                                        w = W, trim = alpha, runs = nstart, points = Inf, 
                                        maxit = 10000)
    }
  }
  if (sparse) {
    Result$oW <- sort(Result$oW)
    if (Nout == 0) {
      Result$oW <- Result$oE <- "undefined"
    }
  }
  else {
    Result <- modified.result.nonsparse(Result, ncl, f)
    if (Nout == 0) {
      Result <- modified.result.kmean(Result)
    }
  }
  # Result$disttom <- Result$ropt <- Result$trim <- Result$scaling <- Result$centers <- Result$criterion <- Result$classification <- Result$means <- Result$ropt <- Result$k <- Result$runs <- NULL
  if (!is.null(colnames(d))) 
    names(Result$weights) <- colnames(d)
  Input <- list(N = n, p = f, ncl = ncl, L1 = L1, nstart = nstart, 
                alpha = alpha, scaling = scaling, correlation = correlation, 
                missing = miss)
  r2 <- c(Input, Result)
  class(r2) <- "rskc"
  return(r2)
}

RSKC.a1.a2.b <-
  function(d,L,ncl,nstart,alpha,n,f,g=f+1,Nout,silent=TRUE)
  {
    iteration<-10;
    ## The robust sparse K-means algorithm terminates when ( newss - oldss ) / newss  < e^-8 at i^th step (b) and 
    ## return C from the i^th step(a) and weights from (i-1)^th step (b)
    ## must be Nout=floor(alpha*n)
    oldss<-(-Inf)
    WBSS<-rep(NA,iteration)
    W<-rep(1/sqrt(f),f)
    for (i in 1:iteration){	
      ## Step a
      w.data<-t(t(d[,W!=0,drop=FALSE])*sqrt(W[W!=0])) # reduce the dimention so that the alg is more efficient
      reduced.f<-sum(W!=0)
      
      ## A1<-RSKC.trimkmeans(w.data,ncl,trim=alpha,runs=nstart,points=NULL,maxit=10000) ## modi Dec 29
      A1<-RSKC.trimkmeans(w.data,ncl,trim=alpha,runs=nstart,points=Inf,maxit=10000)
      
      ## if (Nout!=0){out.a1<-which(A1$classification==(ncl+1))}else{out.a1<-n+1} ## modi Dec 29
      if (Nout!=0){
        out.a1<-A1$oW
      }else{
        out.a1<-n+1
      }
      
      if (i!=1&oldss!=-Inf){
        Ab4<-RSKC.trimkmeans(w.data,ncl,trim=alpha,runs=1,
                             points=cent(w.data[-save.out.a1,,drop=FALSE],C[-save.out.a1],
                                         ncl,reduced.f),maxit=10000)
        
        ##if (Nout!=0) out.a1.b4<-which(Ab4$classification==(ncl+1)) else{out.a1.b4<-n+1} ## modi Dec 29
        if (Nout!=0) out.a1.b4 <- Ab4$oW else{out.a1.b4<-n+1}
        
        ##bfsum<-sum(Ab4$disttom[-out.a1.b4]) ## modi Dec 29
        b4sum<-Ab4$WSS
        
        ##afsum<-sum(A1$disttom[-out.a1]) ## modi Dec 29
        afsum<-A1$WSS
        
        if (!silent) cat("step a: old wwss",b4sum," new wwss",afsum,"\n")	       
        if (b4sum<afsum){A1<-Ab4;out.a1<-out.a1.b4}
      }
      w.mu<-A1$means; # not necessarily ncl by p mostliely ncl by p' ; p'<<p
      
      ## Assign cluster labels to the outliers
      ## new.C<-A1$classification
      ## if(Nout!=0) new.C<-class.trimk(w.data,w.mu,new.C,ncl,Nout) ## modi Dec 29
      
      new.C<-A1$labels
      now.mu<-cent(d[-out.a1,,drop=FALSE],new.C[-out.a1],ncl,f)
      
      if (Nout==0){
        out.a2<-n+1 # skip (a-2) jump to (b)
      }else{
        ## Step a-2 (trimming)
        out.a2<-RSKC.step.a2(d,new.C,ncl,Nout,n,f,now.mu)
        n00<-n-length(union(out.a1,out.a2));
      }
      # If the reduction of out.a2 makes only one cluster nonempty
      # repeat step a and a-2 again from the begining.
      C0<-new.C[-c(out.a2,out.a1)]
      if(length(unique(C0))==1){cat("alpha is too large for n \n"); next;}
      save.out.a1<-out.a1;save.out.a2<-out.a2;C<-new.C
      
      # Step b
      # perform without out.a2 nor out.a1
      B<-RSKC.step.b(L,cbind(d[-c(out.a2,out.a1),],C0),ncl,W,f,n00)
      newss <- B$WBSS 
      
      if(i!=1&oldss!=-Inf){
        diff<-( newss - oldss ) / newss 
        if(!silent){p.n<-round(newss,2);p.d<-round(diff,2)
        cat("step b: wbss",p.n,"wbss dif",p.d,"\n")}
        if( diff < 1e-8 ) break;
      }
      oldss <- newss
      W<-B$We;WBSS[i]<-B$WBSS
    }
    if (i == iteration) cat("Step (a), (a2) and (b) are repeated over the maximum number of iterations. The algorithm might not converge.")
    return(list(centers=now.mu,oE=save.out.a2,oW=save.out.a1,labels=C,weights=W,WBSS=WBSS[1:(i-1)]))
  }

cent<-function(d,C,ncl,f){
  mu<-matrix(NA,ncl,f)
  for ( k in 1 :ncl){
    mu[k,]<-colMeans(d[C==k,,drop=FALSE])
  }
  return(mu)
}


class.trimk<-function(w.data,mu,trimC,ncl,Nout){
  out.WdisC<-matrix(NA,Nout,ncl)
  out.cases<-w.data[trimC==(ncl+1),,drop=FALSE]
  for (ncl.i in 1 : ncl)  
  {out.WdisC[,ncl.i]<-rowSums(scale(out.cases,center=mu[ncl.i,],scale=FALSE)^2)}
  out.belong.index<-max.col(-out.WdisC)
  trimC[trimC==(ncl+1)]<-out.belong.index
  return(trimC)
}	

RSKC.step.a2 <-function(d,gr,ncl,Nouttr,n,f,center){
  # Step (a-2) (trimming) 
  g<-f+1
  dtr<-cbind(d,gr,1:n)
  whoisout<-RSKC.step.a2.trim.Euclidean(dtr,g,ncl,n,center)
  if (Nouttr<1){otr<-n+1}
  else{
    o<-order(whoisout[,1],decreasing=TRUE)[1:Nouttr]
    otr<-whoisout[o,2]
  }
  return(sort(otr))
}


RSKC.step.a2.trim.Euclidean<-function(d2,g,ncl,n,mu)
{ # d2 is n by (p+2)
  old<-1;whoisout<-matrix(nrow=n,ncol=2);f<-(g-1);ind<-d2[,g]
  for ( k in 1:ncl)
  {
    dff<-NULL;indW<-NULL;df<-d2[ind==k,-g,drop=FALSE] # df nk by p+1
    dff<-df[,-g,drop=FALSE] # dff nk by p
    Wmu<-mu[k,]
    nk<-nrow(dff);indW<-rowSums(!is.na(dff))
    di<-rowSums(scale(dff,center=Wmu,scale=FALSE)^2,na.rm=TRUE)*f/indW # nk by 1
    new<-nk+old-1
    whoisout[old:new,1:2]<-cbind(di,df[,g])
    old<-new+1
  }  
  return(whoisout)
}

RSKC.step.b<-
  function(L,d,ncl,W,f,N00){ # d is indexed data
    a<-RSKC.step.b.BSS_js(d,(f+1),ncl,W,N00)
    a<-pmax(a,0)
    de<-BinarySearch.delta(tss.wss=a,L1=L)
    # FINALLY obatin the weights which maximize the weighted MBCSS
    SS<-pmax(a-de,0) # p by 1
    Sn2<-norm2(SS)
    W<-SS/Sn2
    WBSS<-sum(W*a)
    # the position of inf corresponds to the
    # of itr where the ok element corrumption occur 
    return(list(We=W,WBSS=WBSS))
  }


RSKC.step.b.BSS_js<-
  function(dout,g,ncl,W,n){ 
    a<-matrix(NA,nrow=ncl,ncol=g-1);sumW<-sum(W)
    d<-dout[,-g];p<-g-1;
    
    Tmu<-colMeans(d,na.rm=TRUE);Tterm<-scale(d,center=Tmu,scale=FALSE)^2 # n by p
    # unless all entries of a column is NA, this Tmu vector has no missing value 
    
    index<-dout[,g] # p by 1
    for ( k in 1:ncl){ 
      if (sum(index==k)==0){
        ##cat("no obs in sm clusters");
        next;
      }else{
        df<-d[index==k,,drop=FALSE] #nk by p
        # replace from Rcenter ?
        Wmu<-colMeans(df,na.rm=TRUE) # p by 1
        
        nk<-nrow(df)
        Wterm<-scale(df,center=Wmu,scale=FALSE)^2 #nk by p
        Ttermk<-Tterm[index==k,,drop=FALSE] # nk by p
        
        T_W<-Ttermk-Wterm # nk by p
        
        adjust <- as.vector(sumW / ( (!is.na(Wterm)) %*% W ))
        # nk by 1 contains the scaler multiple of each obs
        adjustTerm<-T_W*adjust #nk by p
        # each row is multiplied by corresponding entry of adjust vec
        a[k,]<-colSums(adjustTerm,na.rm=TRUE) #p by 1
      }
    }
    a<-colSums(a,na.rm=TRUE) 
    # in case sm clus contains no obs, na.rm=TRUE
    return(a)
  }



BinarySearch.delta <- function(tss.wss,L1){
  # L1=l1bound
  # tss.wss= tssperfeatre - wssperfeature
  if(norm2(tss.wss)==0 || sum(abs(tss.wss/norm2(tss.wss)))<=L1) return(0)
  lam1 <- 0
  lam2 <- max(abs(tss.wss))-1e-5
  iter <- 1
  while(iter<=15 && (lam2-lam1)>(1e-4)){
    su <- soft(tss.wss,(lam1+lam2)/2)
    if(sum(abs(su/norm2(su)))<L1){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    iter <- iter+1
  }
  return((lam1+lam2)/2)
}

soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}


CER <-
  function(ind, true.ind,nob=length(ind))
  {
    if (length(ind) != length(true.ind))
      stop("length(ind) != length(true.ind)",length(ind) ,"!=", length(true.ind))
    ind.num <- ind
    true.ind.num <- true.ind
    if (is.character(ind)){
      uniInd <- unique(ind)
      for (i in 1 : length(uniInd))
      {
        ind.num[ind == uniInd[i]] <- i
      }
    }
    
    if (is.character(true.ind) ){
      uniInd <- unique(true.ind)
      for (i in 1 : length(uniInd))
      {
        true.ind.num[true.ind == uniInd[i]] <- i
      }
      
    }
    
    rand <- randIndex(as.numeric(ind.num),as.numeric(true.ind.num),correct=FALSE)
    names(rand) <- NULL
    return(1 - rand)
    ##return(sum(abs(one(true.ind)-one(ind)))/choose(nob,2))
  }

## one <-
## function(index){
## 	on <- NULL 
## 	c <- combn(index,2)
## 	c <- t(c)
## 	on <- 1*(c[,1]==c[,2])
## 	return(on)
## 	}

norm1 <-
  function(y){sum(abs(y))}

norm2 <-
  function(y){sqrt(sum(y^2))}


## function specific for the opt digits
## generate bitmap of given observation
showbitmap <-function(index)
{
  bitmapMat = NULL## This line is to avoid the complain about No visible binding from R CMD CHECK 
  bitmapLab = NULL## This line is to avoid the complain about No visible binding from R CMD CHECK
  data(bitmapMat,envir = environment()) ## lazyloading
  data(bitmapLab,envir = environment()) ## lazyloading
  Nbit=32
  for (iindex in 1 : length(index))
  {
    indivindex<- index[iindex]
    for (ibit in 1 : Nbit)
    {
      cat("\n")
      cat(bitmapMat[[indivindex]][ibit])
    }
    cat("\n obs=",indivindex," true digit=",bitmapLab[indivindex]," \n")
  }
}


## function declaration

showDigit <- function(index,cex.main=1)
{
  DutchUtility = NULL## This line is to avoid the complain about No visible binding from R CMD CHECK
  data(DutchUtility,envir = environment()) ## lazyloading
  ## 4. DutchUtility-pix: 240 pixel averages in 2 x 3 windows; 
  ## 16 by 15
  ncols = 15 ## replace to 15 
  nrows = 16 ## replace to 16
  labels <- rep(0:9,each=200)
  plot(NA,xlim=c(0,ncols),ylim=c(0,nrows),axes=FALSE,xlab="",ylab="",cex.main=cex.main,
       main=paste("observation",index," True digit",labels[index],sep=""))
  abline(h=0:ncols,v=0:nrows,lty=2,col="gray70")
  axis(1,(1:ncols)-0.5,1:ncols,lty=0,cex=0.5)
  axis(2,(1:nrows)-0.5,nrows:1,padj=1,lty=0,cex=0.5)
  ##
  cols <- gray.colors(n=6,start=0.9,end=0.3)
  for (iy in 1 : nrows)
  {
    for (ix in 1 : ncols)
    {
      ## each row vector of a matrix DutchUtility contains the bitmap X; 
      ## DutchUtility[1,] = c(X[1,],X[2,],...) 
      Pickedcolor <- cols[DutchUtility[index,ix+ncols*(iy-1)]]
      polygon(x=c(ix-1,ix-1,ix,ix),
              y=c(nrows-iy,nrows-iy+1,nrows-iy+1,nrows-iy),
              col=Pickedcolor,
              border=FALSE)
    }
  }
}



## sensitivity the other way around for the first example
## sensitivity for the second example both way

Sensitivity <- function (label1, label2
                         ##, Alpha=FALSE, which.Alpha="label1"
)
{
  ## Given two partitions, label1 and label2, compute its "sensitivity"
  ## max{ (label1[label2==k])/sum(label2==k) }, and the corresponding label of label1 that achieves max 
  uni1 <- sort(unique(label1))
  uni2 <- sort(unique(label2))
  
  ## if (length(uni1) >= 27 & Alpha & which.Alpha=="label1")
  ##   stop("The alphabet option is valid only if the number of clusters from label1 is less than 27")
  ## if (length(uni2) >= 27 & Alpha & which.Alpha=="label2")
  ##   stop("The alphabet option is valid only if the number of clusters from label2 is less than 27")
  
  ## ## change the labels of clusters in which.Alpha
  ## if (Alpha & which.Alpha=="label1"){
  ##   temp <- label1
  ##   for (iuniL in 1 : length(uni1))
  ##     {
  ##       temp[label1==uni1[iuniL]] <- letters[iuniL]
  ##     }
  ##   label1 <- temp
  ##   uni1 <- sort(unique(label1))
  
  ## }else if (Alpha & which.Alpha=="label2"){
  ##   temp <- label2
  ##   for (iuniL in 1 : length(uni1))
  ##     {
  ##       temp[label2==uni2[iuniL]] <- letters[iuniL]
  ##     }
  ##   label2 <- temp
  ##   uni2 <- sort(unique(label2))
  ## }
  
  K <- length(uni2)
  senst <- correspondClass <- rep(NA, K)
  ## this will return table both columns and rows are in sorted order.
  tbl <- table(label1, label2)
  
  trueTot <- colSums(tbl)
  ## scale each column of the table matrix by its column sum
  ##  colSums(prMat)= 1,1,...
  prMat <- scale(tbl, scale = trueTot, center = FALSE)
  senst <- apply(prMat, 2, max)*100
  
  correspondClass <- rownames(prMat)[apply(prMat, 2, which.max)]
  senst <- sprintf("%1.0f",senst)
  
  re <- data.frame(rbind(senst, correspondClass))
  
  names(re) <- uni2
  rownames(re) <- c("Sensitivity. (%)", "Class label by label1.")
  return(list(prob=re,table=tbl,marginal=prMat))
}
