##Program Name: Lung_Image_Figure_helperfunctions.R
##Purpose: To create the plot of the various images from
##          the raw but masked imaging files.
###########################################################

f_pickid <- function(Order,Variable){ # automatically pick IDs
  rnm = c(Vble = Variable)
  rad <- dat %>%
    select(GRADSID,all_of(Variable)) %>%
    rename(all_of(rnm))
  if(Order=="min"){
    gid <- filter(rad,Vble == min(Vble))$GRADSID
  }else if(Order=="max"){
    gid <- filter(rad,Vble == max(Vble))$GRADSID
  }else if(Order=="median"){
    gid <- filter(rad,abs(Vble-median(Vble))==min(abs(Vble-median(Vble))))$GRADSID[1]
  }else{
    gid <- NA}
  return(gid)
}


clean_array <- # remove empty exterior slices of masked images
  function(arr){
    cleanarr <- arr
    for(d in 1:length(dim(arr))){
      zeros <- apply(cleanarr,d,function(subarr){sum(!is.na(subarr))!=0})
      cumsumzeros <- cumsum(zeros)
      cumsumzerosrev <- rev(cumsum(rev(zeros)))
      indies <- which(apply(cbind(cumsumzeros,cumsumzerosrev),1,min)>0)
      cleanarr <- 
        eval(parse(text=
                     paste0(
                       "cleanarr[",
                       paste(rep(",",d-1),collapse=""),
                       "c(",
                       paste(indies,collapse=","),
                       ")",
                       paste(rep(",",length(dim(arr))-d),collapse=""),
                       "]"
                     )))
    }
    return(cleanarr)
  }


select1000voxels <- # remove top and bottom slices with fewer than 1000 lung voxels
  function(arr){
    zeros <- apply(arr,3,function(subarr){sum(!is.na(subarr))>=1000})
    cumsumzeros <- cumsum(zeros)
    cumsumzerosrev <- rev(cumsum(rev(zeros)))
    indies <- which(apply(cbind(cumsumzeros,cumsumzerosrev),1,min)>0)
    return(arr[,,indies])
  }



f_readcropreduce <- function(gid,OldPipelineTF,path_img,
                             path_mask=NA,slice = NA){ # right lungs only
  if(OldPipelineTF){
    right <- readRDS(path_img)
  }else{
    right <- as.array(antsImageRead(path_img))
    mask <- as.array(antsImageRead(path_mask))
    right[mask!=1] <- NA
  }
  right <- right %>% 
    clean_array %>% 
    select1000voxels
  z <- slice
  if(is.na(z)){
    z = round(dim(right)[3]/2)
    print(paste0(gid,": selected slice ",z))}
  right <- right[,,z]
  return(right)
}



f_expanddimension <- function(arr,dim1,dim2){
  dim1_old <- nrow(arr)
  dim2_old <- ncol(arr)
  a1 <- round((dim1-dim1_old)/2)
  a2 <- dim1-dim1_old-a1
  b1 <- round((dim2-dim2_old)/2)
  b2 <- dim2-dim2_old-b1
  arr_new <- rbind(matrix(NA,nrow = a1,ncol=dim2_old),
                   arr,
                   matrix(NA,nrow = a2,ncol=dim2_old))
  arr_new <- cbind(matrix(NA,nrow = dim1,ncol=b1),
                   arr_new,
                   matrix(NA,nrow = dim1,ncol=b2))
  return(arr_new)
}



f_converttolongformat <- function(arr){
  arr_long <- cbind(1:nrow(arr),arr)
  colnames(arr_long) <- c("x",paste0("y",1:ncol(arr)))
  arr_long <- arr_long %>%
    as_tibble %>%
    pivot_longer(
      cols = paste0("y",1:ncol(arr)),
      names_to = "y",
      names_prefix = "y",
      values_to = "HU"
    ) %>%
    mutate(y = as.numeric(y))
  return(arr_long)
}



f_doAllTheThings <- function(df_pickids,OldPipelineTF){
  # read in all raw images, crop them down, 
  # and reduce to the axial slice of interest
  df_img_raw <- df_pickids %>% 
    apply(1,function(v){
      img <- f_readcropreduce(v["GRADSID"],OldPipelineTF,v["pathIMG"],v["pathMASK"],v["slice"])
      return(img)
    },simplify = F)
  # determine the largest image dimensions
  newdim <- df_img_raw %>%
    sapply(dim) %>%
    apply(1,max)
  # minimally expand all images out to have the same dimensions
  # then convert to long format in preparation for ggplot
  df_img_expand <- df_img_raw %>%
    lapply(function(arr){
      img <- arr %>%
        f_expanddimension(dim1=newdim[1],dim2=newdim[2]) %>%
        f_converttolongformat
    })
  # joint with df_pickids data and concatenate
  df_img <- 1:nrow(df_pickids) %>%
    lapply(function(indie){
      df_img_expand[[indie]] %>%
        mutate(Order = df_pickids$Order[indie],
               Variable = df_pickids$Variable[indie],
               Variable2 = df_pickids$Variable2[indie],
               GRADSID = df_pickids$GRADSID[indie])
    }) %>%
    do.call(bind_rows,.) %>%
    as_tibble()
  return(df_img)
}


















