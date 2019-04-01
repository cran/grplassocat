gen_dum <-
function(dat){
  x=list()
  for (i in seqWrapper(1,ncol(dat))){
    lvls=as.character(sort(unique(dat[,i])))
    cn=colnames(dat)[i]
    dums=data.frame(matrix(0, nrow=nrow(dat), ncol=length(lvls)))
    dat[,i]=as.character(dat[,i])
    colnames(dums)=paste(cn,lvls,sep=".")
    for (j in 1:nrow(dat)) {
      dums[j,paste(cn,dat[j,i],sep=".")]=1
    }
    x[[i]] <- dums
  }
  
  if(length(x)>0){return(x)} else {return(data.frame()[1:nrow(dat),])}
}
