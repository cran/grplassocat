library(grplasso)

seqWrapper = function(lb, ub, by=1) {
  s = c()
  if(!ub < lb) s = seq(lb,ub, by=by)
  return(s)
}

func_int=function(a){
  if(a[length(a)]!=0){a[1:(length(a)-1)]=-a[length(a)]} 
  return(a)}

gen_dum=function(dat){
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

fit_grp=function(eqn, dat, lambda, model=LinReg() , nonpen=c(), standardize = TRUE, ...){
  n=nrow(dat)
  
  any.nonpen=!is.null(nonpen)
  eqn_prev=eqn
  if(any.nonpen){
    eqn = as.formula(paste(c(deparse(eqn[[2]]), "~",
                            deparse(eqn[[3]]), "+",
                            deparse(nonpen[[length(nonpen)]])),
                          collapse = ""))
  }
  vars=all.vars(eqn)
  
  if(any(!(vars %in% colnames(dat))) & vars[2] !="."){stop("At least one variable from equation is not included in data")}
  
  a=names(Filter(is.factor, dat))
  
  inddp_vars=vars[2:length(vars)] #first variable in formula dependent variable (outcome)
  if(inddp_vars[1] =="."){ #if equation as outcome ~ . take all columns from data as features
    inddp_vars=colnames(dat)
    inddp_vars=inddp_vars[!inddp_vars==vars[1]]
    if(any.nonpen){
      trmlabels=attributes(terms(nonpen))$term.labels
    } else{
      trmlabels=c()
    }
  } else{
    trmlabels=attributes(terms(eqn))$term.labels
  }
 
  intcs=trmlabels[grepl(":",trmlabels)] #interactions
  
  fctrs=inddp_vars[inddp_vars %in% a] #variables that are factors (--> categorical variables)
  cnts=inddp_vars[!(inddp_vars %in% a)] #continuous variables
  

  transf=c() #transformed variables 
  for(j in seqWrapper(1,length(cnts))){
    transf=c(transf,trmlabels[grepl(cnts[j], trmlabels)])
  }
  transf=transf[!(transf %in% cnts)]
  transf=transf[!(transf %in% intcs)]
  eqn_tmp=sprintf("%s ~ %s ", vars[1], transf[1])
  for(j in seqWrapper(2,length(transf))){
    eqn_tmp=sprintf("%s + %s", eqn_tmp, transf[j])
  }
  
  for (j in seqWrapper(1, length(cnts))) { #allows for inclusion of transferred feat. only
    if(vars[2]=="."){ #case where equ. specified as outcome ~ .
      vrs_tmp=inddp_vars
    } else{
      vrs_tmp=as.character(attributes(terms(eqn))$variables) 
    }
    if(!(cnts[j]%in% vrs_tmp)) cnts=cnts[-j]
    
  }
  
  #create feature matrix including dummy variables for categorical variables
  if(length(transf)>0){
    dat_transf=model.matrix(as.formula(eqn_tmp), dat)
    tab=cbind.data.frame(dat[,cnts, drop=FALSE],dat_transf[,-1, drop=FALSE],gen_dum(as.data.frame(dat[,fctrs, drop=FALSE]))) 
  } else {
    tab=cbind.data.frame(dat[,cnts, drop=FALSE],gen_dum(as.data.frame(dat[,fctrs, drop = FALSE]))) 
  }

  #check which interaction terms between cat + cat or cat + cont variables 
  intcs_cnt=c()
  intcs_tmp=intcs
  for(j in seqWrapper(1,length(intcs))){ 
    vrs=unlist(strsplit(intcs[j], ":"))
    if(any(vrs %in% c(cnts,transf))){
      intcs_cnt=c(intcs_cnt, intcs[j])
      intcs_tmp=intcs_tmp[!intcs_tmp==intcs[j]]
    }
  }
  intcs=intcs_tmp
  
  #check with interactions between two cont. variables 
  intcs_tcnt=c()
  intcs_tmp=intcs_cnt
  for(j in seqWrapper(1,length(intcs_cnt))){ 
    vrs=unlist(strsplit(intcs_cnt[j], ":"))
    if(all(vrs %in% c(cnts,transf))){
      intcs_tcnt=c(intcs_tcnt, intcs_cnt[j])
      intcs_tmp=intcs_tmp[!intcs_tmp==intcs_cnt[j]]
    }
  }
  intcs_cnt=intcs_tmp

  nr_cont=length(cnts)+length(transf) + length(intcs_tcnt) #number continuous variables
  nr_cat=length(fctrs) #number categorical variables
  
  #generate vector with labels for all variables
  nr_lvls=rep(NA, nr_cat) #levels of each categorical variable
  for (j in seqWrapper(1,nr_cat)) {
    #nr_lvls[j]=nlevels(dat[,fctrs[j]])
    nr_lvls[j]=length(unique(dat[,fctrs[j]]))
  }
  
  lbls=c(cnts, transf)
  lbls2=lbls #for non-pen list
  lbls_prev=c(lbls,  rep(fctrs, times=nr_lvls))
  
  #interactions between two continuous variables
  x_cont=data.frame(tab[,c(cnts,transf), drop=FALSE])  
  for (j in seqWrapper(1, length(intcs_tcnt))) {
    cvrs=unlist(strsplit(intcs_tcnt[j], ":"))
    x_cont=cbind.data.frame(x_cont, tab[,cvrs[1]] * tab[,cvrs[2]])
    colnames(x_cont)[ncol(x_cont)]=paste(cvrs[1], cvrs[2], sep="_")
    lbls=c(lbls, paste(cvrs[1], cvrs[2], sep="_"))
    lbls2=c(lbls2, paste(cvrs[2], cvrs[1], sep="_"))
  }
  
  #standardization of continuous variables and interactions between cont. vars
  if(standardize){
    mtr=colMeans(x_cont)
    norce=rep(NA,ncol(x_cont))
    for(j in seqWrapper(1,ncol(x_cont))){
      x_cont[,j] = x_cont[,j] - mean(x_cont[,j])
      norce[j]=sqrt(sum(x_cont[,j]^2))
      x_cont[,j] = (x_cont[,j]/norce[j])
    }
  } else {
    norce=rep(1,nr_cont)
    mtr=rep(0,nr_cont)
    if(nr_cont>=1) warning("Continuous variables not standardized")
  }
  scs=mtr/norce #for re-scaling intercept
  x=x_cont
  
  lbls=c(lbls, rep(fctrs, times=nr_lvls))
  lbls2=c(lbls2, rep(fctrs, times=nr_lvls))
  
  #standardize categorical variables
  ns=list() #frequencies for each categorical variable
  for (j in seqWrapper(1,nr_cat)) {
    Xdummy=tab[,which(lbls_prev==fctrs[j])]
    standX <- apply(Xdummy, 2, function(z) (z-mean(z))/sqrt(sum(z)))
    x=cbind.data.frame(x, standX)
    ns=c(ns,list(colSums(Xdummy)))
    scs=c(scs, colMeans(Xdummy))
  }
  
  #generate matrices for interaction terms with categorical vars. + prepare for group lasso fit
  n_ints=list()
  svs_resc=list() #for re-sacling of coefficients based on interaction terms
  lbls_intc=list()
  lvls_ints=c()  #number of levels for interaction terms
  Ls=Ms=c()
  scs_bet=list()
  cnames_intc=c() #names of coefficients for later 
  for (k in seqWrapper(1,length(intcs))) {
    cvrs=unlist(strsplit(intcs[k], ":"))
    X1dummy=tab[,lbls_prev==cvrs[1]]
    X2dummy=tab[,lbls_prev==cvrs[2]]
    
    L=ncol(X1dummy)
    M=ncol(X2dummy)
    X12dummy = matrix(nrow = n, ncol = L*M)
    for(i in 1:L){
      for(j in 1:M){
        X12dummy[,(i-1)*M+j] = X1dummy[,i]*X2dummy[,j]
      }
    }
    n12 = colSums(X12dummy)
    
    
    # start preparing for group lasso fit
    
    svd_main = svd(cbind(X1dummy, X2dummy))
    U = svd_main$u[,1:(L+M-1)]
    
    standX12 = apply(X12dummy, 2, function(z) z/sqrt(sum(z)))
    standX12 = standX12 - U %*% (t(U) %*% standX12)
    
    colnames(standX12)=t(outer(colnames(X1dummy), colnames(X2dummy), FUN = "paste" ))
    
    x=cbind.data.frame(x, standX12)
    
    lbls=c(lbls,rep(paste(cvrs[1], cvrs[2], sep="_" ), ncol(standX12)))
    lbls2=c(lbls2,rep(paste(cvrs[2], cvrs[1], sep="_" ), ncol(standX12)))
    
    n_ints=c(n_ints, list(n12))
    
    resc=t(t(svd_main$v[,1:(L+M-1)]) * (1/svd_main$d[1:(L+M-1)])) %*% t(U) %*% X12dummy 
    svs_resc=c(svs_resc, list(resc))
    lbls_intc=c(lbls_intc, list(c(cvrs[1], cvrs[2])))
    
    lvls_ints=c(lvls_ints, length(n12))
    Ls=c(Ls,L)
    Ms=c(Ms,M)
    
    cnames_intc=c(cnames_intc, colnames(standX12))
    
  }
  
  #names of coefficients for later
  cnames=colnames(x)
  
  #generate matrices for interaction terms with cont. vars. + prepare for group lasso fit
  #save information for re-scaling of coefficients
  us=list()
  vs=list()
  svs=list()
  
  for (k in seqWrapper(1,length(intcs_cnt))) {
    cvrs=unlist(strsplit(intcs_cnt[k], ":"))
   
    Xint=t(apply(tab, 1, function(x) c(as.vector(t(outer(as.numeric(x[lbls_prev==cvrs[1]]),as.numeric(x[lbls_prev==cvrs[2]])))))))
    catvr=ifelse(cvrs[1] %in% cnts, 2, 1)
    convr=ifelse(cvrs[1] %in% cnts, 1, 2)
    colnames(Xint)=as.character(t(outer(colnames(tab)[lbls_prev==cvrs[convr]], colnames(tab)[lbls_prev==cvrs[catvr]], FUN = "paste")))
    
    I=t(apply(Xint, 1, function(x) func_int(x)))
    I=I[,1:(ncol(I)-1)]
    
    PI=apply(data.frame(I), 2, function(x) x-mean(x))
    svd_cn=svd(as.matrix(PI))
    
    u_tmp=svd_cn$u
    us=c(us, list(u_tmp))
    vs=c(vs, list(svd_cn$v))
    svs=c(svs, list(svd_cn$d))
    
    scs=c(scs, colMeans(data.frame(I))[colMeans(data.frame(I))!=0])
    
    colnames(u_tmp)=colnames(Xint)[1:(ncol(Xint)-1)]
    x=cbind.data.frame(x,u_tmp)
    
    lbls=c(lbls,rep(paste(cvrs[1], cvrs[2], sep="_" ), ncol(u_tmp)))
    lbls2=c(lbls2,rep(paste(cvrs[2], cvrs[1], sep="_" ), ncol(u_tmp)))
    cnames=c(cnames, colnames(Xint))
  }
  
  nzcols=!is.nan(colSums(x)) #exclude zero-columnss
  ncfs_all=colnames(x) #names of coefficients (later NA for zero-columns)
  x_all=x
  nr_all=ncol(x_all)
  x=x[,nzcols]
  lbls_all=lbls
  lbls=lbls[nzcols]
  lbls2=lbls2[nzcols]
  
  groups=c() #indices for group lasso fit
  for (i in 1:length(unique(lbls))) {
    ll=unique(lbls)[i]
    groups=c(groups, rep(i, length(which(lbls==ll))))
    
  }

  #groups=c(as.integer(unclass(factor(lbls)))) # 
  if(length(nonpen)>0){ #NAs for non-penalized variables
    vrs_nonpen=gsub(":", "_" ,attributes(terms(nonpen))$term.labels)
    groups[lbls %in% unique(vrs_nonpen)]=NA
    groups[lbls2 %in% unique(vrs_nonpen)]=NA
  }
  
  #print(groups)

  resp=vars[1] #dependent variable

  incl_intc=FALSE

  if(model@name=="Logistic Regression Model"){
    if(is.numeric(dat[,resp])){
      y=dat[,resp]
    } else{
      outc=as.character(sort(unique(dat[,resp]))[1])
      y=rep(0,nrow(dat))
      y[which(dat[,resp]==outc)]=1
    }
    phi=as.matrix(cbind(x,1))
    groups=c(groups, NA)
    incl_intc=TRUE
    nr_all=nr_all+1
  } else if (model@name=="Poisson Regression Model"){
    y=dat[,resp]
    phi=as.matrix(cbind(x,1))
    groups=c(groups, NA)
    incl_intc=TRUE
    nr_all=nr_all+1
  } else{ #Linear Regression
    y=dat[,resp]
    y=y-mean(y)
    phi=as.matrix(x)
  }

  #group lasso fit
  mod<- grplasso(phi, y, index = groups, weights = rep(1, length(y)),
                 offset = rep(0, length(y)), lambda = lambda,
                 penscale = sqrt, model = model, center = FALSE,
                 standardize = FALSE, ...)

  #----rescale coefficients for continuous variables and transfer coefficients for categorical ones----
  
  betahat=rep(0, nr_all)
  betahat[nzcols] = coef(mod)
  betahat_vect=betahat
  
  #compute linear predictors before re-scaling (for testing)
  lps=phi %*% coef(mod)
  
  #stack coefficients into list
  betahat=list()
  ind=0
  lbl_un=unique(lbls_all)
  tab_r = count=sapply(lbl_un,function(x)sum(lbls_all==x))
  for (j in 1:length(unique(lbls_all))) {
    inds=(ind+1):(ind+tab_r[j])
    #print(inds)
    betahat[j]=list(betahat_vect[inds])
    ind=ind+(tab_r[j])
  }

  #continuous variables
  coef_cnts=betahat_vect[seqWrapper(1,nr_cont)]/norce
  cfs_intc=betahat_vect[seqWrapper(1,nr_cont)] #for re-scaling intercept

  #categorical variables
  coef_cat=c()
  ind=nr_cont+1
  for (j in seqWrapper(1,nr_cat)) {
    bhat_cat=betahat[ind]
    n1=ns[j]
    theta1 = unlist(bhat_cat) / sqrt(unlist(n1))
    theta1 = theta1 - mean(theta1)
    coef_cat=c(coef_cat, theta1)
    cfs_intc=c(cfs_intc, theta1)
    ind=ind+1
  }

  #interaction terms (categorical variables)
  coef_int=c()
  deltas=list()
  for (j in seqWrapper(1, length(intcs))) {
    theta12 = unlist(betahat[ind]) / sqrt(unlist(n_ints[j]))
    indsNaN=is.nan(theta12)
    #theta12[indsNaN]=0
    
    L=Ls[j]
    M=Ms[j]
    theta12 = matrix(theta12, nrow = L, ncol = M, byrow = TRUE)

    while(any(abs(c(rowSums(theta12, na.rm = TRUE), colSums(theta12, na.rm = TRUE))) > 1E-12)){
      theta12 = theta12 - matrix(rep(rowMeans(theta12, na.rm = TRUE), each = M), nrow = L, ncol = M, byrow = TRUE)
      theta12 = theta12 - matrix(rep(colMeans(theta12, na.rm = TRUE), each = L), nrow = L, ncol = M)
    }

    theta12 = as.numeric(t(theta12))
    theta12[indsNaN]=NA
    ind=ind+1
    coef_int=c(coef_int, theta12)
    
    #save information to adjust coefficients of main effects
    theta12_tmp=theta12
    theta12_tmp[is.na(theta12_tmp)]=0
    delta_tmp=matrix(unlist(svs_resc[j]), nrow=L+M) %*% theta12_tmp
    delta1=list(delta_tmp[1:L])
    delta2=list(delta_tmp[(L+1):(L+M)])
    names(delta1)=unlist(lbls_intc[j])[1]
    names(delta2)=unlist(lbls_intc[j])[2]
    deltas=c(deltas,delta1, delta2)
  }
  
  #adjust coefficients for main effects of interaction terms
  corr_terms=0 #adjustment of intercept after re-scaling
  for (j in seqWrapper(1, length(deltas))) {
    inds_cat=which(lbls[lbls %in% fctrs] == names(deltas)[j])
    coef_cat[inds_cat]=coef_cat[inds_cat]-unlist(deltas[j])
    #center again
    corr_terms=corr_terms+mean(coef_cat[inds_cat])
    coef_cat[inds_cat]=coef_cat[inds_cat]-mean(coef_cat[inds_cat])
  }
  
  #interaction terms (categorical + continuous variables)
  coef_int_cont=c()
  for (j in seqWrapper(1, length(intcs_cnt))) {
    alpha=unlist(betahat[ind])/unlist(svs[j])
    beta_int=matrix(unlist(vs[j]), ncol=length(alpha))%*%alpha
    cfs_intc=c(cfs_intc, beta_int)
    coef_int_cont=c(coef_int_cont,beta_int,-sum(beta_int))
    ind=ind+1
    }
  
  #intercept 
  if(incl_intc) {
    intc = betahat_vect[length(betahat_vect)]
    intc = intc-sum(scs * cfs_intc)+corr_terms
  } else{intc=-sum(scs * cfs_intc)+corr_terms}
  
  coefs=data.frame(c(intc,coef_cnts, coef_cat, coef_int, coef_int_cont))
  rownames(coefs)=c("Intc", cnames)
  #if(incl_intc){rownames(coefs)=c("Intc", cnames)}else{rownames(coefs)=cnames}
  
  
  if(model@name=="Linear Regression Model"){
    coefs["Intc",]=coefs["Intc",]+mean(dat[,resp])
  }
  
  colnames(coefs)="Coefficient"
  
  return(coefs)
}