########### Preparation ##########
library(nnet)
library(MASS)
library(combinat)
library(Hmisc)

ReadBasis<-function(BasisFun, DATA)
{
  tmp <- "group~"
  y1=c()
  y2=c()
  y3=c()
  group=c()
  x=read.DATA(DATA)
  for(i in 1:length(DATA)){
    group=c(group, rep(i-1,nrow(DATA[[i]])))
  }
  for(i in 1:length(BasisFun)){
    y1=c(y1,paste0("x",i))
    y2=c(y2,paste0(paste0("x",i),"=0"))
    y3=c(y3,paste0(paste0("x",i),"=",BasisFun[i],sep=""))
  }
  for(i in 1:length(BasisFun)){
    if(i!=length(BasisFun)){
      tmp <- paste(tmp,y1[i], "+",sep="")  
    }
    if(i==length(BasisFun)){
      tmp <- paste(tmp,y1[i],sep="")
    }
  }
  for(i in 1:length(BasisFun)){
    eval(parse(text = y2[i]))
  }
  for(i in 1:length(BasisFun)){
    eval(parse(text = y3[i]))
  }
  y1
  y2
  tmp
  out=multinom(as.formula(tmp),trace=F)
  out
}

read.DATA <- function(DATA)
{
  x = c()
  m = length(DATA)
  for(i in 1:m){
    x = c(x, DATA[[i]][,2])
  }
  x
}

quan.ecdf=function(tt, prob, quan.alpha)
{
  ## input is probability function, supports: tt; prob mass: prob.
  ## quan.alpha: the level of the desired quantile.
  st=sort(tt,index.return=T)
  st.value=st$x
  st.index=st$ix
  sortprob=prob[st.index]
  cdf=cumsum(sortprob)
  n=length(tt)
  if (cdf[1]<=quan.alpha)
  {
    i=which(cdf>quan.alpha)[1] ## make it smooth
    #st.value_1=st.value[i-1]+(st.value[i]-st.value[i-1])*(quan.alpha-cdf[i-1])/(cdf[i]-cdf[i-1])
    st.value_1=st.value[i]
  }
  else {st.value_1=NA}
  st.value_1
}

estmar=function(DATA, quan.alpha, BasisFun, Comp)  ## BasisFun: the choice of basis function 
{
  a=Comp[1]
  b=Comp[2]
  out.info=c()
  m=length(DATA)
  x=read.DATA(DATA)
  s=0
  for(i in 1:m)
  {
    s=s+nrow(DATA[[i]])
  }
  p=matrix(0,nrow=s,ncol=m)
  out=ReadBasis(BasisFun,DATA=DATA)
  hatG=fitted(out)
  if(m==2)
  {
    hatG=cbind(1-hatG,hatG)
    for(i in 1:m)
    {
      p[,i]=hatG[,i]/sum(hatG[,i])      # p[,i] is the empirical probability of ith data set
    }
  }
  if(!m==2)
  {
    for(i in 1:m)
    {
      p[,i]=hatG[,i]/sum(hatG[,i])
    }
  }
  out1=quan.ecdf(x,p[,a],quan.alpha)
  out2=quan.ecdf(x,p[,b],quan.alpha)
  out.info=c(out1,out2,out1-out2,-out$deviance/2,out$AIC)
  out.info
}

sboot=function(dat)
{
  index=dat[,1]
  ci=unique(index)                      # return a vector of the index of populations
  nci=length(ci)                        # amount of populations
  newci=sample(ci,nci,replace=T)
  newdata=c()
  for(i in 1:nci)
  {
    newdata=rbind(newdata,dat[(dat[,1]==newci[i]),])	
  }
  newdata
}

boots=function(dat)
{
  index=dat[,1]
  ci=unique(index)                      # return a vector of the index of populations
  nci=length(ci)                        # amount of populations
  newci=sample(ci,nci,replace=T)
  newdata=c()
  for(i in 1:nci)
  {
    newdata=rbind(newdata,dat[(dat[,1]==newci[i]),])	
  }
  newdata
}

###############################################
########## Data Generating Function ###########
###############################################

rmixn=function(d,mu,sig.epi,sig.gam) 
{
  ran=rnorm(1,0,sig.gam)
  err=rnorm(d,0,sig.epi)
  dat=mu+ran+err
  c(dat,ran)
}

ngendat=function(d,nc,mu,sig.epi,sig.gam,r)
{
  m <- length(mu) # the number of populations
  tmp1 <- c()
  tmp2 <- c()
  index <- c()
  group <- c() # the index of population
  data <- c() # data
  for (j in 1:nc)
  {
    data=rbind(data,rmixn(d,mu[1],sig.epi,sig.gam)[-(d+1)])
    tmp1=c(tmp1,rmixn(d,mu[1],sig.epi,sig.gam)[d+1])
    tmp1=as.numeric(tmp1)
    index=sample(1:nc,r,replace=F)
  }
  group=rep(1,nc)
  for (i in 2:m)
  {
    group=c(group,rep(i,nc))
    for (j in (1:nc)[-index])
    {
      data=rbind(data,rmixn(d,mu[i],sig.epi,sig.gam)[-(d+1)])
      tmp2=c(tmp2,rmixn(d,mu[i],sig.epi,sig.gam)[d+1])
    }
    for (j in index)
    {
      err=rnorm(d,0,sig.epi)
      dat=mu[i]+tmp1[j]+err
      data=rbind(data,dat)
    }
    tmp1=c(tmp2,tmp1[index])
    tmp1=as.numeric(tmp1)
    index=sample(1:nc,r,replace=F)
  }
  cbind(group,data)
}


####################################
########## Main Functions ##########
####################################

CEL_est=function(DATA, quan, BasisFun)
{
  m=length(DATA)
  x=read.DATA(DATA)
  out=0
  output=c()
  s=0
  rho <- c()
  for(i in 1:m)
  {
    s=s+nrow(DATA[[i]])
    rho <- c(rho,nrow(DATA[[i]]))
  }
  p=matrix(0,nrow=s,ncol=m)
  out=ReadBasis(BasisFun,DATA=DATA)
  hatG=fitted(out)
  theta_hat=summary(out)$coefficients
  if (identical(nrow(theta_hat),NULL))
  {
    theta_hat[1] <- theta_hat[1]-log(rho[-1]/rho[1])
    theta_hat <- as.matrix(theta_hat)
    rownames(theta_hat) <- NULL
    colnames(theta_hat) <- paste('theta_hat',2:m)   
  }
  else 
  {
    theta_hat[,1] <- theta_hat[,1]-log(rho[-1]/rho[1])
    colnames(theta_hat) <- NULL
    rownames(theta_hat) <- paste('theta_hat',2:m)   
  }
  if(m==2)
  {
    hatG=cbind(1-hatG,hatG)
    for(i in 1:m)
    {
      p[,i]=hatG[,i]/sum(hatG[,i])
    }
  }
  if(!m==2)
  {
    for(i in 1:m)
    {
      p[,i]=hatG[,i]/sum(hatG[,i])
    }
  }
  for (i in 1:m)
  {
    out=quan.ecdf(x,p[,i],quan)
    output=c(output,out)
  }
  cnames=quan
  rnames=c(1:m)
  result=list(quan_est=output,density=p,theta=theta_hat)
}

#### Simulation Functions

### gamma_DRM_1
drm_1=function(DATA, quantiles, BasisFun,Comp)
{
  DATA <- DATA
  BasisFun <- BasisFun
  Comp <- Comp
  a <- Comp[1]
  b <- Comp[2]
  quan <- 0.05
  output <- CEL_est(DATA,quan,BasisFun)
  dat <- c()
  for (i in 1:length(DATA))
  {
    dat=c(dat,DATA[[i]][,2])
  }
  dens1 <- output$density[,a]
  dens1[dens1 < 1e-100] <- 0 # prevent overflow
  dens1 <- dens1/min(dens1[dens1!=0])
  dens2 <- output$density[,b]
  dens2[dens2 < 1e-100] <- 0 # prevent overflow
  dens2 <- dens2/min(dens2[dens2!=0])
  quan1 <- wtd.quantile(dat,dens1,quantiles,type='i/(n+1)')
  quan2 <- wtd.quantile(dat,dens2,quantiles,type='i/(n+1)')
  gamma <- (sum(quan1>quan2)+1)/(length(quantiles)+2)
  return(gamma)
}

### gamma_DRM_2
drm_2 <- function(DATA,BasisFun,Comp)
{
  DATA <- DATA
  ld <- length(DATA)
  BasisFun <- BasisFun
  Comp <- Comp
  a <- Comp[1]
  b <- Comp[2]
  quan <- 0.05
  output <- CEL_est(DATA,quan,BasisFun)
  dat <- c()
  lc <- 0
  for (i in 1:ld)
  {
    lc <- lc+length(DATA[[i]][,2])
  }
  
  for (i in 1:length(DATA))
  {
    dat=c(dat,DATA[[i]][,2])
  }
  dat=sort(dat,index.return=TRUE)
  dat.index=dat$ix
  dens1 <- output$density[,a][dat.index]
  dens2 <- output$density[,b][dat.index]
  cdf1 <- cumsum(dens1)
  cdf2 <- cumsum(dens2)
  cdf1[lc] <- 1
  cdf2[lc] <- 1
  dcdf <- (cdf1-cdf2)<0
  tmp <- rle(dcdf)
  tmp1 <- tmp$lengths
  tmp2 <- tmp$values
  l <- length(tmp1)
  bindex <- 1
  eindex <- cumsum(tmp1)
  if (l>=2) 
  {
    for (i in 2:l)
    {
      bindex <- c(bindex,sum(tmp1[1:i])+1)
    }
    count <- rbind(tmp2,bindex,eindex)
    index <- as.matrix(count[c(2,3),count[1,]==1])
    ll <- ncol(index)
    gamma <- 0
    if (ll>=1) 
    {
      for (i in 1:ll)
      {
        gamma <- gamma + cdf2[index[2,i]] - cdf1[index[1,i]]
      }
      
      if (is.na(gamma)|identical(gamma,numeric(0))) gamma <- 0
    }
    
    if (ll==0) gamma <- 0
  }
  if (l==1) gamma <- 0
  return(gamma)
}

### empirical method
emp <- function(DATA,Comp)
{
  a <- Comp[1]
  b <- Comp[2]
  x <- DATA[[a]][,2]
  y <- DATA[[b]][,2]
  gamma <- mean(sort(x)>sort(y))
  return(gamma)
}
