#### Simulation Setting
library(readxl)
source('Functions.R')
# options(warn = -1)

### 1.Input of simulation functions

# For DRM method
quantiles <- seq(from=0.001,to=0.999,length.out=999) # quantiles to be compared for DRM method
BasisFun <- c('x','x^2','log(x)') # basis function for DRM method

# Other input
s <- 7 # dimension
B <- 1000 # bootstrap
n <- 50 # data length
gamma0 <- 0.5 # null hypothesis
cols_n_rows <- list( # colnames and rownames
  Row = c("SegRNN_96_96_AIT","ModernTCN_96_96_AIT","iTransformer_96_96_AIT","MTST_96_96_AIT","DLinear_96_96_AIT","NLinear_96_96_AIT","MTSMixer_96_96_AIT"), 
  Col = c("SegRNN_96_96_AIT","ModernTCN_96_96_AIT","iTransformer_96_96_AIT","MTST_96_96_AIT","DLinear_96_96_AIT","NLinear_96_96_AIT","MTSMixer_96_96_AIT")
)

### 2.Data Preparation

data <- read_excel("results_all.xlsx")

x0 <- as.vector(as.matrix(data['SegRNN_96_96_ms/sample']))
x1 <- as.vector(as.matrix(data['ModernTCN_96_96_ms/sample']))
x2 <- as.vector(as.matrix(data['iTransformer_96_96_ms/sample']))
x3 <- as.vector(as.matrix(data['MTST_96_96_ms/sample']))
x4 <- as.vector(as.matrix(data['DLinear_96_96_ms/sample']))
x5 <- as.vector(as.matrix(data['NLinear_96_96_ms/sample']))
x6 <- as.vector(as.matrix(data['MTSMixer_96_96_ms/sample']))

DATA <- list()
dat_1 <- cbind(c(1:n),x0)
dat_2 <- cbind(c(1:n),x1)
dat_3 <- cbind(c(1:n),x2)
dat_4 <- cbind(c(1:n),x3)
dat_5 <- cbind(c(1:n),x4)
dat_6 <- cbind(c(1:n),x5)
dat_7 <- cbind(c(1:n),x6)
DATA[[1]] <- dat_1
DATA[[2]] <- dat_2
DATA[[3]] <- dat_3
DATA[[4]] <- dat_4
DATA[[5]] <- dat_5
DATA[[6]] <- dat_6
DATA[[7]] <- dat_7


### 3.gamma_hat calculation

gamma_hat <- array(NA,dim=c(s,s))
dimnames(gamma_hat) <- cols_n_rows

for (i in 2:s){
  for (j in 1:(i-1)){
    Comp <- c(i,j) # the indice of populations to be compared for DRM method
    gamma_hat[i,j]<-drm_1(DATA,quantiles,BasisFun,Comp)
    # gamma_hat[i,j]<-drm_2(DATA,BasisFun,Comp)
  }
}


### 4.gamma_hat_star calculation

gamma_hat_star <- array(NA,dim=c(s,s,B))
dimnames(gamma_hat_star) <- cols_n_rows

for (k in 1:B)
{ 
  # resampling
  set.seed(k+1000)
  boot_DATA <- list()
  boot_dat_1 <- cbind(c(1:n),sample(x0, n, replace = TRUE))
  boot_dat_2 <- cbind(c(1:n),sample(x1, n, replace = TRUE))
  boot_dat_3 <- cbind(c(1:n),sample(x2, n, replace = TRUE))
  boot_dat_4 <- cbind(c(1:n),sample(x3, n, replace = TRUE))
  boot_dat_5 <- cbind(c(1:n),sample(x4, n, replace = TRUE))
  boot_dat_6 <- cbind(c(1:n),sample(x5, n, replace = TRUE))
  boot_dat_7 <- cbind(c(1:n),sample(x6, n, replace = TRUE))
  boot_DATA[[1]] <- boot_dat_1
  boot_DATA[[2]] <- boot_dat_2
  boot_DATA[[3]] <- boot_dat_3
  boot_DATA[[4]] <- boot_dat_4
  boot_DATA[[5]] <- boot_dat_5
  boot_DATA[[6]] <- boot_dat_6
  boot_DATA[[7]] <- boot_dat_7
  
  for (i in 2:s){
    for (j in 1:(i-1)){
      Comp <- c(i,j)
      gamma_hat_star[i,j,k]<-drm_1(boot_DATA,quantiles,BasisFun,Comp)
      # gamma_hat_star[i,j,k]<-drm_2(boot_DATA,BasisFun,Comp)
    }
  }
  
  if (k%%100 == 0){
    cat(sprintf("%d epochs have finished.\n", k))
  }
}


### 5.p value calculation

p_value <- array(NA, dim=c(s,s))
dimnames(p_value) <- cols_n_rows

for (i in 2:s){
  for (j in 1:(i-1)){
    temp <- mean(gamma_hat_star[i,j,] - gamma_hat[i,j] > gamma_hat[i,j] - gamma0)
    p_value[i,j] <- 2 * min(temp,1-temp)
  }
}


### 6.save results

write.csv(gamma_hat, file = "AIT_gamma_hat.csv")
write.csv(p_value, file = "AIT_p_value.csv")

