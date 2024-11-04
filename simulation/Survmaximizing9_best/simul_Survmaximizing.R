###----------- Load Libraries ---------------###
library(dplyr)
library(tidyverse)
library(Rglpk) 
library(glue)
library(survival)
library(doParallel)


######---------------------------------------------------------------------#####
######                     Simulation Settings                             #####
######---------------------------------------------------------------------#####

## SLURM Batch array id
SLURMid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Number of simulations per each array
Simul_Num = 10


######---------------------------------------------------------------------#####
######                   Data generation function                          #####
######---------------------------------------------------------------------#####

## Step 0. Cluster size (N) distribution
N.dist = function(m) sample(x = c(5,10,15), size = m, replace = T)

## Step1. Covariate generation
X.gen <- function(N){
  X1 <- rnorm(N, 0, 3)
  X2 <- runif(N, 0, 1)
  X3 <- rbinom(N, 1, 0.5)
  X = data.frame(X1, X2, X3)
  return(X)
}

## All covariates
X.names = c("X1", "X2", "X3")

## Step2. Treatment model (experiment setting)
p.A = function(X,N) with(X, plogis(-1 + 0.5*X1 - 0.5*X2))

## Covariates for A model
X.A.names = c("X1", "X2")

A.gen <- function(X, N){
  A <- rbinom(N, 1, p.A(X,N))
  return(A)
}

## Step 3. T model (uniform distribution)
upper_T <- function(j, A, X, N){
  with(X,
       1/(2*abs(X1[j]) + 1*X2[j] + 
            (1*X1[j] - 0.5*X2[j])*A[j] + 
            6*1/(N-1)*( sum((X3)*A) - (X3[j])*A[j] ))
  )
}

## R(T) function - defines value function

### Time of interest
tau = 0.3

R <- function(t){
  return(I(t > tau))
}

## E[R(T)|A,X,N] function
## For R(t) = I(t > tau), E[R(T)|A,X,N] = S^T(tau|A,X,N) = max(0, 1 - tau/u^T(A,X,N))
mu_RT <- function(j, A, X, N){
  u = upper_T(j, A, X, N)
  return(max(0, 1 - tau / u))
}

T.gen <- function(A, X, N){
  T <- sapply(1:N, function(j) runif(1, 0, upper_T(j, A, X, N)))
}

### beta.Oracle from indirect learning - by maximizing \sum_j mu_RT(j,A,X,N)
beta.Oracle.init = c(0, -1, 0.5, -6)

## Step 4. C model (exponential distribution)
C.gen <- function(A, X, N){
  lambda_C <- with(X, 0.5*A + 0.5*X1 + 0.5*X2)
  C <- rexp(N, exp(lambda_C))     
}

## Covariates for C model
X.C.names = c("X1", "X2")

## Step 4. Cluster data generation
DGP <- function(N){
  
  ## Covariate generation
  X = X.gen(N)
  
  ## Treatment model
  A = A.gen(X, N)
  
  ## T model
  T_ = T.gen(A, X, N)
  
  ## C model
  C = C.gen(A, X, N)
  
  ## Combine data
  Y <- pmin(T_,C)
  D <- as.numeric(T_ <= C)
  
  data = cbind(data.frame(Y = Y, T = T_, C = C, D = D, A = A), X)
  
  return(data)
  
}

## Step 5. All data generation
data.sim = function(m){
  N.list = N.dist(m)
  data = lapply(N.list, DGP)
  data = dplyr::bind_rows(data,.id = "id") %>% mutate(id = as.numeric(id))
  return(data)
}


### Example ###
m = 100
data = data.sim(m)
summary(data)

######---------------------------------------------------------------------#####
######                         Help functions                              #####
######---------------------------------------------------------------------#####


## Normalize vector function
normalize_vector <- function(vec) {
  
  # Compute the second norm (Euclidean norm)
  norm <- sqrt(sum(vec^2))
  
  # Check to avoid division by zero
  if (norm == 0) return(vec)
  
  # Normalize the vector
  normalized_vec <- round(vec / norm, 2)
  
  return(normalized_vec)
}

###--- A (treatment) function fitting ---###
#' @input Training data O_train = {(Y, D, A, X, N)}_train
#' @output censoring function fit object and training time points (for stieltjes integral)
#' @example A.fit = A.train(data_train, X.A.names)
A.train = function(data_train, X.A.names){
  
  A.fit <- glm(A ~ ., data = data_train %>% select(all_of(c("A", X.A.names))), family = "binomial")
  
  return(A.fit)
}
  
###--- G (censoring) function fitting ---###
#' @input Training data O_train = {(Y, D, A, X, N)}_train
#' @output censoring function fit object
#' @example G.fit = G.train(data_train, X.C.names)
G.train = function(data_train, X.C.names){
  
  ## Fit Cox PH regression
  G.fit <- coxph(Surv(Y,1-D) ~ ., 
                 data = data_train %>% select(all_of(c("Y","D","A", X.C.names))))
  
  return(G.fit)
}

###--- IPCW.i function ---###
#' @input cluster level data O.i, nuisance estimators G.fit, 
#' @output IPCW.i = {D_ij / G_ij(Y_ij)}_{j=1}^{N_i}

IPCW = function(O.i, G.fit){
  
  Y.i = O.i$Y
  D.i = O.i$D
  N.i = nrow(O.i)
  
  G.i.Y.i = rep(1, N.i)
  
  for(j in 1:N.i){
    
    # Use survfit to compute the survival curve for j-th individual
    # and extract survival probabilities at time Y_ij, using the 'extend' argument
    G.i.Y.i[j] <- summary(
      survfit(G.fit, newdata = O.i[j,]), 
      times = Y.i[j], 
      extend = TRUE)$surv
    
  }
  
  IPCW.i = D.i / G.i.Y.i
  
  return(IPCW.i)
  
}

###--- Mixed Integer Linear Programming function ---###
#' @input w: weight of MIP formulation, X_all: all covariates (including 1 column for intercept)
#' @output beta: optimizer

solveMILP <- function(w, X_all){
  
  ## Total observations
  nn_all = nrow(X_all) 
  
  ## dimension of beta
  d = ncol(X_all)
  
  ## Objective function coefficients for (beta, p_{ij}) 
  # First d elements are coefficients of beta in the objective function, which are 0 (since they do not appear)
  # and remaining elements are coefficients of p_{ij}, which is \overline{Y}_i * (Aij / ej(1|Xi) - (1-Aij) / ej(0|Xi))
  f <- c(rep(0,d), w)  
  
  ## bounds on coefficients, i.e., all |beta_k| <= B = 1'
  B <- 1  
  
  ## maximum values of x'beta
  C <- B * apply(abs(X_all), 1, sum)  
  
  ## prevent non-integer numbers the integrality constraint of integers from being counted as integers
  minmargin <- max(1, C) * (1e-8)  
  
  ## Constraints matrix A (LHS)
  Aineq <- rbind(cbind(-X_all, diag(C)), cbind(X_all, -diag(C)) )
  
  ## Direction of constraints
  dir =  c(rep("<=", 0.5*nrow(Aineq)),rep("<", 0.5*nrow(Aineq))) 
  
  ## Constraints constant b (RHS)
  bineq <- c(C - minmargin, 0*C-minmargin)
  
  ## First d elements are lower bound on \beta_k values, and remaining are lower bound on p_{ij} (which is 0)
  lb <- c(-B * rep(1,d), rep(0, nn_all))  
  
  ## First d elements are upper bound on \beta_k values, and remaining are upper bound on p_{ij} (which is 1)
  ub <- c(B * rep(1, d), rep(1, nn_all))  
  
  bounds= list(lower = list(ind = 1:length(lb), val =lb),
               upper = list(ind = 1:length(ub), val =ub))
  
  ## Variable type string: "C" for continuous (\beta_k) and "B" for binary (p_ij)
  ctype <- c(rep("C", d), rep("B", nn_all))
  
  ## Mixed Integer Programming
  sol <- Rglpk_solve_LP(obj = f, mat = Aineq, dir = dir, rhs = bineq, bounds = bounds, types = ctype,max = T)
  
  ## First d elements are beta
  beta = normalize_vector(sol$solution[c(1:d)])
  
  return(beta)
  
}


######---------------------------------------------------------------------#####
######                    Policy Learning Methods                          #####
######---------------------------------------------------------------------#####


### Scenario indicates 
### (i) KN: Known PS, No Censoring
### (ii) UN: Unknown PS (estimated), No Censoring
### (iii) UC: Unknown PS (estimated), Censoring

###----------- Method 1: Additive IPW ---------------###

learning.addIPW <- function(data){
  
  ## Number of clusters
  m = length(unique(data$id))
  
  ## Propensity score estimation
  A.fit <- A.train(data, X.A.names)
  
  ## C model estimation (Cox PH regression)
  G.fit <- G.train(data, X.C.names)
  
  ## Compute optimization weight at each cluster
  w.addIPW = function(i){
    
    O.i = data %>% filter(id == i)
    Y.i = O.i$Y
    T.i = O.i$T
    A.i = O.i$A
    N.i = nrow(O.i)
    
    ### PS predict
    p.A.fit.i = predict(A.fit, O.i, type = "response")
    
    ### IPCW
    IPCW.fit.i = IPCW(O.i, G.fit)
    
    ### (i) Known PS, No Censoring
    p.A.i = p.A(O.i, N.i)
    w.i.KN = mean( R(T.i) )*(A.i/p.A.i - (1-A.i)/(1-p.A.i))
    
    ### (ii) Unknown PS (estimated), No Censoring
    w.i.UN = mean( R(T.i) )*(A.i/p.A.fit.i - (1-A.i)/(1-p.A.fit.i))
    
    ### (iii) Unknown PS (estimated), Censoring
    w.i.UC = mean( IPCW.fit.i*R(Y.i) )*(A.i/p.A.fit.i - (1-A.i)/(1-p.A.fit.i))
    
    # Combine weights into a data frame for each id
    w.i <- data.frame(KN = w.i.KN, UN = w.i.UN, UC = w.i.UC)
    
    return(w.i)
    
  }
  
  w_all <- do.call(rbind, lapply(1:m, w.addIPW))
  row.names(w_all) <- NULL  # Reset row names
  
  ## 1 is for intercept 
  X_all = cbind(1, data %>% select(all_of(X.names)))  
  
  beta.addIPW = list()
  
  ### (i) Known PS, No Censoring
  beta.addIPW[["KN"]] = solveMILP(w_all$KN, X_all)
  
  ### (ii) Unknown PS (estimated), No Censoring
  beta.addIPW[["UN"]] = solveMILP(w_all$UN, X_all)
  
  # Combine weights into a data frame for each id
  beta.addIPW[["UC"]] = solveMILP(w_all$UC, X_all)
  
  return(beta.addIPW)
  
}

### Example ###
# learning.addIPW(data)


###----------- Method 2: No Interference IPW ---------------###

learning.NoIntIPW <- function(data){
  
  ## Number of clusters
  m = length(unique(data$id))
  
  ## Propensity score estimation
  A.fit <- A.train(data, X.A.names)
  
  ## C model estimation (Cox PH regression)
  G.fit <- G.train(data, X.C.names)
  
  ## Compute optimization weight at each cluster
  w.NoInt = function(i){
    
    O.i = data %>% filter(id == i)
    Y.i = O.i$Y
    T.i = O.i$T
    A.i = O.i$A
    N.i = nrow(O.i)
    
    ### PS predict
    p.A.fit.i = predict(A.fit, O.i, type = "response")
    
    ### IPCW
    IPCW.fit.i = IPCW(O.i, G.fit)
    
    ### (i) Known PS, No Censoring
    p.A.i = p.A(O.i, N.i)
    w.i.KN = 1/N.i * R(T.i) *(A.i/p.A.i - (1-A.i)/(1-p.A.i))
    
    ### (ii) Unknown PS (estimated), No Censoring
    w.i.UN = 1/N.i * R(T.i) *(A.i/p.A.fit.i - (1-A.i)/(1-p.A.fit.i))
    
    ### (iii) Unknown PS (estimated), Censoring
    w.i.UC = 1/N.i * IPCW.fit.i * R(Y.i) *(A.i/p.A.fit.i - (1-A.i)/(1-p.A.fit.i))
    
    # Combine weights into a data frame for each id
    w.i <- data.frame(KN = w.i.KN, UN = w.i.UN, UC = w.i.UC)
    
    return(w.i)
    
  }
  
  w_all <- do.call(rbind, lapply(1:m, w.NoInt))
  row.names(w_all) <- NULL  # Reset row names
  
  ## 1 is for intercept 
  X_all = cbind(1, data %>% select(all_of(X.names)))  
  
  beta.NoInt = list()
  
  ### (i) Known PS, No Censoring
  beta.NoInt[["KN"]] = solveMILP(w_all$KN, X_all)
  
  ### (ii) Unknown PS (estimated), No Censoring
  beta.NoInt[["UN"]] = solveMILP(w_all$UN, X_all)
  
  # Combine weights into a data frame for each id
  beta.NoInt[["UC"]] = solveMILP(w_all$UC, X_all)
  
  return(beta.NoInt)
  
}

### Example ###
# learning.NoIntIPW(data)

###----------- Method 3: Cluster IPW ---------------###

# Optimization with classical cluster-IPW value estimator involves smooth stochastic approximation
# at the computation of objective_fun 
## Define objective function under classical IPW value estimator since optim can only minimize

obj_fun_IPW = function(beta, data, p.A.list, IPCW.list) { 
  
  m = length(unique(data$id))  # number of clusters
  
  dev_wts=NULL
  
  for(i in 1:m){
    
    O.i = data %>% filter(id == i)
    Y.i = O.i$Y
    A.i = O.i$A
    X.i = O.i %>% select(all_of(X.names))
    N.i = nrow(O.i)
    
    ### PS predict
    p.A.i = p.A.list[[i]]
    
    ### IPCW
    IPCW.i = IPCW.list[[i]]
    
    ### Smooth Stochastic approximation
    all_xb = as.matrix(cbind(1,X.i)) %*% beta  
    all_prob = (A.i)*plogis(all_xb) + (1-A.i)*(1-plogis(all_xb))
    pd_prob = rep(prod(all_prob), N.i)
    
    ### PS predict
    ps_prob = (A.i)*p.A.i + (1-A.i)*(1-p.A.i)
    cluster_prob = rep(prod(ps_prob), N.i)
    
    dev_wts = rbind(
      dev_wts,  
      mean( (pd_prob/cluster_prob) * IPCW.i * R(Y.i) ) 
    )
  }
  
  mean_val = sum(dev_wts)/m
  
  return(-mean_val)
}

learning.ClusIPW <- function(data){
  
  ## Number of clusters
  m = length(unique(data$id))
  
  ## Propensity score estimation
  A.fit <- A.train(data, X.A.names)
  
  ## C model estimation (Cox PH regression)
  G.fit <- G.train(data, X.C.names)
  
  ## Save PS and IPCW to save time, which are repeatedly used in obj function
  p.A.true.list = list()
  p.A.fit.list = list()
  IPCW.noCensoring.list = list()
  IPCW.fit.list = list()
  
  for(i in 1:m){
    
    O.i = data %>% filter(id == i)
    N.i = nrow(O.i)
    
    ### PS true
    p.A.true.list[[i]] = p.A(O.i, N.i)
    
    ### PS predict
    p.A.fit.list[[i]] = predict(A.fit, O.i, type = "response")
    
    ### Noncensoring case => IPCW = 1 and change Y = T
    IPCW.noCensoring.list[[i]] = rep(1,N.i)
    
    ### IPCW predict
    IPCW.fit.list[[i]] = IPCW(O.i, G.fit)
    
  }
  
  ## 1 is for intercept 
  X_all = cbind(1, data %>% select(all_of(X.names)))  
  
  ## dimension of beta
  d = ncol(X_all)
  
  beta.ClusIPW = list()
  
  ### (i) Known PS, No Censoring
  beta.ClusIPW[["KN"]] = optim(rep(0,d), 
                       obj_fun_IPW, 
                       data = data %>% mutate(Y = T),         ## To account for noCensoring
                       p.A.list = p.A.true.list,
                       IPCW.list = IPCW.noCensoring.list,
                       method="L-BFGS-B", 
                       lower = -1, upper = 1, 
                       control = list(maxit = 1))$par
  
  beta.ClusIPW[["KN"]] = normalize_vector(beta.ClusIPW[["KN"]])
  
  
  ### (ii) Unknown PS (estimated), No Censoring
  beta.ClusIPW[["UN"]] = optim(rep(0,d), 
                       obj_fun_IPW, 
                       data = data %>% mutate(Y = T),         ## To account for noCensoring
                       p.A.list = p.A.fit.list,               ## Fitted PS
                       IPCW.list = IPCW.noCensoring.list,
                       method="L-BFGS-B", 
                       lower = -1, upper = 1, 
                       control = list(maxit = 1))$par
  
  beta.ClusIPW[["UN"]] = normalize_vector(beta.ClusIPW[["UN"]])
  
  ### (iii) Unknown PS (estimated), Censoring
  beta.ClusIPW[["UC"]] = optim(rep(0,d), 
                               obj_fun_IPW, 
                               data = data,         
                               p.A.list = p.A.fit.list,               ## Fitted PS
                               IPCW.list = IPCW.fit.list,     ## Estimated IPCW
                               method="L-BFGS-B", 
                               lower = -1, upper = 1, 
                               control = list(maxit = 1))$par
  
  beta.ClusIPW[["UC"]] = normalize_vector(beta.ClusIPW[["UC"]])

  return(beta.ClusIPW)
  
}

### Example ###
# learning.ClusIPW(data)


###----------- Method 4: Oracle Indirect Learning ---------------###

# Optimization with mu_RT function involves smooth stochastic approximation
# at the computation of objective_fun 

obj_fun_Indirect = function(beta, data) { 
  
  m = length(unique(data$id))  # number of clusters
  
  mean_outcome = sapply(1:m, function(i){
    
    O.i = data %>% filter(id == i)
    N.i = nrow(O.i)
    X.i = cbind(1, O.i %>% select(all_of(X.names)))
    
    ## Treatment assignment based on the linear rule
    A = as.numeric((as.matrix(X.i) %*% beta >= 0))
    
    value.i = mean(sapply(1:N.i, function(j) mu_RT(j,A,X.i,N.i)))
    
    return(value.i) 
    
  })
  
  value = mean(mean_outcome)
  
  return(-value)
}

learning.Oracle <- function(data){
  
  beta.Oracle = optim(beta.Oracle.init, 
                      obj_fun_Indirect, 
                      data = data,
                      method="L-BFGS-B")$par
  
  return(beta.Oracle)
  
}

### Example ###
# learning.Oracle(data)




######---------------------------------------------------------------------#####
######                  True Policy Value computation                      #####
######---------------------------------------------------------------------#####

###------- Compute value V(\pi) given for each \pi, defined by \beta -------###

## eval_performance function generatesM copies of clusters and compute value functions 
## under policy \pi. M needs to be large to well-approximate the value
## (each sampled cluster serves as a Monte Carlo sample)

eval_performance <- function(beta, M = 10000){
  
  results = sapply(1:M, function(u){
    
    N = N.dist(1)
    
    ## Step1. Covariate generation
    X = as.data.frame(cbind(1, X.gen(N)))
    
    ## Step2. Treatment assignment based on the linear rule
    A = as.numeric((as.matrix(X) %*% beta >= 0))
    mean_treated = mean(A)
    
    value = mean(sapply(1:N, function(j) mu_RT(j,A,X,N)))
    
    return(cbind(value = value, proportion = mean_treated )) 
  })
  return(round(apply(as.matrix(results),1,mean),3))
  
}

######---------------------------------------------------------------------#####
######           Parallel Computing to Run Multiple Simulations            #####
######---------------------------------------------------------------------#####

## Register Parallel backend
no_cores <- min(detectCores() - 1, 10)
registerDoParallel(cores = no_cores)

## Number of clusters in one simulation
generate_simul_result <- function(idx, m){
  
  print(glue("Simulation No: {idx} || m = {m}"))
  
  ### Simulation Data ###
  data = data.sim(m)
  
  ###----------- Method 1: Additive IPW ---------------###
  beta.addIPW = learning.addIPW(data)
  print(glue(". AddIPW done"))
  
  ###----------- Method 2: No Interference IPW ---------------###
  beta.NoIntIPW = learning.NoIntIPW(data)
  print(glue(". No Interference IPW done"))
  
  ###----------- Method 3: Cluster IPW ---------------###
  beta.ClusIPW = learning.ClusIPW(data)  
  print(glue(". Cluster IPW done"))
  
  ###----------- Method 4: Oracle Indirect ---------------###
  beta.Oracle = learning.Oracle(data)
  print(glue(". Oracle Indirect done"))
  Value.Oracle   = eval_performance(beta = beta.Oracle  ) # Oracle
  
  Result.list = list()
  
  for(Scenario in c("KN", "UN", "UC")){
    
    ###----------- Evaluate learned policy by V(\pi) ---------------###
    Value.NoIntIPW = eval_performance(beta = beta.NoIntIPW[[Scenario]]) # IPW without interference
    Value.ClusIPW  = eval_performance(beta = beta.ClusIPW [[Scenario]]) # IPW with interference
    Value.addIPW   = eval_performance(beta = beta.addIPW  [[Scenario]]) # AddIPW
    
    Result = as.data.frame(
      rbind(
        c(Value.NoIntIPW, normalize_vector(beta.NoIntIPW[[Scenario]])),
        c(Value.ClusIPW,  normalize_vector(beta.ClusIPW [[Scenario]]) ),
        c(Value.addIPW,   normalize_vector(beta.addIPW  [[Scenario]])  ),
        c(Value.Oracle,   normalize_vector(beta.Oracle)  ))
    )
    
    rownames(Result) = c("NoIntIPW", "ClusIPW", "addIPW", "Oracle")
    colnames(Result) = c("Value", "Tx_Prop", paste0("beta", 0:(length(beta.Oracle)-1)))
    
    write.csv(Result, glue("data/Output_{Scenario}_m{m}_{SLURMid}_{idx}.csv"))
    
    Result.list[[Scenario]] = Result
    
  }
  
  return(Result.list)
  
}

# generate_simul_result(idx=1, m=50)

# Run the parallel loop
m = 50
results <- foreach(idx = 1:Simul_Num, .combine = 'rbind', .packages = 'dplyr') %dopar% {
  generate_simul_result(idx, m)
}

m = 100
results <- foreach(idx = 1:Simul_Num, .combine = 'rbind', .packages = 'dplyr') %dopar% {
  generate_simul_result(idx, m)
}

m = 200
results <- foreach(idx = 1:Simul_Num, .combine = 'rbind', .packages = 'dplyr') %dopar% {
  generate_simul_result(idx, m)
}

m = 400
results <- foreach(idx = 1:Simul_Num, .combine = 'rbind', .packages = 'dplyr') %dopar% {
  generate_simul_result(idx, m)
}

# m = 800
# results <- foreach(idx = 1:Simul_Num, .combine = 'rbind', .packages = 'dplyr') %dopar% {
#   generate_simul_result(idx, m)
# }


# Stop the cluster
stopImplicitCluster()


