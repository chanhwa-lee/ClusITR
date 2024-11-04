###----------- Load Libraries ---------------###
library(dplyr)
library(tidyverse)
library(Rglpk) 
library(glue)
library(doParallel)

######---------------------------------------------------------------------#####
######                   Data generation function                          #####
######---------------------------------------------------------------------#####

## Step 0. Cluster size (N) distribution
N.dist = function(m) sample(x = c(5,10,15), size = m, replace = T)

## Step1. Covariate generation
X.gen <- function(N){
  X1 <- rnorm(N, 0, 1)
  X2 <- rnorm(N, 0, 1)
  X3 <- rnorm(N, 0, 1)
  X4 <- rnorm(N, 0, 1)
  X = data.frame(X1, X2, X3, X4)
  return(X)
}

## Step2. Treatment model (experiment setting)
p.A = 0.3
A.gen <- function(X, N){
  A <- rbinom(N, 1, p.A)
  return(A)
}

## Step 3. Outcome model (Generate Y from Normal)
mu_Y <- function(j, A, X, N){
  with(X,
       (0.2*X2[j] + 0.2*X3[j] + 
          (1*X1[j] + 0.5*X2[j] - 1*X3[j] - 0.5*X4[j])*A[j] + 
          1.5*1/(N-1)*( sum((X3+X4)*A) - (X3[j]+X4[j])*A[j] ))
       )
}

Y.gen <- function(A, X, N){
  sigma = 1
  Y <- sapply(1:N, function(j) rnorm(1, mu_Y(j, A, X, N), sigma))
}

### beta.Oracle from indirect learning - by maximizing \sum_j mu_Y(j,A,X,N)
beta.Oracle = c(0, 1, 0.5, 0.5, 1)

## Step 4. Cluster data generation
DGP <- function(N){
  
  ## Covariate generation
  X = X.gen(N)
  
  ## Treatment model
  A = A.gen(X, N)
  
  ## Outcome model
  Y = Y.gen(A, X, N)
  
  data = cbind(data.frame(Y = Y, A = A), X)
  
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
m = 10
data = data.sim(m)
summary(data)


######---------------------------------------------------------------------#####
######                    Policy Learning Methods                          #####
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

###----------- Method 1: Additive IPW ---------------###

learning.addIPW <- function(data){
  
  ## Number of clusters
  m = length(unique(data$id))
  
  ## Total observations
  nn_all = nrow(data) 
  
  ## 1 is for intercept 
  X_all = cbind(1, data %>% select(-c(id,Y,A)))  
  
  ## dimension of beta
  d = dim(X_all)[2] 
  
  ## Compute optimization weight at each cluster
  w.i = function(i){
    
    O.i = data %>% filter(id == i) 
    Y.i = O.i$Y
    A.i = O.i$A
    
    w.i = mean( Y.i )*(A.i/p.A - (1-A.i)/(1-p.A))
    
    return(w.i)
    
  }
  
  w_all = unlist(sapply(1:m, w.i))
  
  ## Objective function coefficients for (beta, p_{ij}) 
  # First d elements are coefficients of beta in the objective function, which are 0 (since they do not appear)
  # and remaining elements are coefficients of p_{ij}, which is \overline{Y}_i * (Aij / ej(1|Xi) - (1-Aij) / ej(0|Xi))
  f <- c(rep(0,d), w_all)  
  
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
  
  bounds= list(lower = list(ind = 1:length(lb),val =lb),
               upper = list(ind =  1:length(lb),val =ub))
  
  ## Variable type string: "C" for continuous (\beta_k) and "B" for binary (p_ij)
  ctype <- c(rep("C", d), rep("B", nn_all))
  
  ## Mixed Integer Programming
  sol.addIPW <- Rglpk_solve_LP(obj = f, mat = Aineq, dir = dir, rhs = bineq, bounds = bounds, types = ctype,max = T)
  
  ## First d elements are beta
  beta.addIPW = sol.addIPW$solution[c(1:d)]
  
  return(normalize_vector(beta.addIPW))
  
}

### Example ###
learning.addIPW(data)


###----------- Method 2: No Interference IPW ---------------###

learning.NoIntIPW <- function(data){
  
  ## Number of clusters
  m = length(unique(data$id))
  
  ## Total observations
  nn_all = nrow(data) 
  
  ## 1 is for intercept 
  X_all = cbind(1, data %>% select(-c(id,Y,A)))  
  
  ## dimension of beta
  d = dim(X_all)[2] 
  
  ## Compute optimization weight at each cluster
  w.NoInt.i = function(i){
    
    O.i = data %>% filter(id == i) 
    Y.i = O.i$Y
    A.i = O.i$A
    N.i = nrow(O.i)
    
    ## This part is the only difference from addIPW
    w.i = 1/N.i * Y.i * (A.i/p.A - (1-A.i)/(1-p.A))
    
    return(w.i)
    
  }
  
  w.NoInt_all = unlist(sapply(1:m, w.NoInt.i))
  
  ## Objective function coefficients for (beta, p_{ij}) 
  # First d elements are coefficients of beta in the objective function, which are 0 (since they do not appear)
  # and remaining elements are coefficients of p_{ij}, which is 1/N.i * Y_i * (Aij / ej(1|Xi) - (1-Aij) / ej(0|Xi))
  f <- c(rep(0,d), w.NoInt_all)  
  
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
  
  bounds = list(lower = list(ind = 1:length(lb),val =lb),
                upper = list(ind =  1:length(lb),val =ub))
  
  ## Variable type string: "C" for continuous (\beta_k) and "B" for binary (p_ij)
  ctype <- c(rep("C", d), rep("B", nn_all))
  
  ## Mixed Integer Programming
  sol.NoIntIPW <- Rglpk_solve_LP(obj = f, mat = Aineq, dir = dir, rhs = bineq, bounds = bounds, types = ctype,max = T)
  
  ## First d elements are beta
  beta.NoIntIPW = sol.NoIntIPW$solution[c(1:d)]
  
  return(normalize_vector(beta.NoIntIPW))
  
}

### Example ###
learning.NoIntIPW(data)

###----------- Method 3: Cluster IPW ---------------###

# Optimization with classical cluster-IPW value estimator involves smooth stochastic approximation
# at the computation of objective_fun 

Value_est_IPW = function(beta, data) { 
  
  m = length(unique(data$id))  # number of clusters
  
  dev_wts=NULL
  for(i in 1:m){
    
    O.i = data %>% filter(id == i) %>% select(-id)
    Y.i = O.i$Y
    A.i = O.i$A
    X.i = O.i %>% select(-c(Y, A))
    N.i = nrow(O.i)
    
    all_xb = as.matrix(cbind(1,X.i)) %*% beta  
    all_prob = (A.i)*plogis(all_xb) + (1-A.i)*(1-plogis(all_xb))
    pd_prob = rep(prod(all_prob), N.i)
    
    ### PS predict
    ps_prob = (A.i)*p.A + (1-A.i)*(1-p.A)
    cluster_prob = rep(prod(ps_prob), N.i)
    
    dev_wts = rbind(
      dev_wts,  
      mean( (pd_prob/cluster_prob) * Y.i ) 
    )
  }
  
  mean_val = sum(dev_wts)/m
  
  return(mean_val)
}

learning.ClusIPW <- function(data){
  
  ## Define objective function under classical IPW value estimator since optim can only minimize
  obj_fun_IPW = function(beta, data){
    
    obj_value = -Value_est_IPW(beta, data)
    
    return(obj_value)
  }
  
  ## dimension of beta
  d = ncol(data %>% select(-c(id,Y,A))) + 1 ## 1 is for intercept 
  
  beta.ClusIPW = optim(rep(0,d), 
                   obj_fun_IPW, 
                   data = data, 
                   method="L-BFGS-B", 
                   lower = -1, upper = 1, 
                   control = list(maxit = 1))$par
  
  return(normalize_vector(beta.ClusIPW))
  
}

### Example ###
learning.ClusIPW(data)


######---------------------------------------------------------------------#####
######                  True Policy Value computation                      #####
######---------------------------------------------------------------------#####

###------- Compute value V(\pi) given for each \pi, defined by \beta -------###

## eval_performance function generatesM copies of clusters and compute value functions 
## under policy \pi. M needs to be large to well-approximate the value
## (each sampled cluster serves as a Monte Carlo sample)

eval_performance <- function(beta, M = 1000){
  
  results = sapply(1:M, function(u){
    
    N = N.dist(1)
    
    ## Step1. Covariate generation
    X = as.data.frame(cbind(1, X.gen(N)))
    
    ## Step2. Treatment assignment based on the linear rule
    A = as.numeric((as.matrix(X) %*% beta >= 0))
    mean_treated = mean(A)
    
    value = mean(sapply(1:N, function(j) mu_Y(j,A,X,N)))
    
    return(cbind(value = value, proportion = mean_treated )) 
  })
  return(apply(as.matrix(results),1,mean))
  
}

######---------------------------------------------------------------------#####
######           Parallel Computing to Run Multiple Simulations            #####
######---------------------------------------------------------------------#####

## Register Parallel backend
no_cores <- min(detectCores() - 1, 8)
registerDoParallel(cores = no_cores)

## Number of simulations
Simul_Num = 50

## Number of clusters in one simulation
generate_simul_result <- function(idx, m){
  
  print(glue("Simulation No: {idx}"))
  
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
  
  
  ###----------- Evaluate learned policy by V(\pi) ---------------###
  Value.NoIntIPW = eval_performance(beta = beta.NoIntIPW) # IPW without interference
  Value.ClusIPW  = eval_performance(beta = beta.ClusIPW ) # IPW with interference
  Value.addIPW   = eval_performance(beta = beta.addIPW  ) # AddIPW
  Value.Oracle   = eval_performance(beta = beta.Oracle  ) # Oracle
  
  Result = as.data.frame(
    rbind(
      c(Value.NoIntIPW, normalize_vector(beta.NoIntIPW)),
      c(Value.ClusIPW,  normalize_vector(beta.ClusIPW) ),
      c(Value.addIPW,   normalize_vector(beta.addIPW)  ),
      c(Value.Oracle,   normalize_vector(beta.Oracle)  ))
  )
  
  rownames(Result) = c("NoIntIPW", "ClusIPW", "addIPW", "Oracle")
  colnames(Result) = c("Value", "Tx_Prop", paste0("beta", 0:4))
  
  write.csv(Result, glue("Output_Zhang/Output_m{m}_{idx}.csv"))
  
  return(Result)
  
}

generate_simul_result(idx = 3, m = 50)

# Run the parallel loop
m = 400
results <- foreach(idx = 1:Simul_Num, .combine = 'rbind', .packages = 'dplyr') %dopar% {
  generate_simul_result(idx, m)
}

# Stop the cluster
stopImplicitCluster()


######---------------------------------------------------------------------#####
######                  Visualization of Simulation Result                 #####
######---------------------------------------------------------------------#####

# Initialize empty data frames for both plots
policy_value_plot_data = data.frame()  # For policy value plot
beta_bias_plot_data = data.frame()  # For beta bias comparison plot

beta.Oracle <- normalize_vector(beta.Oracle)  # These are the true beta values

# Loop over different sample sizes
for(m in c(50, 100, 200, 400)){  # Include multiple m values as needed
  
  # Get a list of all CSV files in the folder for the given sample size m
  csv_files <- list.files(path = "./Output_Zhang", pattern = glue("Output_m{m}.*.csv"), full.names = TRUE)
  
  # Read all CSV files into a list of data frames
  data_list <- lapply(csv_files, read.csv)
  
  # Combine all data frames into one
  combined_data <- do.call(rbind, data_list)
  colnames(combined_data) = c("Method", "Value", "TxProp", paste0("beta", 0:4))
  combined_data$Method = factor(combined_data$Method, 
                                levels = c("NoIntIPW", "ClusIPW", "addIPW", "Oracle"), 
                                ordered = TRUE)
  
  combined_data$m = m  # Add the sample size column
  
  # Add to policy value plot data
  policy_value_plot_data = rbind(policy_value_plot_data, combined_data)
  
  ###----------- Beta Bias Comparison Plot ---------------###
  
  # Step 1: Reshape the data to long format
  combined_long <- combined_data %>%
    pivot_longer(cols = starts_with("beta"), names_to = "beta", values_to = "value")
  
  # Step 2: Add the true beta values from beta.Oracle based on the 'beta' column
  combined_long <- combined_long %>%
    mutate(beta_num = as.numeric(gsub("beta", "", beta)),  # Extract beta number (0 to 4)
           oracle_value = beta.Oracle[beta_num + 1])  # Add 1 since R is 1-based index
  
  # Step 3: Calculate the bias
  combined_long <- combined_long %>%
    mutate(bias = value - oracle_value) %>%
    filter(Method != "Oracle")  # Remove Oracle rows since we are using beta.Oracle for true values
  
  # Step 4: Add sample size `m` to the long data frame
  beta_bias_plot_data = rbind(beta_bias_plot_data, combined_long)
}

###----------- Policy Value Plot ---------------###

ggplot(policy_value_plot_data, 
       aes(x = as.factor(m), y = Value, fill = Method)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("NoIntIPW" = "darkgreen", 
                               "ClusIPW" = "blue", 
                               "addIPW" = "red", 
                               "Oracle" = "purple")) +
  labs(title = "Learned Policy Value by Method",
       x = "Sample Size", y = "Policy Value")

###----------- Beta Bias Comparison Plot ---------------###

# Draw the boxplot of bias for each beta value by Method and sample size
ggplot(beta_bias_plot_data, aes(x = beta, y = bias, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("NoIntIPW" = "darkgreen", 
                               "ClusIPW" = "blue", 
                               "addIPW" = "red")) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
  labs(title = "Bias of Beta Estimates by Method and Sample Size",
       x = "Beta",
       y = "Bias") +
  facet_wrap(~m, scales = "free") +  # Facet by sample size `m`
  theme_bw()
