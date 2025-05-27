library(doParallel)
library(foreach)
library(tscount)

# Setup parallel backend
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
clusterExport(cl, "tsglm")
# Outer and inner loop ranges
p <- 1:7
q <- 1:7

# Flatten into a grid of combinations
grid <- expand.grid(p = p, q = q)

# Parallel double loop
results <- foreach(k = 1:nrow(grid), .combine = rbind) %dopar% {
  p <- grid$p[k]
  q <- grid$q[k]
  
  # Your computation
  model <- tsglm(d_1_data[,2], model = list(past_obs = 1:p, past_mean = 1:q), 
                 xreg = d_1_data[,3:6], link='log', distr = "nbinom") 
  
  aic <- summary(model)$AIC
  
  c(p = p, q = q, aic = aic)
}

stopCluster(cl)