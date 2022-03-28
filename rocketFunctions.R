# Infers time series length and number of channels / dimensions 
# from input DataFrame, and generates random kernels

generate_kernels <- function(n_timepoints, num_kernels, n_columns, seed)
{
  num_kernels = 10
  n_timepoints = 20
  n_columns = 5
  seed = 124
  
  set.seed(seed)
  candidate_lengths <- c(7, 9, 11)
  lengths <- sample(candidate_lengths, num_kernels, replace = TRUE)
  
  num_channel_indices = integer(num_kernels)
  
  for (i in 1:num_kernels)
  {
    limit <- min(n_columns, lengths[i])
    num_channel_indices[i] <- 2 ** runif(1, 0, log2(limit + 1))
  }
  
  channel_indices = integer(sum(num_channel_indices))
  
  weights <- integer(as.integer(lengths %*% num_channel_indices))
  
  biases <- integer(num_kernels)
  dilations <- integer(num_kernels)
  paddings <- integer(num_kernels)
  
  a1 <- 0
  a2 <- 0
  
  for (i in 1:num_kernels)
  {
    normal_weights <- rnorm(num_channel_indices[i] * lengths[i], 0, 1)
    
    b1 <- a1 + (num_channel_indices[i] * lengths[i])
    b2 <- a2 + num_channel_indices[i]
    
    a3 <- 0
    
    for (j in 1:num_channel_indices[i])
    {
      b3 <- a3 + lengths[i]
      normal_weights[a3+1:b3+1] <- normal_weights[a3+1:b3+1] - mean(normal_weights[a3+1:b3+1])
      
      a3 <- b3
    }
    
    weights[a1+1:b1+1] <- normal_weights
    
    channel_indices[a2+1:b2+1] <- sample(seq(0, n_columns - 1), num_channel_indices[i])
    
    biases[i] <- runif(1, -1, 1)
    
    dilation <- 2 ** runif(1, 0, log2((n_timepoints - 1)/(lengths[i] - 1)))
    dilation <- as.integer(dilation)
    
    dilations[i] <- dilation
    
    if (sample(c(0,1), 1) == 0)
    {
      padding <- as.integer(((lengths[i] - 1) * dilation) / 2)
    }
    else
    {
      padding <- 0
    }
    
    paddings[i] <- padding
    
    a1 <- b1
    a2 <- b2
  }
  
  
}

apply_kernel_univariate <- function(
    X, weights, length, bias, dilation, padding)
{
  # X <- c(1,2,3)
  # weights <- c(0.2, 0.1, 0.6)
  # length <- 3
  # bias <- 3
  # dilation <- 5
  # padding <- 3
  
  n_timepoints <- length(X)
  
  output_length <- (n_timepoints + (2 * padding)) - ((length - 1) * dilation)
  
  ppv <- 0
  max1 <- -Inf
  
  end <- (n_timepoints + padding) - ((length - 1) * dilation)
  
  for (i in -padding:end-1)
  {
    i <- -padding
    sum1 <- bias
    
    index <- i
    
    for (j in 1:length)
    {
      if (index > -1 & index < n_timepoints)
      {
        sum1 <- sum1 + weights[j] * X[index+1]
      }
      
      index <- index + dilation
    }
    
    if(sum1 > max1)
    {
      max1 <- sum1
    }
    
    if (sum1 > 0)
    {
      ppv <- ppv + 1
    }
  }
  
  c(ppv/output_length, max1)
}

fit <- function(X)
{
  dimensions <- dim(worms)
  kernels <- generate_kernels()
  
}