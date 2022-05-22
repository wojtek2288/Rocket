# Infers time series length and number of channels / dimensions 
# from input DataFrame, and generates random kernels

generate_kernels <- function(n_timepoints, num_kernels, n_columns, seed)
{
  set.seed(seed)
  candidate_lengths <- c(7, 9, 11)
  lengths <- sample(candidate_lengths, num_kernels, replace = TRUE)
  
  num_channel_indices = integer(num_kernels)
  
  for (i in 1:num_kernels)
  {
    limit <- min(n_columns, lengths[i])
    num_channel_indices[i] <- as.integer(2 ** runif(1, 0, log2(limit + 1)))
  }
  
  channel_indices = integer(sum(num_channel_indices))
  
  weights <- integer(as.integer(lengths %*% num_channel_indices))
  
  biases <- integer(num_kernels)
  dilations <- integer(num_kernels)
  paddings <- integer(num_kernels)
  
  a1 <- 1
  a2 <- 1
  
  for (i in 1:num_kernels)
  {
    normal_weights <- rnorm(num_channel_indices[i] * lengths[i], 0, 1)
    
    b1 <- a1 + (num_channel_indices[i] * lengths[i])
    b2 <- a2 + num_channel_indices[i]
    
    a3 <- 1
    
    for (j in 1:num_channel_indices[i])
    {
      b3 <- a3 + lengths[i]
      normal_weights[a3:(b3-1)] <- normal_weights[a3:(b3-1)] - mean(normal_weights[a3:(b3-1)])
      
      a3 <- b3
    }
    
    weights[a1:(b1-1)] <- normal_weights
    
    channel_indices[a2:(b2-1)] <- sample(seq(0, n_columns - 1), num_channel_indices[i])
    
    biases[i] <- runif(1, -1, 1)
    
    logVal <- log2((n_timepoints - 1)/(lengths[i] - 1))
    
    if (logVal > 0)
    {
      dilation <- 2 ** runif(1, 0, logVal)
    }
    else
    {
      dilation <- 2 ** runif(1, logVal, 0)
    }
    
    dilation <- as.integer(dilation)
    dilations[i] <- dilation
    
    if (sample(c(0,1), 1) == 0)
    {
      paddings[i] <- as.integer(((lengths[i] - 1) * dilation) / 2)
    }
    else
    {
      paddings[i] <- 0
    }
    
    a1 <- b1
    a2 <- b2
  }
  
  list("weights" = weights,
       "lengths" = lengths,
       "biases" = biases,
       "dilations" = dilations,
       "paddings" = paddings,
       "num_channel_indices" = num_channel_indices,
       "channel_indices" = channel_indices)
}

apply_kernel_univariate <- function(
    X, weights, length, bias, dilation, padding)
{
  n_timepoints <- length(X)
  
  output_length <- (n_timepoints + (2 * padding)) - ((length - 1) * dilation)
  
  ppv <- 0
  max1 <- -Inf
  
  end <- (n_timepoints + padding) - ((length - 1) * dilation)
  
  for (i in -padding:end)
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
  
  c(as.numeric(ppv)/as.numeric(output_length), as.numeric(max1))
}

apply_kernel_multivariate <- function(
  X, weights, length, bias, dilation,
  padding, num_channel_indices, channel_indices)
{
  dimensions <- dim(X)
  
  n_columns <- dimensions[1]
  n_timepoints <- dimensions[2]
  
  output_length <- (n_timepoints + (2 * padding)) - ((length - 1) * dilation)
  
  ppv <- 0
  max1 <- -Inf
  
  end <- (n_timepoints + padding) - ((length - 1) * dilation)
  
  for (i in -padding:end-1)
  {
    sum1 <- bias
    
    index <- i
    
    for (j in 1:length)
    {
      if (index > -1 & index < n_timepoints)
      {
        for (k in 1:num_channel_indices)
        {
          sum1 <- sum1 + weights[k, j] * X[channel_indices[k] + 1, index + 1]
        }
      }
      
      index <- index + dilation
    }
    
    if (sum1 > max1)
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

apply_kernels <- function(
  X,
  weights,
  lengths,
  biases,
  dilations,
  paddings,
  num_channel_indices,
  channel_indices)
{
  dimensions <- dim(X)
  
  n_instances <- dimensions[1]
  n_columns <- dimensions[2]
  
  num_kernels = length(lengths)
  
  X1 <- array(c(numeric(n_instances), numeric(num_kernels * 2)), 
              dim = c(n_instances, num_kernels * 2))
  
  for (i in 1:n_instances)
  {
    cat(i, " of ", n_instances, "\n")
    
    a1 <- 1
    a2 <- 1
    a3 <- 1
    
    for (j in 1:num_kernels)
    {
      b1 = a1 + num_channel_indices[j] * lengths[j]
      b2 = a2 + num_channel_indices[j]
      b3 = a3 + 2
      
      if (num_channel_indices[j] == 1)
      {
        vec <- apply_kernel_univariate(
          X[i, channel_indices[a2] + 1, ],
          weights[a1:(b1-1)],
          lengths[j],
          biases[j],
          dilations[j],
          paddings[j])
        
        X1[c(i), c(a3:(b3-1))] <- vec
      }
      else
      {
        weights1 <- array_reshape(weights[a1:(b1 - 1)],
                                  c(num_channel_indices[j], lengths[j]))
        
        vec <- apply_kernel_multivariate(
          X[i, , ],
          weights1,
          lengths[j],
          biases[j],
          dilations[j],
          paddings[j],
          num_channel_indices[j],
          channel_indices[a2:(b2-1)])
        
        X1[i, a3:(b3-1)] <- vec
      }
      
      a1 = b1
      a2 = b2
      a3 = b3
    }
  }
  
  X1
}

transformMatrix <- function(data)
{
  data3d <- array(numeric(), c(nrow(data), 2, ncol(data)))
  
    for (i in 1:nrow(data))
    {
      for (j in 1:ncol(data))
      {
        data3d[i, 1, j] <- data[i, j]
        data3d[i, 2, j] <- data[i, j]
      }
    }
  
  data3d
}

rocket <- function ()
{
  data <- read.arff("C:\\Users\\User\\Desktop\\Rocket\\ArrowHead_TRAIN.arff")
  data <- data[, 1:(ncol(data) - 1)]
  data3d <- transformMatrix(data)
  
  kernels <- generate_kernels(ncol(data), 1000, 2, 124)
  res <- apply_kernels(data3d,kernels[[1]], kernels[[2]], kernels[[3]],
                       kernels[[4]], kernels[[5]], kernels[[6]], kernels[[7]])
  
  write.table(res, file="C:\\Users\\User\\Desktop\\Rocket\\result_train.txt",
              row.names=F, sep=",")
  
  data <- read.arff("C:\\Users\\User\\Desktop\\Rocket\\ArrowHead_TEST.arff")
  data <- data[, 1:(ncol(data) - 1)]
  data3d <- transformMatrix(data)
  
  res <- apply_kernels(data3d,kernels[[1]], kernels[[2]], kernels[[3]],
                       kernels[[4]], kernels[[5]], kernels[[6]], kernels[[7]])
  
  write.table(res, file="C:\\Users\\User\\Desktop\\Rocket\\result_test.txt",
              row.names=F, sep=",")
}
