generate_kernels <- function(input_length, num_kernels)
{
  candidate_lengths <- c(7L, 9L, 11L)
  lengths <- sample(candidate_lengths, num_kernels, replace = TRUE)
  
  weights <- numeric(sum(lengths))
  biases <- numeric(num_kernels)
  dilations <- integer(num_kernels)
  paddings <- integer(num_kernels)
  
  a1 <- 1
  
  for (i in 1:num_kernels)
  {
    length <- lengths[i]
    
    normal_weights <- rnorm(length, 0, 1)

    b1 <- a1 + length
    weights[a1:(b1-1)] <- normal_weights - mean(normal_weights)
    
    biases[i] <- runif(1, -1, 1)
    
    logVal <- log2((input_length - 1)/(length - 1))

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
    
    
    if (sample(c(0,1), 1) == 1)
    {
      paddings[i] <- as.integer(((length - 1) * dilation) / 2)
    }
    else
    {
      paddings[i] <- 0
    }
    
    a1 <- b1
  }
  
  list("weights" = weights,
       "lengths" = lengths,
       "biases" = biases,
       "dilations" = dilations,
       "paddings" = paddings)
}

apply_kernel <- function(
    X, weights, length, bias, dilation, padding)
{
  input_length <- length(X)
  
  output_length <- (input_length + (2 * padding)) - ((length - 1) * dilation)
  
  ppv <- 0
  max <- -Inf
  
  end <- (input_length + padding) - ((length - 1) * dilation)
  
  for (i in -padding:end-1)
  {
    i <- -padding
    sum <- bias
    
    index <- i
    
    for (j in 1:length)
    {
      if (index > -1 & index < input_length)
      {
        sum <- sum + weights[j] * X[index+1]
      }
      
      index <- index + dilation
    }
    
    if(sum > max)
    {
      max <- sum
    }
    
    if (sum > 0)
    {
      ppv <- ppv + 1
    }
  }
  
  list("ppv" = ppv/output_length, "max" = max)
}

apply_kernel_mulivariate <- function(
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
          sum1 <- sum1 + weights[k, j] * X[channel_indices[k], index+1]
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