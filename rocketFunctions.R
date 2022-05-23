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
      normal_weights[a3:(b3-1)] <-
        normal_weights[a3:(b3-1)] - mean(normal_weights[a3:(b3-1)])
      
      a3 <- b3
    }
    
    weights[a1:(b1-1)] <- normal_weights
    
    channel_indices[a2:(b2-1)] <-
      sample(seq(0, n_columns - 1), num_channel_indices[i])
    
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
    cat(i, "rows out of", n_instances, "\n")
    
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

readData <- function(dataName, count, trainOrTest)
{
  setwd(paste0(getwd(), "/data"))
  
  data <- read.arff(paste0(dataName, "Dimension1_", trainOrTest, ".arff"))
  classes <- data[, ncol(data)]
  data <- data[, 1:(ncol(data) - 1)]
  
  dataSets <- array(numeric(), c(nrow(data), count, ncol(data)))
  
  dataSets[, 1, ] <-
    array(data = c(unlist(data)), dim = c(nrow(data), ncol(data)))
  
  if (count > 1)
  {
    for (i in 2:count)
    {
      data <-
        read.arff(paste0(dataName, "Dimension", i, "_", trainOrTest, ".arff"))
      data <- data[, 1:(ncol(data) - 1)]
      
      dataSets[, i, ] <-
        array(data = c(unlist(data)), dim = c(nrow(data), ncol(data)))
    }
  }
  
  list("dataSets" = dataSets, "classes" = classes)
}

replaceValues <- function(dataSets, classes, useMean = TRUE)
{
  dims <- dim(dataSets)
  
  uniqueClasses <- unique(classes)
  groupedClasses <- table(classes) 
  
  for (i in 1:length(uniqueClasses))
  {
    val <- as.integer(groupedClasses[[uniqueClasses[i]]] / 3)
    classCount <- 0
    for (j in 1:length(classes))
    {
      if (classes[j] == uniqueClasses[i])
      {
        classCount <- classCount + 1
        if (classCount < val)
        {
          from <- sample(as.integer(dims[3] * 0.1):as.integer(dims[3] * 0.4), 1)
        }
        else if (classCount >= val & classCount < 2 * val)
        {
          from <- sample(as.integer(dims[3] * 0.4):as.integer(dims[3] * 0.7), 1)
        }
        else
        {
          from <- sample(as.integer(dims[3] * 0.7):dims[3], 1)
        }
        
        for (k in 1:dims[2])
        {
          if (useMean)
          {
            dataSets[j, k, from:dims[3]] <-
              rep(mean(dataSets[j, k, 1:(from - 1)]), dims[3] - from + 1)
          }
          else
          {
            dataSets[j, k, from:dims[3]] <- numeric(dims[3] - from + 1)
          }
        }
      }
    }
  }
  
  dataSets
}

rocket <- function (dataName, count, kernelCount, seed)
{
  library("foreign")
  library("reticulate")
  cat(paste0(getwd(), "/data"))
  setwd(paste0(getwd(), "/data"))
  
  trainDataSets <- readData(dataName, count, "TRAIN")
  testDataSets <- readData(dataName, count, "TEST")
  
  for (i in 1:10)
  {
    cat("Iteration", i, "of 10 \n")
    dir.create(paste0(getwd(), "/results"),
               showWarnings = FALSE)
    
    setwd(paste0(getwd(), "/results"))
    meanDataSets <- replaceValues(trainDataSets[[1]], trainDataSets[[2]], TRUE)
    
    write.table(trainDataSets[[2]], file=paste0(i, "_Y_TRAIN.txt"),
                row.names=F, sep=",")
    
    kernels <- generate_kernels(dim(trainDataSets[[1]])[3], kernelCount,
                                count, seed)
    
    cat("Train data set: \n")
    cat("Mean train data set: \n")
    
    res <- apply_kernels(meanDataSets, kernels[[1]], kernels[[2]], kernels[[3]],
                         kernels[[4]], kernels[[5]], kernels[[6]], kernels[[7]])
    
    write.table(res, file=paste0(i, "_Result_TRAIN.txt"),
                row.names=F, sep=",")
    
    meanDataSets <- replaceValues(testDataSets[[1]], testDataSets[[2]], TRUE)
    
    write.table(testDataSets[[2]], file=paste0(i, "_Y_TEST.txt"),
                row.names=F, sep=",")
    
    cat("Test data set: \n")
    cat("Mean test data set: \n")
    
    res <- apply_kernels(meanDataSets, kernels[[1]], kernels[[2]], kernels[[3]],
                         kernels[[4]], kernels[[5]], kernels[[6]], kernels[[7]])
    
    write.table(res, file=paste0(i, "_Result_TEST.txt"),
                row.names=F, sep=",") 
  }
}

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4)
{
  stop("4 arguments must be supplied", call.=FALSE)
}

rocket(args[1], strtoi(args[2]), strtoi(args[3]), strtoi(args[4]))





