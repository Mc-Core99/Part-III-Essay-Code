#in this experiment we see what happens
#when we put Theorem 3.1 into practice



library(tidyverse)
#used to produce plots
library(Matrix)
#used for some matrix manipulation
library(BB)
#used to perform the optimisation

set.seed(4)

#generating random hyperplanes
random.Hyperplanes <- function(number, X) {
  hyperplanes <<- array(NA, c(length(X[, 1]), length(X[, 1]), number))
  for (j in (1:number)) {
    #forming the normal vectors, creating the
    #corresponding hyperplane arrangement
    g <- as.matrix(rnorm(length(X[1, ])))
    D <- matrix(0, length(X[, 1]), length(X[, 1]))
    for (i in 1:length(X[, 1])) {
      D[i, i] <- sum(((X%*%g)[i] >=0))
    }
    #adding the hyperplane arrangement to 
    #the set of all hyperplane arrangements
    hyperplanes[,,j] <<- D
  }
}

#creating the objective function for the relaxed problem

relaxed.Objective <- function(number, hyperplanes, X, y, lambda) {
  relaxed.objective.function <<- function(x) {
    #calculating sumDX(u-v)
    hyp.sum <- 0
    d <- length(hyperplanes[1, 1, ])
    l <- length(X[1, ])
    for (i in 1:number) {
      hyp.sum <- hyp.sum + hyperplanes[,, i]%*%X%*%(x[(1+(i-1)*l):(i*l)]-x[(d*l+1+(i-1)*l):(d*l + i*l)])
    }
    #calculating 1/2||sumDX(u-v) - y||^2
    first.term <-  0.5*sum((hyp.sum-y)^2)
    #calculating sum ||u_i||_2 + ||v_i||_2
    l2.norm1 <- 0
    for (i in 1:number) {
      l2.norm1 <- l2.norm1 + sqrt(sum(x[(1+(i-1)*l):(i*l)]^2)) + sqrt(sum(x[(d*l+1+(i-1)*l):(d*l + i*l)])^2)
    }
    #calculating the second term
    second.term <- lambda*l2.norm1
    #outputting the objective value at
    #the values u, v
    return(first.term + second.term)
  }
}


#creating the gradient of the objective function
#we limit any issues of dividing 0 by 0 using pmax
#this helps speed up optimisation of the relaxed objective


relaxed.gradient <- function(number, hyperplanes, X, y, lambda) {
  relaxed.objective.derivative <<- function(x) {
    l <- length(X[1, ])
    g <- rep(0, 2*number*l)
    h <- number
    hyp.sum <- 0
    d <- length(hyperplanes[1, 1, ])
    for (i in 1:number) {
      hyp.sum <- hyp.sum + hyperplanes[,, i]%*%X%*%(x[(1+(i-1)*l):(i*l)]-x[(d*l+1+(i-1)*l):(d*l + i*l)])
    }
    #calculating it for the u's
    for (i in 1:h) {
      for (j in 1:l) {
        ent <- (i-1)*l + j
        g[ent] <- lambda*x[ent]/pmax(norm(as.matrix(x[((i-1)*l + 1):(i*l)]), type = "2"), 0.0001) + (t(hyperplanes[, , i]%*%X)%*%(hyp.sum- y))[j]
      }
    }
    #calculating it for the v's
    for (i in 1:h) {
      for (j in 1:l) {
        ent <- h*l + (i-1)*l + j
        g[ent] <- lambda*x[ent]/pmax(norm(as.matrix(x[(h*l +(i-1)*l + 1): (h*l + i*l)]), type = "2"), 0.0001) - (t(hyperplanes[, , i]%*%X)%*%(hyp.sum - y))[j]
      }
    }
    return(g)
  } 
}


#creating the prjection function parameters
#this is used for the constaints in the minimisation

constraints <- function(hyperplanes, X) {
  l <- length(X[1, ])
  d <- length(hyperplanes[1, 1, ])
  k <- length(hyperplanes[, 1, 1])
  k2 <- length(hyperplanes[1, , 1])
  A <- NULL
  #creating the matrix defining the constraint sets K_D_i
  A <- (2*hyperplanes[, , 1] - diag(1, k, k2))%*%X
  for (i in 2:d) {
    A <- bdiag(A, (2*hyperplanes[, , i] - diag(1, k, k2))%*%X)
  }
  for (i in 1:d) {
    A <- bdiag(A, (2*hyperplanes[, , i] - diag(1, k, k2))%*%X)
  }
  AMat <<- as.matrix(A)
}

#perform the minimisation for the relaxed objective
#make trace FALSE to remove output

relaxed.minimum <- function(AMat, relaxed.objective.function, b1, relaxed.objective.derivative) {
  relaxed.solution <<- spg(par = rep(0,2*length(hyperplanes[1, 1, ])*length(X[1, ])), 
                           fn = relaxed.objective.function,
                           gr = relaxed.objective.derivative,
                           project = "projectLinear",
                           projectArgs = list(A = as.matrix(AMat), b=b1, meq = 0), 
                           control = list(ftol = 1e-3, trace = TRUE, checkGrad = FALSE))
}


#making the true objective function of the
#neural network

true.objective.function <- function(width, X, y, lambda) {
  true.objective <<- function(x) {
    #start by calculating the first term 
    #1/2 ||sum(Xu_j)_+a_j - y||_2^2
    term1 <- 0
    p <- length(X[1, ])
    w <- width
    for (i in 1:w) {
      term1 <- term1 + pmax(X%*%x[(1+(i-1)*p):(i*p)], 0)%*%x[i+w*p]
    }
    first.term <-  0.5*sum((term1-y)^2)
    
    #calculating the second term
    term2 <- 0
    for (i in 1:w) {
      term2 <- term2 + sum(x[(1+(i-1)*p):(i*p)]^2) + sum(x[(i+w*p)]^2) 
    }
    second.term <- lambda*0.5*term2
    return(first.term + second.term)
  }
}

#minimising the true neural network objective

true.minimum <- function(true.objective, epochs, learning.rate) {
  x <- rnorm(p*width + width, mean = 0, sd = 0.1)
  for (e in 1:epochs) {
    grad <- c()
    noise <- rnorm(length(x), mean = 0, sd = 0.1)
    for (i in 1:length(x)) {
      h <- 0.000001
      xderv <- x
      xderv[i] <- x[i] + h
      xderv2 <- x
      xderv2[i] <- x[i] - h
      grad[i] <- (true.objective(xderv) - true.objective(xderv2))/(2*h)
    }
    grad <- grad + noise
    x <- x - learning.rate*grad
  }
  true.sol <<- true.objective(x)
}


#setting up parameters

c <- 5
lambda <- 0.3
ratio.temp <- matrix(data = NA, nrow= 3, ncol = 9)

for (q in 1:3) {
  ratio1 <- c()
  k <- 1
  
  #iterating over all different n is slow
  #for a quicker plot, make the end point 80
  #you will also need to change to 80 for the
  #computation of baseline and the code that produces
  #the plot
  
  for (j in seq(30, 70, by =5)) {
    print(j)
    #creating a random design matrix and output y
    width <- ceiling(0.2*j)
    number <- ceiling(width*0.5)
    n <- j
    p <- ceiling(j/c)
    X <- matrix(rnorm(n*p), n, p)
    y <- matrix(rnorm(n, mean = 1/sqrt(n), sd = 1/sqrt(n)), n, 1)
    
    
    #computing the minimum of the convex relaxation
    random.Hyperplanes(number, X)
    relaxed.Objective(number, hyperplanes, X, y, lambda)
    relaxed.gradient(number, hyperplanes, X, y, lambda)
    constraints(hyperplanes, X)
    b1 = as.numeric(rep(0, length(AMat[, 1])))
    relaxed.minimum(AMat, relaxed.objective.function, b1, relaxed.objective.derivative)
    re.sol <- relaxed.solution$value
    
    #calculating the minimum of the true objective function
    #we repeat to try and avoid local minima
    #and find the closest approximation to the global minimum
    true.objective.function(width, X, y, lambda)
    minimum <- c()
    for (w in 1:3) {
      true.minimum(true.objective, epochs = 200, learning.rate = 0.05)
      minimum[w] <- true.sol
    }
    true.sol <- min(minimum)
    ratio1[k] <- re.sol/true.sol
    k <-  k+1
  }
  ratio.temp[q, ] <- ratio1
}
ratio.temp <- pmin(ratio.temp, mean(ratio.temp) + sd(ratio.temp))
std.dev <- apply(ratio.temp, 2, sd)
ratio <- colMeans(ratio.temp)
ymin <- ratio - std.dev
ymax <- ratio + std.dev
#transforming the ratio
ratio <- as.data.frame(ratio)


#theoretical bound when assumption are valid
#but we take the unknown constant to be 1
#so plotting the shape of the upper boundary
baseline <- sqrt(log(2*seq(25, 75, by =0.1)))
baseline <- as.data.frame(baseline)

#plotting the results

Legend <- c("Ratio" = "blue", "Shape-of-Upper-Bound" = "red", "Errors" = "black")

df <- data.frame(ratio,
                 ymax,
                 ymin)

names(df) <- c('ratio', 'ymax', 'ymin')

#to produce a smoothed plot, change line 232 to 
#geom_smooth()


ggplot(df, aes(x = seq(30, 70, by =5), y = ratio, color = "Ratio")) +
  geom_point() +
  geom_smooth(aes(color = "Ratio"), se = FALSE) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = "Errors")) +
  geom_line(data = baseline, mapping = aes(x = seq(25, 75, by =0.1), y = baseline, color = "Shape-of-Upper-Bound")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "purple") +
  labs(
    x = "n",
    y = "Ratio",
    title = "Comparison of Experimental & Theoretical Results For Normally Generated Data"
  ) +
  scale_color_manual(values = Legend) +
  coord_cartesian(ylim = c(0.9, 2.3)) +
  scale_y_continuous(breaks = seq(0.75, 2.5, 0.25))


