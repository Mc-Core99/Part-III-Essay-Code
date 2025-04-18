# Here we test Theorem 4.1 in the case of
# spike free matrices

#to make the code quicker change 
#n.restart to a lower value for each call to
#gosolnp

library(tidyverse)
#used to produce plots
library(Matrix)
#used for some matrix manipulation
library(Rsolnp)
#used to perform the optimisation
set.seed(1)


#first we define the spike-free objective function


spike.objective <- function(X, p) {
  spike.objective.function <<- function(x) {
    norm(x[1:p], type = "2") + norm(x[(1+p):(2*p)], type = "2")
  }
}

#now we define the spike-free inequality constraints

spike.ineq.constraints <- function(X, p) {
  spike.ineq <<- function(x) {
    c(X%*%x[1:p], 
    X%*%x[(p+1):(2*p)])
  }
}

#now we define the equality constraints

spike.equal.constraints <- function(X, p, y) {
  spike.equal <<- function(x) {
    X%*%(x[1:p] - x[(1+p):(2*p)])
  }
}



#optimising the spike-free objective
#make trace 1 to see outputs of the 
#algorithm

spike.optimal <- function(X, y, n, p) {
  spike.minimum <<- gosolnp(
                          fun = spike.objective.function,
                          eqfun = spike.equal,
                          eqB = y,
                          ineqfun = spike.ineq,
                          ineqLB = rep(0, 2*n),
                          ineqUB = rep(1000000000, 2*n),
                          LB = rep(-1000000, 2*p),
                          UB = rep(1000000, 2*p),
                          control = list(tol = 1e-8, delta = 1e-7, trace = 0),
                          distr = rep(1, 2*p),
                          n.restarts = 10
                          )
}



#creating the true objective function
#again we put the equality constraint in the 
#objective function


true.objective <- function(X, width, p) {
  true.objective.function <<- function(x) {
    norm(x[1:width], type = "2")^2 + norm(x[(width + 1): (p*width + width)], type = "2")^2
  }
}


#the equality constraint, with no bias term

true.equal.constraint <- function(X, y, n, p, width) {
  true.equal <<- function(x) {
    pmax(X%*%matrix(x[(width + 1):(p*width + width)], p, width)+ as.matrix(rep(1, n))%*%t(as.matrix(x[(p*width+width+1): ((p+2)*width)])), 0)%*%x[1:width]
  }
}


#optimising the true objective
#make trace 1 to see outputs of the 
#algorithm


true.optimal <- function(X, y, n, p, width){
  true.minimum <<- gosolnp( fun = true.objective.function, 
                            eqfun = true.equal,
                            eqB = y,
                            LB = rep(-1000000, (p+2)*width),
                            UB = rep(1000000, (p+2)*width),
                            control = list(tol = 1e-8, delta = 1e-7, trace = 0),
                            distr = rep(1, (p+2)*width), 
                            n.restarts = 3)
}

ratio <- c()
k <- 1
for (i in 2:6) {
  #forming the parameters
  X <-  diag(rep(1, i), i, i)
  n <-  i
  p <-  i
  width <- 4*i
  y <- 10*rnorm(i)
  
  
  #spike-free optimisation
  spike.objective(X, p)
  spike.ineq.constraints(X, p)
  spike.equal.constraints(X, p, y)
  spike.optimal(X,y, n, p)
  spike.min <- last(spike.minimum$value)
  print(spike.min)
  
  #true optimisation
  true.objective(X, width, p)
  true.equal.constraint(X, y, n, p, width)
  true.optimal(X, y, n, p, width)
  true.min <- last(true.minimum$value)
  print(true.min)
  ratio[k] <-  spike.min/true.min
  k <- k+1
}

print(round(ratio, digits = 1))


