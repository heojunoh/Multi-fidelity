cor.sep <- function(X, x=NULL, theta, nu, derivative=0){
  d <- NCOL(X)
  n <- NROW(X)
  nu <- rep(nu, d)
  if(is.null(x)){
    K <- matrix(1, n, n)
    for(i in 1:d){
      R <- sqrt(distance(X[,i]/theta[i]))
      K <- K * matern.kernel(R, nu=nu[i], derivative=derivative)
    }
  }else{
    n.new <- NROW(x)
    K <- matrix(1, n, n.new)
    for(i in 1:d){
      R <- sqrt(distance(X[,i]/theta[i], x[,i]/theta[i]))
      K <- K * matern.kernel(R, nu=nu[i], derivative=derivative)
    }
  }
  return(K)
}


addcor.sep <- function(X, x=NULL, theta, input, nu, derivative=0){
  n <- NROW(X)
  if(is.null(x)){
    K <- matrix(0, n, n)
    for(i in 1:length(input)){
      K <- K + theta[length(c(unique(unlist(input, use.names=FALSE))))+order(unique(sapply(input, length)))[length(input[[i]])==unique(sapply(input, length))]]*cor.sep(as.matrix(X[,input[[i]]]), theta=theta[input[[i]]], nu=nu[[i]])
    }
  }else{
    n.new <- NROW(x)
    K <- matrix(0, n, n.new)
    for(i in 1:length(input)){
      K <- K + theta[length(c(unique(unlist(input, use.names=FALSE))))+order(unique(sapply(input, length)))[length(input[[i]])==unique(sapply(input, length))]]*cor.sep(as.matrix(X[,input[[i]]]), as.matrix(x[,input[[i]]]), theta=theta[input[[i]]], nu=nu[[i]])
    }
  }
  return(K)
}

