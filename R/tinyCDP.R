#'
#' tinyCDP.sqr: the continuous DP matching using squared Euclidean distance
#' @export
#' @param x: a matrix of the input signal
#' @param y: a matrix of the signal to be detected
#' @return A list contains three elements:
#'  xsize: length (number of rows) of x
#'  ysize: length (number of rows) of y
#'  g: xsize*ysize matrix of the accumulated distance
#'  d: the distance matrix
#'  bp: the backpointer matrix

tinyCDP.sqr <- function(x,y) {
  #x.size <- nrow(x)
  #y.size <- nrow(y)
  #d <- array(Inf,dim=c(x.size,y.size))
  #for (i in 1:x.size) {
  #    for (j in 1:y.size) {
  #        d[i,j] <- sum((x[i,]-y[j,])^2)
  #    }
  #}
  d <- distMatrix(x,y)
  tinyCDP(d)
}

distMatrix <- function(X,Y) {
  Xd <- matrix(apply(X,1,function(x){sum(x^2)}),nrow=nrow(X),ncol=nrow(Y))
  Yd <- matrix(apply(Y,1,function(x){sum(x^2)}),nrow=nrow(X),ncol=nrow(Y),byrow=TRUE)
  Xd+Yd-2*X%*%t(Y)
}


#'
#' tinyCDP: the continuous DP matching
#' @export
#' @param d: the distance matrix (nrow(d) is the length of the database signal, ncol(d) is the length of the query signal)
#' @return A list contains three elements:
#'  xsize: length (number of rows) of x
#'  ysize: length (number of rows) of y
#'  g: xsize*ysize matrix of the accumulated distance
#'  d: the distance matrix
#'  bp: the backpointer matrix
#'
tinyCDP <- function(d) {
  x.size <- nrow(d)
  y.size <- ncol(d)
  g <- array(Inf,dim=c(x.size,y.size))
  bp <- array(0,dim=c(x.size,y.size))
  g[1,1] <- d[1,1]
  bp[1,1] <- 2
  for (i in 2:x.size) {
      g[i,1] <- d[i,1]
      bp[i,1] <- 2
      for (j in 2:y.size) {
          g3 <- Inf
          g1 <- g[i-2,j-1]+d[i-1,j]
          g2 <- g[i-1,j-1]
          if (j > 1) {
              g3 <- g[i-1,j-2]+d[i,j-1]
          }
          dists <- c(g1, g2, g3)
          m <- which.min(dists)
          g[i,j] <- dists[m]+d[i,j]
          bp[i,j] <- m
      }
  }
  return(list(xsize=x.size,ysize=y.size,
              g=g,d=d,bp=bp))
}

#'
#' backtrace: find the optimum path from the result of tinyCDP
#' @export
#' @param dp: result of DP maching calculated by tinyCDP()
#' @param pos: the detection position
#' @return the optimum path:  
backtrace <- function(dp,pos) {
  i <- pos
  j <- ncol(dp$g)
  opt <- t(as.matrix(c(i,j)))
  while (i > 1 && j > 1) {
    b <- dp$bp[i,j]
    if (b == 1) { 
      opt <- rbind(matrix(c(i-1,j),ncol=2,byrow=T),opt)
      i <- i-2
      j <- j-1
    }
    else if (b == 2) { 
      opt <- rbind(matrix(c(i-1,j-1),ncol=2),opt)
      i <- i-1
      j <- j-1 
    }
    else if (b == 3) { 
      opt <- rbind(matrix(c(i,j-1),ncol=2,byrow=T),opt)
      i <- i-1
      j <- j-2 
    }
    else { 
      print("Invalid backpointer")
      break 
    }
  }
  opt
}

#' detect.peaks: detect peaks from the vector
#' @param values a vector
#' @param thres the detection threshold
#' @param pwidth the inhibition width
#' @param type=c("min","max") detection type
#' @return vector of detection positions
detect.peaks <- function(values,thres,pwidth,type="min") {
  len <- length(values)
  p <- values[1:(len-1)]<values[2:len]
  q1 <- c(p,FALSE) # current value is smaller than the next value
  q2 <- c(FALSE,p) # current value is larger than the previous value
  ind <- 1:len
  if (type == "min") {
    q <- q1 & !q2
    ii <- ind[q]
    t <- values[ii] < thres
    ii <- ii[t]
    x <- data.frame(pos=ii,value=values[ii])
    od <- order(x$value)
  } else {
    q <- !q1 & q2
    ii <- ind[q]
    t <- values[ii] > thres
    ii <- ii[t]
    x <- data.frame(pos=ii,value=values[ii])
    od <- order(x$value,decreasing=TRUE)
  }
  n <- length(ii)
  for (i in 1:(n-1)) {
    if (is.na(x$value[od[i]])) {
      next
    }
    for (j in (i+1):n) {
      if (abs(ii[od[i]]-ii[od[j]]) < pwidth) {
        x$value[od[j]] <- NA
      }
    }
  }
  x <- x[!is.na(x$value),]
  x$pos
}
