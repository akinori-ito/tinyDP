#'
#' tinyDP: A fast dynamic time warping
#'
#' @export
#' @param x A matrix
#' @param y A matrix
#' @param window The window width with which the DTW is calculated,
#' @param type If type==1, the DTW has not restriction on the length difference between x and y. If type==2, the length of the one matrix should be more than the half and less than the twice of the other matrix.
#' @return A list contains three elements:
#'  xsize: length (number of rows) of x
#'  ysize: length (number of rows) of y
#'  opt: n times 2 matrix, where a row contains the corresponding indices of x and y
#'  score: the matching distance

tinyDP <-
function(x,y,window=50,type=1) {
    if (type == 1) {
        tinyDP1(x,y,window)
    } else if (type == 2) {
        tinyDP2(x,y,window)
    } else {
        stop(paste("tinyDP: unknown type:",type))
    }
}

 
tinyDP1 <- function(x,y,window=50) {
    x.size <- dim(x)[1]
    y.size <- dim(y)[1]
    wlimit <- 2*window+1
    g <- array(Inf,dim=c(x.size,wlimit))
    d <- array(Inf,dim=c(x.size,wlimit))
    bp <- array(0,dim=c(x.size,wlimit))
    jcenter <- 1+floor((0:(x.size-1))/(x.size-1)*(y.size-1))
    j_ind <- function(i,j) { return(j-jcenter[i]+window+1) }
    g[1,j_ind(1,1)] <- d[1,j_ind(1,1)] <- sum((x[1,]-y[1,])^2)
    bp[1,j_ind(1,1)] <- 2
    for (i in 2:x.size) {
        jmin <- max(jcenter[i]-window,1)
        jmax <- min(jcenter[i]+window,y.size)
        for (j in jmin:jmax) {
            w <- j_ind(i,j)
            w1 <- j_ind(i-1,j)
            d[i,w] <- sum((x[i,]-y[j,])^2)
            g1 <- g2 <- g3 <- Inf
            if (w1 <= wlimit) {
                g1 <- g[i-1,w1]
            }
            if (w1-1 <= wlimit) {
                g2 <- g[i-1,w1-1]
            }
            if (w-1 <= wlimit) {
                g3 <- g[i,w-1]
            }
            dists <- c(g1, g2, g3)
            m <- which.min(dists)
#            message(sprintf("i=%d w=%d m=%d",i,w,m)) 
            g[i,w] <- dists[m]+d[i,w]
            bp[i,w] <- m
        }
    }
    i <- x.size
    j <- y.size
    opt <- t(as.matrix(c(i,j)))
    while (i > 1 && j > 1) {
                                        #print(c(i,j,bp[i,j_ind(i,j)]))
        #assert_that(i >= 1, i <= x.size, j_ind(i,j) >= 1, j_ind(i,j) <= wlimit)
        b <- bp[i,j_ind(i,j)]
        if (b == 1) { 
            opt <- rbind(matrix(c(i-1,j),ncol=2,byrow=T),opt)
            i <- i-1 
        }
        else if (b == 2) { 
            opt <- rbind(matrix(c(i-1,j-1),ncol=2),opt)
            i <- i-1
            j <- j-1 
        }
        else if (b == 3) { 
            opt <- rbind(matrix(c(i,j-1),ncol=2,byrow=T),opt)
            j <- j-1 
        }
        else { 
            print("Invalid backpointer")
            break 
        }
    }
    list(xsize=x.size,ysize=y.size,opt=opt,score=g[x.size,j_ind(x.size,y.size)])
}

tinyDP2 <- function(x,y,window=50) {
    x.size <- dim(x)[1]
    y.size <- dim(y)[1]
    wlimit <- 2*window+1
    g <- array(Inf,dim=c(x.size,wlimit))
    d <- array(Inf,dim=c(x.size,wlimit))
    bp <- array(0,dim=c(x.size,wlimit))
    jcenter <- 1+floor((0:(x.size-1))/(x.size-1)*(y.size-1))
    j_ind <- function(i,j) { return(j-jcenter[i]+window+1) }
    g[1,j_ind(1,1)] <- d[1,j_ind(1,1)] <- sum((x[1,]-y[1,])^2)
    bp[1,j_ind(1,1)] <- 2
    d[2,j_ind(2,2)] <- sum((x[2,]-y[2,])^2)
    g[2,j_ind(2,2)] <- g[1,j_ind(1,1)]+d[2,j_ind(2,2)]
    bp[2,j_ind(2,2)] <- 2
    for (i in 3:x.size) {
        jmin <- max(jcenter[i]-window,1)
        jmax <- min(jcenter[i]+window,y.size)
        for (j in jmin:jmax) {
            w <- j_ind(i,j)
            w1 <- j_ind(i-1,j)
            w2 <- j_ind(i-2,j)
            d[i,w] <- sum((x[i,]-y[j,])^2)
            g1 <- g2 <- g3 <- Inf
            if (w1 <= wlimit && w2-1 <= wlimit) {
                g1 <- g[i-2,w2-1]+d[i-1,w1]
            }
            if (w1-1 <= wlimit) {
                g2 <- g[i-1,w1-1]
            }
            if (w1-2 <= wlimit && w-1 <= wlimit) {
                g3 <- g[i-1,w1-2]+d[i,w-1]
            }
            dists <- c(g1, g2, g3)
            m <- which.min(dists)
#            message(sprintf("i=%d w=%d m=%d",i,w,m)) 
            g[i,w] <- dists[m]+d[i,w]
            bp[i,w] <- m
        }
    }
    i <- x.size
    j <- y.size
    opt <- t(as.matrix(c(i,j)))
    while (i > 1 && j > 1) {
                                        #print(c(i,j,bp[i,j_ind(i,j)]))
        #assert_that(i >= 1, i <= x.size, j_ind(i,j) >= 1, j_ind(i,j) <= wlimit)
        b <- bp[i,j_ind(i,j)]
        if (b == 1) { 
            opt <- rbind(matrix(c(i-2,j-1,i-1,j-1),ncol=2,byrow=T),opt)
            i <- i-2
            j <- j-1 
        }
        else if (b == 2) { 
            opt <- rbind(matrix(c(i-1,j-1),ncol=2),opt)
            i <- i-1
            j <- j-1 
        }
        else if (b == 3) { 
            opt <- rbind(matrix(c(i-1,j-2,i-1,j-1),ncol=2,byrow=T),opt)
            i <- i-1
            j <- j-2 
        }
        else { 
            print("Invalid backpointer")
            break 
        }
    }
    list(xsize=x.size,ysize=y.size,opt=opt,opt=opt,g=g,d=d,bp=bp)
}

toMatrix.tinyDP <- function(g) {
    g.mat <- array(Inf,dim=c(g$xsize,g$ysize))
    d.mat <- array(Inf,dim=c(g$xsize,g$ysize))
    bp.mat <- array(0,dim=c(g$xsize,g$ysize))
    for (i in 1:g$xsize) {
        for (w in -g$window:g$window) {
            j <- g$jcenter[i]+w
            if (j > 0 && j <= g$ysize) {
                g.mat[i,j] <- g$g[i,w+g$window+1]
                d.mat[i,j] <- g$d[i,w+g$window+1]
                bp.mat[i,j] <- g$bp[i,w+g$window+1]
            }
        }
    }
    return(list(g=g.mat,d=d.mat,bp=bp.mat))
}

