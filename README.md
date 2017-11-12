# tinyDP: A fast dynamic time warping for R
## Example
```R:
> library(tinyDP)
> x<-matrix(c(1,2,3,4,5,6,7,8,9,7,8,9),byrow=T,ncol=3)
> x
     [,1] [,2] [,3]
[1,]    1    2    3
[2,]    4    5    6
[3,]    7    8    9
[4,]    7    8    9
> y<-matrix(c(1,2,3,4,5,6,4,5,6,7,8,9),byrow=T,ncol=3)
> y
     [,1] [,2] [,3]
[1,]    1    2    3
[2,]    4    5    6
[3,]    4    5    6
[4,]    7    8    9
> dp<-tinyDP(x,y)
> dp$opt
     [,1] [,2]
[1,]    1    1
[2,]    2    2
[3,]    2    3
[4,]    3    4
[5,]    4    4
>
```

