args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) %in% 3:4)
method_name <- if (length(args)==4) args[4] else "block_gls"
stopifnot(method_name %in% c("block_gls","adaptive_gls"))
X <- as.matrix(read.table(args[1],check.names=FALSE))
stopifnot(identical(dim(X),c(100L,113L)),all(is.finite(X)),all(X>=0))
n <- ncol(X)+1
i <- 1:(n-1)
d <- i*(n-i)/2
S <- diag(-2*d)
idx <- cbind(1:(n-2),2:(n-1))
S[idx] <- S[idx[,2:1]] <- sqrt(d[-1]*d[-length(d)])
eig <- eigen(S,symmetric=TRUE)
ends <- c(1,2,3,4,6,9,14,21,32,46,60,79,99,113)
labels <- cut(i, breaks=c(0,ends),include.lowest=TRUE)
bin <- function(x) as.numeric(rowsum(x,labels))
B <- t(apply(X,1,bin)); even <- B[seq(2,100,2),];odd <- B[seq(1,100,2),]
inv <- function(x) {
    V <- cov(x)
    if (method_name=="block_gls") return(solve((V+diag(diag(V)))/2))
    z <- scale(x)
    m <- nrow(x)
    C <- crossprod(z)/m
    target <- diag(mean(diag(C)),ncol(C))
    delta <- sum((C-target)^2)
    noise <- sum(crossprod(z^2)/m-C^2)/m
    intensity <- if (delta>0) max(0,min(1,noise/delta)) else 0
    sd_x <- apply(x,2,sd)
    solve(((1-intensity)*C+intensity*target)*outer(sd_x,sd_x))
}
P1 <- inv(even); P2 <- inv(odd); mean1 <- colMeans(odd); mean2 <- colMeans(even)
f <- function(par) {
    pop <- exp(par)
    ancestral <- pop[1]; current <- pop[2]; T <- pop[3]
    old <- drop(eig$vectors %*% (exp(eig$values*T/(2*current))*crossprod(eig$vectors,sqrt(d)/i)))/sqrt(d)
    model <- 4*3.42e-8*1e6*bin(current/i+(ancestral-current)*old)
    a <- mean1-model;b <- mean2-model
    drop(crossprod(a,P1%*%a)+crossprod(b,P2%*%b))
}
p <- optim(log(c(20000,13000,3000)),f,method="Nelder-Mead",control=list(maxit=5000,reltol=1e-14))
z <- exp(p$par)
ref <- read.delim(args[2]); ref <- subset(ref,phase==args[3] & method==method_name)
stopifnot(nrow(ref)==1L,all(is.finite(z)),is.finite(ref$generations),is.finite(ref$fraction),p$convergence==0,abs(z[3]-ref$generations)<0.02,abs(z[2]/z[1]-ref$fraction)<1e-6)
cat(args[3],sprintf("%.5f %.8f",z[3],z[2]/z[1]),"\n")
