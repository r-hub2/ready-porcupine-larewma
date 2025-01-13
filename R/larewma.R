larewma <- function(F, S, lambda, arl0, h, only.mean=FALSE, horiz=20000, A=3, B=0.6, csi = 5,
                    burnin = 1000, iter = 50000, seed = 11226457 ) {
    F <- as.matrix(F)
    S <- as.matrix(S)
    n <- NROW(F)
    p <- NCOL(F)
    if ( (NROW(S) != n) || (NCOL(S) != n) ) stop("S has wrong dimension")
    allstat <- 1+(p+2)+p*p+p+max(p,n)+3+7*p+2*p*p
    npar <- 2 + n*n + p*p + n*p + 2*(p+1)
    si <- solve(S)
    sif <- si %*% F
    fsif <- crossprod(F,sif)
    ipar <- c(1,n,allstat,0,6,npar,allstat,p,only.mean)
    par <- c(5,lambda,si,fsif,t(sif),rep(0,2*(p+1)))
    stat0 <- c(0,rep(0,p+1),1,rep(0,allstat-p-3))
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        seed.keep <- get(".Random.seed", envir = .GlobalEnv, 
                         inherits = FALSE)
        on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
    }
    set.seed(seed)
    par <- .C("gglarcompmom",as.integer(horiz),as.double(t(chol(S))),
              as.integer(ipar),par=as.double(par),as.double(stat0),PACKAGE="larewma")$par
    if ( missing(arl0) ) {
        if (missing(h)) stop("Either arl0 or h must be not missing")
        arl0 <- NA
        par[1] <- h
    } else {
        opt <- c(arl0,csi,burnin,iter,A,B)
        par <- .C("ggmatcharl0one",
                  as.integer(ipar), par=as.double(par),as.double(stat0),
                  as.double(t(chol(S))), as.double(opt), DUP=FALSE,PACKAGE="larewma")$par
    }
    ab <- t(matrix(par[(3+n*n+p*p+n*p):length(par)],p+1))
    rownames(ab) <- c("a","b")
    colnames(ab) <- 1:(p+1)    
    list(F=F, S=S, lambda=lambda, arl0=arl0, h=par[1], ab=ab,
         only.mean=only.mean, horiz=horiz,
         A=A, B=B, csi=csi, burnin=burnin, iter=iter, seed=seed,
         ipar=ipar , par=par , stat0=stat0)
}


larewma.simrl <- function(chart , nsim  , tau=0 , delta=0 , omega=chart$S ) {
    n <- chart$ipar[2]
    if ( (NROW(omega) != n) || (NCOL(omega) != n) )
        stop("omega has wrong dimension")
    if ( length(delta) == 1 ) delta <- rep( delta , chart$ipar[2] )
    if ( length(delta) != n ) stop("delta has a wrong length")
    .C("ggsimrlss",
       as.integer(chart$ipar), as.double(chart$par), as.double(chart$stat0),
       as.double(t(chol(chart$S))), as.integer(tau) , as.double(delta),
       as.double(t(chol(omega))), as.integer(c(nsim,.Machine$integer.max)),
       rl = integer(nsim), DUP=FALSE, PACKAGE="larewma")$rl
}




larewma.apply <- function(chart , y, mu=0, plot=TRUE ) {
    y <- t(as.matrix(y))
    ny <- NCOL(y)
    n <- NROW(y)
    if ( chart$ipar[2] != n )
        stop("Number of variables in y and chart do not match")
    if ( !missing(mu) ) y <- y - mu
    ans <- .C("ggcomputecharts",
              as.integer(chart$ipar), as.double(chart$par),
              as.integer(ny), as.double(y) , as.double(chart$stat0),
              ans = double(ny * chart$ipar[3]) , DUP = FALSE , PACKAGE="larewma")$ans
    ans <- ts( matrix( ans , ncol = chart$ipar[3] , byrow = TRUE ) )
    if ( plot ) {
        plot(ans[,1], type="b", ylim=c(min(ans[,1]),1.05*max(ans[,1],chart$par[1])),
             xlab="t",ylab=expression(W[t]))
        abline(h=chart$par[1], lty="dotted")
    }
    p <- chart$ipar[8]
    js <- 2:(p+2)
    jb <- (p+4):(p+3+p*p)
    select <- function(i) {
        w <- as.numeric(ans[i,1])
        d <- cbind(ans[i,js],rbind(matrix(ans[i,jb],p,byrow=TRUE),0))
        rownames(d) <- 1:(p+1)
        colnames(d) <- c("(S-a)/b",paste("beta",1:p,sep=""))
        list(t=i,W=w,stat=d)
    }
    alarms <- which(ans[,1]>chart$par[1])
    if ( length(alarms) ) {
            lapply(alarms,select)
    } else {
        list()
    }
}




