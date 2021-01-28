llikratio<-function(
   theta,  # theta to calculate the log-lik (relative to theta0)
   theta0, # theta0 from sim
   statsmatrix # sim stats matrix
   )
{
   probs <- rep(1/nrow(statsmatrix),nrow(statsmatrix))
#
# Map functions
#
  identity.map<- function(x,y){x} 
  geosdeg.map<- function(x,y){x[1]*(exp(-(x[2])*seq(along=y))-1)} 
  identity.grad<- function(x,y){diag(x-x+1,ncol=length(x))} 
  geosdeg.grad<- function(x,y){
    cbind(exp(-(x[2])*seq(along=y))-1,
          -x[1]*seq(along=y)*exp(-x[2]*seq(along=y))
         )
                              } 
#
# Identify curved parameters
#
  eta.parm <- list(from=matrix(rep(TRUE, length(theta0)),ncol=1),
                     to=matrix(rep(TRUE, ncol(statsmatrix)),ncol=1),
                     map=vector(mode="list"), mapgrad=vector(mode="list")
#                    map="identity.map", mapgrad="identity.grad"
                  )
  eta.parm$map[[1]] <- identity.map
  eta.parm$mapgrad[[1]] <- identity.grad
  names(eta.parm$map) <- "identity"
  names(eta.parm$mapgrad) <- "identity"
  dimnames(eta.parm$from) <- list(names(theta0),"identity") 
  dimnames(eta.parm$to) <- list(colnames(statsmatrix),"identity") 
  geod <- grep("geosdegree#",colnames(statsmatrix))     
  if(length(geod)>0){
    geodf <- rep(FALSE, ncol(statsmatrix))
    geodf[geod] <- TRUE
    eta.parm$from <- cbind(eta.parm$from,names(theta0)=="geosdegree")
    loc.geosdegree <- match("geosdegree",names(theta0))
    colnames(eta.parm$from)[loc.geosdegree] <- "geosdegree"
#
    eta.parm$from[loc.geosdegree,     ncol(eta.parm$from)] <- TRUE
    eta.parm$from[nrow(eta.parm$from),ncol(eta.parm$from)] <- TRUE
    eta.parm$from[nrow(eta.parm$from),                  1] <- FALSE
    eta.parm$to   <- cbind(eta.parm$to,geodf)
    eta.parm$to[nrow(eta.parm$from),                   1] <- FALSE
    eta.parm$to[nrow(eta.parm$from), ncol(eta.parm$from)] <- FALSE
    colnames(eta.parm$to  )[loc.geosdegree] <- "geosdegree"
    eta.parm$map[[2]]  <- geosdeg.map
    eta.parm$mapgrad[[2]] <- geosdeg.grad
    names(eta.parm$map)[2] <- "geosdegree"
    names(eta.parm$mapgrad)[2] <- "geosdegree"
    eta.parm$from[,1] <- eta.parm$from[,1] & !eta.parm$from[,2] 
    eta.parm$to[,1]   <- eta.parm$to[,1] & !eta.parm$to[,2] 
  }
#
  av <- apply(sweep(statsmatrix,1,probs,"*"), 2, sum)
  xsim <- sweep(statsmatrix, 2, av,"-")
  xobs <- - av
#
  trace <- 0
#
# Log-Likelihood and gradient functions
#
#   complete data code
#
    llik.fun <- function(x, x0, xobs, xsim, probs, eta.parm){
#
#    eta transformation
#
     eta <- rep(0,length(xobs))
     for(i in 1:length(eta.parm$map)){
      eta[eta.parm$to[,i]] <- eta.parm$map[[i]](
        x[eta.parm$from[,i]],eta[eta.parm$to[,i]])
     }
#
     eta0 <- rep(0,length(xobs))
     for(i in 1:length(eta.parm$map)){
      eta0[eta.parm$to[,i]] <- eta.parm$map[[i]](
        x0[eta.parm$from[,i]],eta[eta.parm$to[,i]])
     }
     x <- eta-eta0
     basepred <- xsim %*% x
# alternative based on log-normal approximation
     mb <- sum(basepred*probs)
     vb <- sum(basepred*basepred*probs) - mb*mb
     llr <- sum(xobs * x) - (mb + 0.5*vb)
     llr
    }
    llik.grad <- function(x, x0, xobs, xsim, probs, 
                          eta0, eta.parm){
#
#    eta transformation
#
     eta <- rep(0,length(xobs))
     etagrad <- matrix(0,ncol=length(x),nrow=length(xobs))
     for(i in 1:length(eta.parm$map)){
#     eta[eta.parm$to[,i]] <- do.call(eta.parm$map[i],
      eta[eta.parm$to[,i]] <- eta.parm$map[[i]](
        x[eta.parm$from[,i]],eta[eta.parm$to[,i]])
#     etagrad[eta.parm$to[,i],eta.parm$from[,i]] <- 
#       do.call(eta.parm$mapgrad[i],
      etagrad[eta.parm$to[,i],eta.parm$from[,i]] <- eta.parm$mapgrad[[i]](
        x[eta.parm$from[,i]],eta[eta.parm$to[,i]])
     }
     x <- eta-eta0
     prob <- max(xsim %*% x)
     prob <- probs*exp(xsim %*% x - prob)
     prob <- prob/sum(prob)
     E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
# 
#    The next lines are for the Hessian which optim does not use
# 
#    vtmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
#    V <- t(vtmp) %*% vtmp
#    list(gradient=xobs-E,hessian=V)
      (xobs-E)%*%etagrad
    }
#
  llik <- llik.fun(x=theta, 
    x0=theta0, xobs=xobs, xsim=xsim, probs=probs,
    eta.parm=eta.parm
  )
# llik <- llik + sum(xobs * eta0)
#
# llikgradient <- llik.grad(x=theta, x0=theta0, xobs=xobs, xsim=xsim,
#                       probs=probs,
#                       eta.parm=eta.parm
#                          )
  names(theta) <- names(theta0)
#
# Calculate the (global) log-likelihood
# without the + log |\cal X| added in ergm.statseval
#
  c0  <- llik.fun(x= theta, x0=theta0, xobs=xobs-xobs,
                  xsim=statsmatrix, probs=probs,
                  eta.parm=eta.parm
                 )
  c01 <- llik.fun(x=theta0-theta0, x0=theta0, xobs=xobs-xobs,
                  xsim=statsmatrix, probs=probs,
                  eta.parm=eta.parm
                 )
#
  loglik <- c0 - c01

# Output results as ergm-class object
  list(coef=theta, 
      mle.lik=loglik, llik=llik,
      MCMCtheta=theta0, loglikelihoodratio=c0,
      c01=c01
      )
}
simplellikratio<-function(
   theta,  # theta to calculate the log-lik (relative to theta0)
   theta0, # theta0 from sim
   statsmatrix # sim stats matrix
   )
{
  probs <- rep(1/nrow(statsmatrix),nrow(statsmatrix))
  av <- apply(sweep(statsmatrix,1,probs,"*"), 2, sum)
  xsim <- sweep(statsmatrix, 2, av,"-")
  xobs <- - av
#
  trace <- 0
#
# Log-Likelihood and gradient functions
#
#   complete data code
#
    llik.fun <- function(x, x0, xobs, xsim, probs){
     x <- x-x0
#    aaa <- sum(xobs * x) - log(sum(probs*exp(xsim %*% x)))
#    aaa <- (xsim %*% x) - sum(xobs * eta)
     aaa <- (xsim %*% x) - sum(xobs * x)
     if(any(is.na(aaa))){return(-10^8)}
     a1 <- aaa >  650
     a2 <- aaa < -650
     a3 <- !a1 & !a2
     bbb <- 0
     if(sum(a1) > 0){bbb <- bbb + exp( 650)*sum(probs[a1]*exp(aaa[a1]-650),
       na.rm=TRUE) }
     if(sum(a2) > 0){bbb <- bbb + exp(-650)*sum(probs[a2]*exp(aaa[a2]+650),
       na.rm=TRUE) }
     if(sum(a3) > 0){bbb <- bbb + sum(probs[a3]*exp(aaa[a3]),
       na.rm=TRUE) }
# 
#    This is the log-likelihood (and not its negative)
# 
#    -aaa
# 
#    This is the log-likelihood
# 
     bbb <- -log( bbb )
     bbb
    }
    llik.grad <- function(x, x0, xobs, xsim, probs){
     x <- x-x0
     prob <- max(xsim %*% x)
     prob <- probs*exp(xsim %*% x - prob)
     prob <- prob/sum(prob)
     E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
# 
#    The next lines are for the Hessian which optim does not use
# 
#    vtmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
#    V <- t(vtmp) %*% vtmp
#    list(gradient=xobs-E,hessian=V)
      (xobs-E)%*%etagrad
    }
#
  llik <- llik.fun(x=theta, 
    x0=theta0, xobs=xobs, xsim=xsim, probs=probs
  )
# llik <- llik + sum(xobs * eta0)
#
# llikgradient <- llik.grad(x=theta, x0=theta0, xobs=xobs, xsim=xsim,
#                       probs=probs,
#                       eta.parm=eta.parm
#                          )
  names(theta) <- names(theta0)
#
# Calculate the (global) log-likelihood
# without the + log |\cal X| added in ergm.statseval
#
  c0  <- llik.fun(x= theta, x0=theta0, xobs=xobs-xobs,
                  xsim=statsmatrix, probs=probs
                 )
  c01 <- llik.fun(x=theta0-theta0, x0=theta0, xobs=xobs-xobs,
                  xsim=statsmatrix, probs=probs
                 )
#
  loglik <- c0 - c01

# Output results as ergm-class object
  list(coef=theta, 
      mle.lik=loglik, llik=llik,
      MCMCtheta=theta0, loglikelihoodratio=c0,
      c01=c01
      )
}
meanvalue<-function(
   theta,  # theta to calculate the mean (based on a sample from theta0)
   theta0, # theta0 from sim
   statsmatrix # sim stats matrix
   )
{
#
# Map functions
#
  identity.map<- function(x,y){x} 
  geosdeg.map<- function(x,y){x[1]*(exp(-(x[2])*seq(along=y))-1)} 
  identity.grad<- function(x,y){diag(x-x+1,ncol=length(x))} 
  geosdeg.grad<- function(x,y){
    cbind(exp(-(x[2])*seq(along=y))-1,
          -x[1]*seq(along=y)*exp(-x[2]*seq(along=y))
         )
                              } 
#
# Identify curved parameters
#
  eta.parm <- list(from=matrix(rep(TRUE, length(theta0)),ncol=1),
                     to=matrix(rep(TRUE, ncol(statsmatrix)),ncol=1),
                     map=vector(mode="list"), mapgrad=vector(mode="list")
#                    map="identity.map", mapgrad="identity.grad"
                  )
  eta.parm$map[[1]] <- identity.map
  eta.parm$mapgrad[[1]] <- identity.grad
  names(eta.parm$map) <- "identity"
  names(eta.parm$mapgrad) <- "identity"
  dimnames(eta.parm$from) <- list(names(theta0),"identity") 
  dimnames(eta.parm$to) <- list(colnames(statsmatrix),"identity") 
  geod <- grep("geosdegree#",colnames(statsmatrix))     
  if(length(geod)>0){
    geodf <- rep(FALSE, ncol(statsmatrix))
    geodf[geod] <- TRUE
    eta.parm$from <- cbind(eta.parm$from,names(theta0)=="geosdegree")
    loc.geosdegree <- match("geosdegree",names(theta0))
    colnames(eta.parm$from)[loc.geosdegree] <- "geosdegree"
#
    eta.parm$from[loc.geosdegree,     ncol(eta.parm$from)] <- TRUE
    eta.parm$from[nrow(eta.parm$from),ncol(eta.parm$from)] <- TRUE
    eta.parm$from[nrow(eta.parm$from),                  1] <- FALSE
    eta.parm$to   <- cbind(eta.parm$to,geodf)
    eta.parm$to[nrow(eta.parm$from),                   1] <- FALSE
    eta.parm$to[nrow(eta.parm$from), ncol(eta.parm$from)] <- FALSE
    colnames(eta.parm$to  )[loc.geosdegree] <- "geosdegree"
    eta.parm$map[[2]]  <- geosdeg.map
    eta.parm$mapgrad[[2]] <- geosdeg.grad
    names(eta.parm$map)[2] <- "geosdegree"
    names(eta.parm$mapgrad)[2] <- "geosdegree"
    eta.parm$from[,1] <- eta.parm$from[,1] & !eta.parm$from[,2] 
    eta.parm$to[,1]   <- eta.parm$to[,1] & !eta.parm$to[,2] 
  }
#
  eta <- rep(0,ncol(statsmatrix))
  for(i in 1:length(eta.parm$map)){
   eta[eta.parm$to[,i]] <- eta.parm$map[[i]](
     theta[eta.parm$from[,i]],eta[eta.parm$to[,i]])
  }
#
  eta0 <- rep(0,ncol(statsmatrix))
  for(i in 1:length(eta.parm$map)){
   eta0[eta.parm$to[,i]] <- eta.parm$map[[i]](
     theta0[eta.parm$from[,i]],eta[eta.parm$to[,i]])
  }
  x <- eta-eta0
  basepred <- statsmatrix %*% x
  mbasepred <- max(basepred)
  prob <- exp(basepred - mbasepred)
  prob <- prob/sum(prob)
  E <- apply(sweep(statsmatrix, 1, prob, "*"), 2, sum)
# mb <- sum(basepred*probs)
# vb <- sum(basepred*basepred*probs) - mb*mb
# llr <- (mb + 0.5*vb)
# delta <- x-x
# delta[3] <- 0.0001
# basepred <- statsmatrix %*% (x+delta)
# mb <- sum(basepred*probs)
# vb <- sum(basepred*basepred*probs) - mb*mb
# llr <- (llr - (mb + 0.5*vb))/0.0001
# Output results as ergm-class object
  E
}
ergm.mahalanobis <- function(x, center, cov, inverted=FALSE, ...)
{
    x <- matrix(x, ncol=length(x))
    x <- sweep(x, 2, center)
    cov <- robust.inverse(cov, ...)
    retval <- rowSums((x%*%cov) * x)
    names(retval) <- rownames(x)
    retval
}
