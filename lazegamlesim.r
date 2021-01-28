# The next sets the number of simulations to do.
# It is set to 2 for illustration purposes
nsims <- 20

library(ergm)

model.f <- as.formula(g ~ edges + gwesp(alpha.mle,fixed=TRUE)+
                      nodecov("seniority") +
                      nodecov("practice") + match("practice") + 
                      match("gender") + 
                      match("office"))
#
# Lazega's lawyers mutual collaboration network
#

# Next for fit to original data
load(file="lazega.fit.RData")
# Next for fit to increased transitivity version
#load(file="lazega.fit.high.RData")

theta0.mle <- work.fit$coef
theta0.mcmc <- work.fit$MCMCtheta

alpha.mle
theta0.mle
theta0.mcmc
#
# Make curved.parms
#
g <- work.fit$network
trms <- ergm.getterms(model.f)
xobs0<-summary.statistics.network(model.f)
xobs0
model.f
g <- try(as.network(eval(trms[[2]],sys.parent())))                  
m <- ergm.getmodel(model.f, g, drop=FALSE)
Clist <- ergm.Cprepare(g, m)                                            
m.expanded <- ergm.getmodel(model.f, g, drop=FALSE, expanded=TRUE)    
Clist.expanded <- ergm.Cprepare(g, m.expanded)                
ecurved <- ergm.curved(theta0.mle,m,m.expanded,g,NULL)

statsmatrix <- work.fit$sample
#
sm <- statsmatrix
print(apply(sm,2,summary.statsmatrix.ergm),scipen=7)
pdf(file="lazega.mle.pdf")
#
source("lazega.llik.r")
#
# Calculate the MLEs
#
lmle <- array(NA,dim=c(nsims,6*length(theta0.mle)+3,2))
isim <- round(seq(from=1, to=nrow(sm), length=nsims))
dimnames(lmle) <- list(NULL,c("llratio",
                       names(theta0.mle),
                       paste("se",names(theta0.mle),sep="."),
                       paste("mv",names(theta0.mle),sep="."),
                       paste("mv.se",names(theta0.mle),sep="."),
                       "mahalanobis", "mv.mahalanobis",       
                       paste("obs",names(theta0.mle),sep="."),
                       paste("obs0",names(theta0.mle),sep=".")
                             ),
                       c("MLE","MPLE"))
l0<-ergm.estimate(model=m, theta0=theta0.mcmc, statsmatrix=sm,
                  parms.curved=ecurved$parms.curved, calc.mcmc.se=F)
print(l0$coef)
print(xobs0)
ll0 <- llikratio(theta=l0$coef, theta0=theta0.mcmc, statsmatrix=sm)

times<-NULL 
for( i in (1:(nsims))){
 g <- lsim$networks[[i]]
 xobs <- summary.statistics.network(model.f)
 print(xobs)
 t1<-proc.time()
 fitmple<-ergm(model.f, MPLEonly=TRUE, algorithm.control=list(drop=FALSE))
 t2<-proc.time()
 if(all(xobs>0) & xobs["edges"] > 10){
  xsim <- sweep(sm, 2, xobs0-xobs,"-")
  t3<-proc.time()   
  fitmle <- ergm(model.f,
              theta0=0.9*theta0.mle+0.1*fitmple$coef,
              maxit=2,
              algorithm.control=list(steplength=0.7, drop=FALSE),
#             parallel=2,
#             burnin=150000, MCMCsamplesize=5000, interval=10000
              burnin=500000, MCMCsamplesize=5000, interval=2000
             )
  fitmle <- ergm(model.f,
              theta0=fitmle$coef,
              maxit=2,
              algorithm.control=list(steplength=0.7, drop=FALSE),
#             burnin=150000, MCMCsamplesize=50000, interval=15000
              burnin=500000, MCMCsamplesize=10000, interval=3000
             )
  t4<-proc.time()
  times<-rbind(times,c(t2-t1,t4-t3))
  ll <- llikratio(theta=fitmle$coef, theta0=theta0.mcmc, statsmatrix=sm)
  lmv <- apply(fitmle$sample,2,mean)+xobs-xobs0
  lmvse <- sqrt(diag(solve(fitmle$covar)))
  mlese <- sqrt(diag(fitmle$covar))
  mlemvcov<-cov(fitmle$sample)
  mlemah<-ergm.mahalanobis(fitmle$coef,theta0.mle,fitmle$covar,inverted=F)
  mlemvmah<-ergm.mahalanobis(lmv,xobs0,mlemvcov,inverted=F)

  mplese <- sqrt(diag(fitmple$covar))
  lmple<-fitmple$coef
  llmple <- llikratio(theta=lmple, theta0=theta0.mcmc, statsmatrix=sm)
  fit <- ergm(model.f,
              theta0=lmple,
              maxit=1,
              algorithm.control=list(drop=FALSE),
              parallel=2,
              burnin=1500000, MCMCsamplesize=10000, interval=150,
              estimate=FALSE
             )
  lmplemv <- apply(fit$sample,2,mean)+xobs-xobs0
  lmplemvse <- sqrt(diag(solve(fitmple$covar)))
  mplemvcov <- cov(fit$sample) 
  mplemah <- ergm.mahalanobis(lmple, theta0.mle, fitmple$covar, inverted=F)
  mplemvmah <- ergm.mahalanobis(lmplemv, xobs0, mplemvcov, inverted=F) 
  if(all(xobs>0) & xobs["edges"] > 10){
   lmle[i,,1] <- c(ll$llik-ll0$llik, fitmle$coef, mlese, lmv, lmvse, mlemah,
                   mlemvmah, xobs, xobs0)
  }
  lmle[i,,2] <- c(llmple$llik-ll0$llik, lmple, mplese, lmplemv,
                  lmplemvse, mplemah, mplemvmah, xobs, xobs0)
  cat("\n")
  print(paste("Simulation:",
  paste(c(i,format(c(ll$llik,llmple$llik),digits=2)),collapse=" ")))
 }else{
  cat("\n")
  print(paste("Simulation:",i,"degenerate",collapse=" "))
 }
}
save(theta0.mle, lmle, file="lazega.mle.new.RData")
print(apply(times,2,mean))
