#' @title bdribs (Bayesian Detection of Risk using Inference on Blinded Safety data)
#'
#' @description Bayesian Detection of Risk using Inference on Blinded Safety data
#' @author Saurabh Mukhopadhyay
#'
#' @param y observed pooled events (combined active and control group) e.g., y = 20
#' @param pyr total payr exposure (combined active and control group) e.g., pyr = 2000
#' @param bg.events background (historical) events for the control group e.g., bg.events = 5
#' @param bg.pyr background (historical) pyr exposure for the control group e.g., bg.pyr = 1000
#' @param bg.rate when specified used as the true background rate for the control group and ignores bg.events and bg.pyr, default: bg.rate=NULL
#' @param k allocation ratio of treatment vs. control group, default: k=1
#' @param p.params paramaters of beta prior of p (used only when inf.type = 1 or = 2); default: p.param= list(a=1,b=1). See details below.
#' @param r.params paramaters of log-normal prior of r (used only when p.params=NULL); default: r.param= list(mu=0,sd=2). See details below.
#' @param adj.k when TRUE adjusts the prior specification for k >1 (or for k <1), default: adj.k = FALSE . See deatil below.
#' @param mc.params contains detials of MCMC parameters, default: mc.params=list(burn=1000, iter=10000, nc=2)
#' @param inf.type indicate inference type, default: inf.type =1 (gives conditional inference for fixed background rate). See deatil below.
#' @param plots indicates whether standard plots to be generated, default: plots= TRUE
#' @param prnt indicates whether inputs  to be printed, default: prnt= TRUE
#'
#' @return returns a dataframe of MCMC output from the posterior distribution for parameters of interests
#'
#' @details This 'bdribs' package obtains Bayesian inferences on blinded pooled safety data  ...
#' @details Values of p.params are used to specify a beta prior for p - default is Jeffreys non-informative prior: Beta(a=0.5,b=0.5).
#' @details If inf.type=1, then conditional posterior inference on r is obtained for
#'      a given fixed values of del0 = bg.rate = bg.events/bg.pyr.
#' @details If inf.type=2, then an average (marginal) Bayesian inference on r is obtained with respect to a prior on del0, where
#'     del0 ~ Gamma(bg.events, bg.pyr).
#' @details If prior on r must be specified directly it can be done by using a log-normal prior. To do that, p.params must be set to NULL and
#'    and then r.params should be specifed as a list to supply mean and sd of the lognormal. For example, to have a lognormal prior
#'    with log-mean 0 and log-sd = 2, we should set r.params = list(mu=0, sd=2) and p.params=NULL.
#' @details when adj.k = TRUE, and k is not 1 (that is, allocation ratio is not 1:1), then
#' a non-informative prior such as (beta(.5, .5) is first specified on p, assuming equal allocation ratio and then adjusted for the give k.
#' When adj.k = F, then no such adjustment is made on the prior for p. Note that no such adjustments
#' needed if prior on r is directly specified (as discussed above). However, it is always difficult to
#' specify a non-informative prior on r and therefore a a prior on p with adj.k =T is recommended in most cases.
#'
#' @examples ## Sample calls
#'     #run 1: simple case with a fixed background rate of 0.45 per 100 pyr.
#'     bdribs(y=5,pyr=500, bg.rate=0.0045,k=2)
#'
#'     #run 2: same as run 1; here bg.rate gets computed as bg.events/bg.pyr
#'     bdribs(y=5,pyr=500, bg.events = 18, bg.pyr = 4000, k=2)
#'
#'     # run3: when inf.type = 2, uses a Gamma distribution for del0; e.g. here Gamma(18, 4000)
#'     bdribs(y=5,pyr=500, bg.events = 18, bg.pyr = 4000, k=2, inf.type = 2)
#'
#'     #run4: similar to run1, but instead of default p~u(0,1) using p~beta(.5,.5)
#'     bdribs(y=5,pyr=500, , bg.rate=0.0045,k=2, p.params=list(a=.5,b=.5))
#'
#'     #run5: similar to run1, but instead of default p ~ beta(.5,.5) using r ~ lognormal(mu=0,sd=2)
#'     bdribs(y=5,pyr=500, , bg.rate=0.0045,k=2, p.params= NULL, r.params=list(mu=0,sd=2))
#'
#' @import rjags
#' @importFrom  stats quantile density update dbeta dlnorm pbeta plnorm qgamma
#' @importFrom graphics abline axis box contour plot title arrows par points
#' @export

bdribs = function(y, pyr, bg.events, bg.pyr, bg.rate=NULL, k=1,
                  p.params=list(a=1, b=1), r.params=list(mu=0,sd=2), adj.k=FALSE,
                  mc.params=list(burn=1000, iter=10000, nc=2), inf.type=1, plots=TRUE, prnt=TRUE)
{
  # require(rjags) # should be laoded when the package is being laoded

  if (y < 0 | y-round(y) > 1.e-7 ) stop("value of y must be a non-negative integer")
  else y = round(y)
  if (pyr <=0) stop("value of pyr must be positive")
  if (k <=0) stop("value of k must be positive")

  if (!is.null(bg.rate)) # when background rate is directly provided then assuming that is the true rate and using conditional inference given that rate
  {
    if (bg.rate<=0) stop ("bg.rate must be a positive value indicating historical event rate per patient year")
    if (inf.type!=1) {
      print("switching to conditional inference (inf.type=1) as bg.rate is provided")
      inf.type=1
    }
    bg.pyr=1000
    bg.events=bg.rate*bg.pyr
  }
  else if (inf.type==1)
  {
    if (bg.events <=0 | bg.pyr <=0 ) stop ("both bg.events and bg.pyr must be positive")
    bg.rate = bg.events/bg.pyr
  }

  #if (bg.pyr<1000 & inf.type =2) warning("background rate is not precise enough: consider using inf.type=1")
  if (mc.params$nc !=2) stop ("currently the jags are implemented with nc=2 only - i.e. number chains =2" )
  if (!(inf.type %in% c(1,2))) stop ("valid input values of inf.type are 1 or 2 - see more in 'Details' in bdribs help." )

  if (!is.null(p.params))
  {
    if (p.params$a <= 0 | p.params$b <= 0) stop("values of (a,b) in p.params mist be both positive")
    params = p.params
    prior.string = paste("prior is on p ~ beta(a=",p.params$a,", b=",p.params$a,")",sep = "")
    rl = seq(0,7, length.out = 100)
    fr = dbeta(rl*k/(1+rl*k),p.params$a,p.params$b)/(1+rl*k)^2
    if(adj.k) fr = dbeta(rl/(1+rl),p.params$a,p.params$b)/(1+rl)^2 # if adj.k = T

  }
  else
  {
    if (is.null(r.params$mu) | r.params$sd <= 0) stop("mu in r.params must be a real number AND sd must be positive")
    params = r.params
    prior.string = paste("prior is on r ~ log-normal(mu=",r.params$mu,", sd=",r.params$sd,")",sep = "")
    rl = seq(0,7, length.out = 100)
    fr = dlnorm(rl,r.params$mu,r.params$sd)
  }

  if (adj.k & k != 1 & !is.null(p.params) )
    prior.string = paste(prior.string,"(with k =1 and then adjusted for k = ",k,")", sep="")


  para = core.bdribs(y, pyr, bg.events, bg.pyr, k=k, params=params, adj.k=adj.k, mc.params=mc.params, inf.type=inf.type)
  if(plots){
    para.names = names(para)

    #plot(para$r, type="l",col="gray", xlab="Iterations", ylab="", cex.main=0.8,
    #  main=paste("Trace  of ", para.names[3], "\n (MCMC Samples from Posterior Dsitribution)", sep=""))
    if (inf.type==1) {
      pttl = paste("Prior and posterior densities of r \n(y = ",y,", pyr = ", pyr, " and fixed bg.rate = ",bg.rate,")\n",sep="")
      ftnt = prior.string
    }
    else {
      pttl = paste("Prior and posterior densities of r \n(y = ",y,", pyr = ", pyr, " and a mean bg.rate = ",round(bg.events/bg.pyr,4),")", sep="")
      ftnt = paste(prior.string, "\nPrior for control group event rate ~ Gamma (d0, d1) with (d0, d1) = (",bg.events,",", bg.pyr,")",sep="")
    }
    drplot(x=para$r, xgrd=0:7, xlm=c(0,7), ygrd=seq(0,.9,.1), clr="brown", rlst=rl, frl=fr,
           ttl = pttl, ftn = ftnt)
  }
  if(prnt){
    xc = para$r
    cuts = c(1,1.5,2)
    pprobs=c()
    for (j in 1:length(cuts))  pprobs = c(pprobs, mean(xc>cuts[j]))
    prprobs=c()
    for (j in 1:length(cuts))  prprobs = c(prprobs,comp.prior(c=cuts[j], k=k, adj.k=adj.k, p.params=p.params, r.params=r.params))
    bf = round((pprobs/(1-pprobs))/ (prprobs/(1-prprobs)),2)

    smry = round(c(mean(xc), quantile(xc, probs = c(0.05, 0.95)), pprobs, bf),3)
    smry1 = data.frame(smry[1], smry[2], smry[3], smry[4], smry[5], smry[6], bf[1], bf[2], bf[3])
    names(smry) = c("mean", "5%", "95%", paste("P[r",">",cuts,"]",sep=""),paste("BF[r",">",cuts,"]",sep=""))
    names(smry1) = c("mean", "5%", "95%", paste(" P[r",">",cuts,"]",sep=""),paste(" BF[r",">",cuts,"]",sep=""))
    row.names(smry1)=""


    for(j in 1:79) cat("=")
    cat("\n")
    cat("Printing posterior inference of relative risk (r):\n")
    cat(prior.string, "\n")
    cat(paste("y = ", y, ", pyr = ", pyr, ", k = ", k, ", bg.rate = ", round(bg.events/bg.pyr,5),", inf.type = ", inf.type,"\n",sep=""))
    for(j in 1:79) cat("=")
    cat("\n")
    print(smry1)
    for(j in 1:79) cat("=")
    cat("\n")
  }
  return(invisible(para))
}

comp.prior = function(c, k, adj.k=FALSE, p.params, r.params) # computes prior prob of r > c
{
  if(!is.null(p.params))
  {
    if (!adj.k) ans = 1-pbeta(1-1/(1+k*c), p.params$a, p.params$b)
    else ans = 1-pbeta(1-1/(1+c), p.params$a, p.params$b)
  }
  else ans = 1- plnorm(c, r.params$mu, r.params$sd)
  return(round(ans,4))
}

core.bdribs = function(y, pyr, bg.events, bg.pyr, k, params, mc.params, inf.type, adj.k)
{
  # define the model
  model1="model {
    y ~dpois(ld)
    ld <- del0*pyr/(1-pk)/(k+1)
    r <- pk/(1-pk)/k # note that this is = p/(1-p) as it should be since r is free of k
    pk <- k*p/(k*p+ (1-p))
    p ~ dbeta(a,b)
    del0 ~ dgamma(d0+0.001,d1+0.001)
  }"

  # define the model where prior on p is not adjusted for k
  model1.noadj="model {
    y ~dpois(ld)
    ld <- del0*pyr/(1-p)/(k+1)
    r <- p/(1-p)/k
    p ~ dbeta(a,b)
    del0 ~ dgamma(d0+0.001,d1+0.001)
  }"

  # define an alternative model -- when prior is specified via a lognormal distribution of parameter r
  model2="model {
    y ~dpois(ld)
    ld <- del0*pyr/(1-pk)/(k+1)
    pk <- k*p/(k*p+ (1-p))
    p <- r/(1+r)
    r ~ dlnorm(a,b)
    del0 ~ dgamma(d0+0.001,d1+0.001)
  }"

  if (!is.null(params$a)) {
    modelstring = model1 # uses beta prior on p
    if (!adj.k)  modelstring = model1.noadj # uses beta prior on p - but not adjusting for k (old method)
    if (inf.type==1) { #  ... conditional inference with a beta prior on p
      jdata = list(y=y, pyr = pyr, k=k, a=params$a, b=params$b, d0=bg.events, d1=bg.pyr, del0=bg.events/bg.pyr)
      inits1 <- list(p=.5,  .RNG.name="base::Super-Duper", .RNG.seed=1)
      inits2 <- list(p=.45, .RNG.name="base::Wichmann-Hill", .RNG.seed=2)
    }
    if (inf.type==2) {#... unconditional inference with a beta prior on p and gamma prior on del0
      jdata = list(y=y, pyr = pyr,  k=k, a=params$a, b=params$b, d0=bg.events, d1=bg.pyr)
      inits1 <- list(p=.5, del0=bg.events/bg.pyr,  .RNG.name="base::Super-Duper", .RNG.seed=1)
      inits2 <- list(p=.45, del0=bg.events/bg.pyr, .RNG.name="base::Wichmann-Hill", .RNG.seed=2)
    }
  }
  if (is.null(params$a)) {
    modelstring=model2 # uses lognormal prior on r
    if (inf.type==1) { #  ... conditional inference with a beta prior on r
      jdata = list(y=y, pyr = pyr, k=k, a=params$mu, b=params$sd, d0=bg.events, d1=bg.pyr, del0=bg.events/bg.pyr)
      inits1 <- list(r=1,  .RNG.name="base::Super-Duper", .RNG.seed=1)
      inits2 <- list(r=0.95, .RNG.name="base::Wichmann-Hill", .RNG.seed=2)
    }
    if (inf.type==2) { # ... unconditional inference with a log-normal prior on r and gamma prior on del0
      jdata = list(y=y, pyr = pyr,  k=k, a=params$mu, b=params$sd, d0=bg.events, d1=bg.pyr)
      inits1 <- list(r=1, del0=bg.events/bg.pyr,  .RNG.name="base::Super-Duper", .RNG.seed=1)
      inits2 <- list(r=0.95, del0=bg.events/bg.pyr, .RNG.name="base::Wichmann-Hill", .RNG.seed=2)
    }
  }

  tmpModelFile=tempfile()
  tmps=file(tmpModelFile,"w")
  cat(modelstring,file=tmps)
  close(tmps)

  model=jags.model(tmpModelFile, data=jdata, n.chains=mc.params$nc,  inits=list(inits1,inits2), quiet=TRUE)
  update(model,n.iter=mc.params$burn, progress.bar = "none")
  sink ("dump.txt", append=T)
  output=coda.samples(model=model,variable.names=c("p", "r", "del0" ), n.iter=mc.params$iter, thin=1)
  sink()
  #invisible(file.remove("dump.txt"))

  para.names = labels((output[[1]]))[[2]] # getting the parameter names

  para=c()
  for (i in 1: length(para.names)){
    tmp = c(); for (j in 1:mc.params$nc) tmp = c(tmp, output[[j]][,i])
    para= cbind(para, tmp)
  }
  para = data.frame(para)
  names(para) = para.names
  return(para)
  }

#--- only used internally within bdribs package calls
drplot= function(x, ...)
{
  z = list(...)
  clr = "black";
  if (!is.null(z$clr)) clr = z$clr
  {
    op= par(mar=c(5.1, 2.5,3.1, .5))
    if (!is.null(z$xlm)) xlm =z$xlm
    else xlm = range(x)
    if(is.null(z$xlb)) z$xlb=""
    if(is.null(z$ttl)) z$ttl=""
    if (!is.null(z$ygrd)) ylm = range(z$ygrd)
    else ylm = NULL
    ttl0 = "Bayesian Detection of Risk using Inference on Blinded Safety data (BDRIBS)"
    plot(density(x), type="n", col=clr,lwd =2, xlim=xlm,ylim=ylm, axes=F, xlab="", ylab="",main="")
    points(z$rlst, z$frl, lwd=2, col="gray", type="l")
    points(density(x), type="l", col=clr,lwd =2)
    axis(1, at = z$xgrd, cex.axis=0.7,tck=-0.02,mgp=c(3, .5, 0))
    axis(2, at = z$ygrd, cex.axis=0.7,las=2, tck=-0.02,mgp=c(3, .5, 0))
    abline(h=z$ygrd, v=z$xgrd, col="gray", lty=3)
    title(xlab="r", cex.lab=0.8,line=1.5)
    title (main= z$ttl, cex.main=0.8,line=0.5)
    title (sub=z$ftn, col.sub="gray", cex.sub=0.7, line=3)
    title (sub=ttl0, col.sub="gray", cex.sub=0.7, line=4)
    box()
    par(op)
  }
}

