#' @title sensitivity plot - plot of range of posterior probability corresponding to a range of background rate
#'
#' @description plot of range of posterior probability corresponding to a range of background rate
#'
#' @param Y range on number of AESI events for which sensitivity range to be drawn (default = 5:9)
#' @param pyr total patient-year exposure where  AESI events occurred (default =800)
#' @param k allocation ratio (T:C) (default =1)
#' @param tol clinically meaningful relative risk (default =1.2)
#' @param bg.evnt background (historical) number of events in the control group (default =18)
#' @param bg.pyr background (historical) patient-year exposure in the control group (default =4000)
#' @param bg.ci.coef range of background rate estimate to be obtained from bg.ci.coef*100\% CI (default =0.9); takes any value between 0.5 and 0.999.
#' @param bg.rng range of background rate - if specified then bg.evnt, bg.pyr, and bg.ci.coef will be ignored (default =NULL)
#' @param add.mid indicator variable to plot P(r>tol | Y, pyr) under inf.type=2 - requires related parameters to be supplied (default =F)
#' @param ... to supply remaining parameters of bdribs call (other than y, pyr, k, bg.events, bg.pyr) for bdribs call when supplied will override the default values
#' @return returns a plot of P(r>tol | Y, pyr) over the range of background rate
#' @examples ## Sample calls
#'      #run 1: The sensitivity plot
#'      bdribs.sensitivity(Y=5:9,pyr=800,k=1, tol=1.2, bg.evnt=18, bg.pyr=4000,bg.ci.coef=0.90)
#'      #run 2: The sensitivity plot
#'      bdribs.sensitivity(Y=5:9,pyr=800,k=1, tol=1.2, bg.evnt=18, bg.pyr=4000,bg.ci.coef=0.90,
#'      add.mid=TRUE)
#'      #run 3: Using bg.rng parameter
#'      bdribs.sensitivity(Y=5:9,pyr=800,k=1, tol=1.2, bg.rng = c(0.0030, 0.0045, 0.0065))
#' @export
bdribs.sensitivity <- function(Y=5:9, pyr=800, k=1, tol=1.2, bg.evnt=18,
                               bg.pyr=4000,bg.ci.coef=0.90, bg.rng=NULL, add.mid=FALSE, ...)
{
  z = list(...) # setting up defaults values when those are not supplied
  if (is.null(z$p.params)) p.params=list (a=1, b=1)
  else p.params=z$p.params
  if (is.null(z$r.params)) r.params=list (mu=0, sd=2)
  else r.params=z$r.params
  if (is.null(z$mc.params)) mc.params=list(burn=1000, iter=10000, nc=2)
  else mc.params=z$mc.params
  if (is.null(z$adj.k)) adj.k=FALSE
  else adj.k=z$adj.k


  if(tol<1){
    if (tol <0) stop ("tolerance parameter 'tol' must be greater than 0 and typically should be 1 or more")
    else warning("tolerance parameter 'tol' typically should be 1 or more")
  }
  val=c(); lo=c(); up= c(); mid = c();
  if(is.null(bg.rng)) {
    if (bg.evnt <=0 | bg.pyr <=0) stop("Values of bg.evnt or bg.pyr must be both positive")
    if (bg.ci.coef <0.5 | bg.ci.coef > 0.999) stop("Value of bg.ci.coef is out of range - see more in Details")
    alp = (1-bg.ci.coef)/2
    bg.rng = c(qgamma(alp, bg.evnt+0.5, bg.pyr+0.5), bg.evnt/bg.pyr, qgamma(1-alp, bg.evnt+0.5, bg.pyr+0.5))
  }
  else add.mid=F # note when bg.rng is provided then we cannot run inf.type=2
  for (j in Y){
    ans= bdribs(y=j,pyr=pyr,k=k, bg.rate=bg.rng[2],plots=F, prnt=F)
    val = c(val, mean(ans$r>tol))

    ans= bdribs(y=j,pyr=pyr,k=k, bg.rate=bg.rng[1],plots=F, prnt=F,
                p.params = p.params,r.params = r.params, mc.params=mc.params, adj.k=adj.k)
    up = c(up, mean(ans$r>tol))

    ans= bdribs(y=j,pyr=pyr,k=k, bg.rate=bg.rng[3],plots=F, prnt=F,
                p.params = p.params,r.params = r.params, mc.params=mc.params, adj.k=adj.k)
    lo = c(lo, mean(ans$r>tol))

    if(add.mid)
    {
      ans= bdribs(y=j,pyr=pyr,k=k, bg.events=bg.evnt, bg.pyr=bg.pyr, inf.type=2,plots=F, prnt=F,
                  p.params = p.params,r.params = r.params, mc.params=mc.params, adj.k=adj.k)
      mid = c(mid, mean(ans$r>tol))
    }
  }

  df = data.frame(tol = tol, events=Y, pyr = pyr, est=val, lo, up)
  if(add.mid)  df = data.frame(df, mid)

  sm.error.bar <- function(x, upper, lower, wd=0.05,...){
    if(length(x) !=length(lower) | length(lower) != length(upper))
      stop("dimension error")
    arrows(x,upper, x, lower, angle=90, code=3, lwd=2, col="gray" , length=wd, ...)
  }

  xlb = paste( "Y (Total events at E =", df$pyr[1], "patient-years)")
  ylb = paste("P(r >", df$tol[1], " | Y, E =",df$pyr[1], ")\n" )
  ttl = "Sensitivity analyses to account for uncertainties in the true background rate"

  op= par(mar=c(3.1, 3.5,1.7, .5))
  plot(df$events, df$est, type="l", col="cyan", lwd =2,
       ylim=c(0, ceiling(max(df$up))), axes=F, xlab="", ylab="",main=ttl, cex.main=.8)
  title( ylab = ylb, xlab=xlb, cex.lab=0.8,line=1.5)
  sm.error.bar(df$events,df$up, df$lo)
  points(df$events, df$est, lwd=1, pch=21, bg="white", cex=1.5, type="p")
  axis(1, at = df$events, cex.axis=0.9,col.axis="gray2",mgp=c(3, .5, 0))
  axis(2, at = seq(0, ceiling(max(df$up)),.1),col.axis="gray2", cex.axis=0.9,las=2, tck=-0.02,mgp=c(3, .5, 0))

  abline(h=seq(0, ceiling(max(df$up)), 0.1), v=df$events, col="lightgray", lty=3)
  if(add.mid) points(df$events, df$mid, pch=17, type="p", col="blue", cex=1.5)
  box()

  invisible(df)
}
