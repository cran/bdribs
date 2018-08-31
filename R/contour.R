#' @title contour plot - draws plot (optional) and returns a matrix/grid of posterior values
#'
#' @description Contour plot of posterior probabilities on a range of (y, E) values
#'
#' @param ymax maximum number of AESI event for which contour plot to be drawn
#' @param pyrmax maximum risk exposure (in patient-year)
#' @param eincr increment of patient-year exposures (default = 50)
#' @param tol the maximum tolerance value of relative risk r (default =1)
#' @param k allocation ratio (T:C)
#' @param bg.rate estimated background rate (historical control rate) per patient-year (using inf.type=1)
#' @param ... to supply remaining parameters for bdribs call when supplied will override the default values
#' @param plt whether a contour plot to be drawn (default = TRUE)
#' @return returns contour plot matrix over the grid specified
#' @examples ## Sample calls
#'      #run 1: The contour plot
#'      \donttest{
#'      bdribs.contour(ymax=15,pyrmax=2000,eincr=250,tol=1.5,k=2, bg.rate=0.0045)
#'      #run 2: Monitoring blinded AE over time using contour plot
#'      bdribs.contour(ymax=15,pyrmax=2000,eincr=250,tol=1.5,k=2, bg.rate=0.0045)
#'      obs.pyr=c(300,570,650,800, 1200, 1500)
#'      obs.y=c(2,4,5,6,10,12)
#'      points(obs.pyr, obs.y,type="p",pch=16, cex=1.4,col="maroon")
#'      if (length(obs.y)>1) points(c(0,obs.pyr), c(0,obs.y), type="s", lty=3, lwd=2,
#'      col="black")
#'      }
#' @importFrom  stats quantile density update dbeta dlnorm pbeta plnorm qgamma
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline axis box contour plot title arrows par points
#' @export
bdribs.contour<-function(ymax, pyrmax,eincr, tol, k, bg.rate,plt=TRUE, ...)
{
  zz = list(...) # setting up defaults values when those are not supplied
#  if (!exists("cntr.dta")) initiate.contr()
  cntr.dta=bdribs.contour.data(ymax, pyrmax,eincr, k, bg.rate, z=zz)

  y = 0:ymax
  pyr = seq(0, pyrmax,eincr)
  pyrd = pyr; pyrd[1]=0.1

  Mi = length(y)
  Mj =length(pyr)
  post = matrix(0, Mi, Mj)
  cnt=0

  for (i in 1: Mi) for (j in 1:Mj){
    post[i,j] = mean(cntr.dta$post[i,j,] > tol)
  }

  if (plt){
    op= par(mar=c(3.1, 3.5,1.7, .5))
    ttl=paste("Contour Chart for Bayesian Monitoring of Blinded AEs\n P[ r > ",
              tol,"| y,E,k=",k,", bg.rate=", bg.rate,"]",sep="")

    plot(range(pyr), range(y), type="n", xlab = "",main=ttl,
         ylab = "", cex.main=0.8,  axes=F)
    axis(1, at=pyr, cex.axis=0.7,mgp=c(3, .5, 0),tck=-0.02)
    axis(2, at=y, cex.axis=0.7, las=2,mgp=c(3, .5, 0),tck=-0.02)
    abline(v=pyr, h=y, col="lightgray",lty=3)
    title( ylab = "Total no of AEs (y)",
           xlab="Total risk exposure in patient-years (E)", cex.lab=0.8,line=1.5)
    box()
    cmap = colorRampPalette(c("green2","blue","red2")) (9)
    contour (pyr,y,t(post), add=T, col=cmap,levels=1:9/10)
    par(op)
  }
  return(invisible(post))
}

#initiate.contr = function()
#{
  #cat("initiation done")
#  assign("cntr.dta", NULL, .GlobalEnv)
  #lockBinding("cntr.dta", .GlobalEnv)
#  cntr.dta <<- list(post=NULL, params= c(0, 0,0, 0, 0))
#}

bdribs.contour.data<-function(ymax, pyrmax,eincr, k, bg.rate, z)
{
    #if (is.null(cntr.dta$post) |  prod(cntr.dta$params==c(ymax, pyrmax,eincr, k, bg.rate))==0) {
    #attributes(cntr.dta)<<- list(post=NULL, params=c(ymax, pyrmax,eincr, k, bg.rate))
    #prep.contour.data(ymax, pyrmax,eincr, k, bg.rate)
    y = 0:ymax
    pyr = seq(0, pyrmax,eincr)
    pyrd = pyr; pyrd[1]=0.1

    Mi = length(y)
    Mj =length(pyr)

    # z = list(...) # setting up defaults values when those are not supplied
    if (is.null(z$p.params)) p.params=list (a=1, b=1)
    else p.params=z$p.params
    if (is.null(z$r.params)) r.params=list (mu=0, sd=2)
    else r.params=z$r.params
    if (is.null(z$mc.params)) mc.params=list(burn=1000, iter=10000, nc=2)
    else mc.params=z$mc.params
    if (is.null(z$adj.k)) adj.k=FALSE
    else adj.k=z$adj.k
    #print(p.params)

    post = array(0, dim=c(Mi, Mj, 20000) )# default size of posterior sample is 20000
    cnt=0
    for (i in 1: Mi) for (j in 1:Mj){
      bd = bdribs(y[i], pyrd[j], bg.rate = bg.rate, k=k,
                  p.params = p.params,r.params = r.params, mc.params=mc.params, adj.k=adj.k,
                  plots=F, prnt=F)
      post[i,j,] = bd$r
      cnt = cnt+1
      if (cnt %% 10 ==0) cat(paste("... ", round(cnt/Mi/Mj*100,0), "% done\n",sep=""))
    }
#    cntr.dta$post<<- post
#    cntr.dta$params <<- c(ymax, pyrmax,eincr, k, bg.rate)
#      }
    cntr.dta=list(post=post)
  return(cntr.dta)
}
