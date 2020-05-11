
plot_functional_forms <- function(figuredir){

  par(pty="s")
  frac <- sort(runif(10000))
  
  for(alpha in c(1,5)){
    for(beta in c(1,5)){
      png(paste0(figuredir,"fform_",alpha,"_",beta,".png"))
      f <- pbeta(frac, alpha, beta)
      # plot(f~frac,type='n',xlab="",ylab="",xaxt="n",yaxt="n",lwd=3)
      plot(f~frac,type='n',xlab="Fraction of resource extinct",
           ylab="Probability of \n extinction of resource",
           xaxt="n",yaxt="n",lwd=3,cex.lab=3)
      
      ticks <- seq(0, 1, by=0.5)
      axis(side=1, at=ticks, labels = FALSE)
      text(x=xtick,  par("usr")[3], 
           labels = xtick, pos = 1, xpd = TRUE, cex = 3)
      axis(side=2, at=ticks, labels = FALSE)
      text(par("usr")[1], ticks,  
           labels = ytick, pos = 2, xpd = TRUE, cex =2.3)
      
      lines(f~frac,col="grey21",lwd=8)
      abline(h=max(f), lty=3,lwd=6)
      abline(h=min(f), lty=3,lwd=6)
      dev.off()
    }
  }
}
