
library(ggplot2)

ggplot_fform <- function(alpha, beta)
{
  fform <- data.frame(
    x <- sort(runif(10000)),
    y = pbeta(x, alpha, beta)
  )
  theme_bw() 
  par(pty="s")
  ggplot(fform, aes(x, y)) +
    geom_hline(yintercept=0, lty=2, lwd=0.8, color="black", alpha = 0.5) +
    geom_hline(yintercept=1, lty=2, lwd=0.8, color="black", alpha = 0.5) +
    geom_line(lwd=1.1) +
    xlab("Fraction of resources \n extinct") +
    ylab("Prob. of extinction \n of resources") +
    # coord_cartesian(0:1, 0:1) + 
    scale_x_continuous(breaks=c(0, 0.5, 1),
                       labels = c("0","0.5","1")) +
    scale_y_continuous(breaks=c(0, 0.5, 1),
                       labels = c(expression(pi[i]),"0.5","1")) +
    theme(axis.text = element_text(color="black", size = 16), 
          axis.title = element_text(color="black", size=18))  
}

plot_functional_form <- function(figuredir, alpha, beta){
      ggsave(
        paste0(figuredir, "fform_",alpha,"_",beta,".png"),
        ggplot_fform(alpha, beta),
        width = 3.25,
        height = 3.25,
        dpi = 1200
      )
}


