
require(deSolve)
require(tidyverse)

cpal <- c("#999999","#E69F00","#56B4E9","#009E73","#0072B2","#CC79A7","#D55E00")

## eq <- function(T, L, pars) { return(list(-log((1-pars$p)*(1-L))))}
## dat <- tibble(p=c(10^seq(-4.5,-1,by=0.05), seq(0.11,0.999,by=0.01)), maxTL=0)
## tls <- seq(0, 10, by=0.01)
## for (i in 1:nrow(dat)[1]) {
##     params <- list(p=dat$p[i])
##     dat$maxTL[i] <- ode(func=eq, y=0, times=tls, parms=params) %>%
##         as.data.frame %>%
##         as_tibble %>%
##         rename(tl=time, L=`1`) %>%
##         filter(!is.nan(L)) %>%
##         filter(L==max(L)) %>%
##         pull(tl)
## }

dat <- tibble(p=c(10^seq(-4.3, -1, by=0.01),
                  seq(0.11, 0.99999, by=0.0001))) %>%
    mutate(exact=0, approximate=log(1/p))
for (i in 1:nrow(dat)) {
    dat$exact[i] <- (integrate(function(r) -1/log(r), 0, 1-dat$p[i],
                               rel.tol=1e-12, abs.tol=1e-12,
                               stop.on.error=FALSE)$value)/(1-dat$p[i])
}

dat %>%
    gather("type", "TL", c("exact", "approximate")) %>%
    ggplot +
    geom_line(aes(x=p, y=TL, colour=type), alpha=0.8) +
    scale_x_continuous(name=expression(paste(pi))) +
    scale_y_continuous(name=expression(paste(tau))) +
    scale_colour_manual(values=c(cpal[2], cpal[5])) +
    theme_bw()

## ggsave("../SI/TL_diffeq.pdf", width=4.2, height=2.7)

