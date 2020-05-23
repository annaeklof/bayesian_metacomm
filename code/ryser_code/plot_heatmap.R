
require(tidyverse)

dat <- read_csv("../../data/ryser_data/results/tl_dat.csv")

dat %>%
  ggplot +
  aes(x=log10(xi), y=N, fill=tl) +
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0, 0), name=expression(paste(log[10](xi)))) +
  scale_y_continuous(expand=c(0, 0), name="number of patches") +
  scale_fill_gradient2(name="max trophic level", midpoint=3.5) +
  facet_grid(beta~alpha, labeller=label_bquote(rows=beta==.(beta),
                                               cols=alpha==.(alpha))) +
  theme_bw()
##ggsave("../../figures/Ryser_reproduced.pdf", width=5.5, height=3.8)



getlambda <- function(lsc, xi) {
  M <- read_csv(lsc, col_types="dd") %>%
    as.matrix(col_types=cols()) %>%
    dist %>%
    as.matrix
  M <- exp(-M/xi)
  diag(M) <- 0
  rs <- rowSums(M)
  lambda <- sum(rs*rs)/sum(rs)
  return(lambda)
}

dat2 <- expand_grid(N=seq(16, 70, by=2),
            xi=10^(seq(-2, 0.5, l=11)),
            lsc=Sys.glob("../../data/ryser_data/landscapes/*.csv")) %>%
  rowwise %>%
  mutate(lM=getlambda(lsc, xi)) %>%
  ungroup %>%
  mutate(tl=lM)

dat2 %>%
  ggplot +
  aes(x=log10(xi), y=N, fill=tl) +
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0, 0), name=expression(paste(log[10](xi)))) +
  scale_y_continuous(expand=c(0, 0), name="number of patches") +
  scale_fill_gradient2(name="max trophic level", midpoint=25,
                       labels=c(expression(paste(10%*%log(pi^-1))),
                                expression(paste(30%*%log(pi^-1))),
                                expression(paste(50%*%log(pi^-1)))),
                       breaks=c(10, 30, 50)) +
  theme_bw()
##ggsave("../../figures/Ryser_semianalytic.pdf", width=5.5, height=3.8)

