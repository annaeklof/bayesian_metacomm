
require(tidyverse)
require(RColorBrewer)
theme_set(theme_bw())
dpalette <- brewer.pal(9, "Blues")[c(4, 5, 6, 7, 8)]


expand.grid(piR=c(seq(0, 0.99, l=5001), seq(0.99, 1, l=15001)),
                  lR=c(0.5, 1, 2, 3, 5)) %>%
    as_tibble %>%
    mutate(lMR=-lR*log(1-piR)) %>%
    filter(lMR>=0, lMR<=5) %>%
    ggplot +
    aes(x=piR, y=lMR, colour=factor(lR), size=factor(lR), linetype=factor(lR)) +
    geom_path() +
    scale_x_continuous(name=expression(paste(pi[R])), labels=abbreviate) +
    scale_y_continuous(name=expression(paste(lambda[M[R]]))) +
    scale_colour_manual(name=expression(paste(lambda[R])), values=dpalette) +
    scale_size_manual(name=expression(paste(lambda[R])),
                      values=c(0.5, 0.6, 0.5, 0.5, 0.5)) +
    scale_linetype_manual(name=expression(paste(lambda[R])),
                          values=c(2, 1, 2, 2, 2))
## ggsave("../SI/R_persist.pdf", width=3.2, height=2.4)

expand.grid(s=seq(0, 1, l=20002), lMC=c(0.5, 1, 2, 3, 4)) %>%
    as_tibble %>%
    mutate(lR=s/(s-exp(-lMC))) %>%
    ggplot +
    aes(x=s, y=lR, colour=factor(lMC)) +
    geom_path(size=0.5) +
    scale_x_continuous(name=expression(paste((1-pi[C])(1-pi[R]))),
                       labels=abbreviate) +
    scale_y_continuous(name=expression(paste(lambda[R])), limits=c(1, 5)) +
    scale_colour_manual(values=dpalette, name=expression(paste(lambda[M[C]])))
## ggsave("../SI/C_persist.pdf", width=3.2, height=2.4)


res <- 301

chain <- expand.grid(lambdaM=seq(0, 8, l=res), pi0=seq(0, 1, l=res)) %>%
    as_tibble %>%
    filter(pi0!=0) %>%
    mutate(lambda=lambdaM/log(1/(1-pi0)), TL=0)
for (i in 1:nrow(chain)) {
    if (chain$lambda[i]<1) next
    while ((chain$lambda[i]>1)&(chain$TL[i]<8)) {
        chain$lambda[i] <- chain$lambdaM[i]/(chain$lambdaM[i]/chain$lambda[i]-log(
        (1-chain$pi0[i])*(1-1/chain$lambda[i])))
        chain$TL[i] <- chain$TL[i] + 1
    }
}

chain2 <- expand.grid(lambdaM=seq(0, 8, l=res), pi0=seq(0, 1, l=res)) %>%
    as_tibble %>%
    filter(pi0!=0) %>%
    mutate(TL=round(lambdaM*log(1/pi0))) %>%
    mutate(TL=ifelse(TL>8, 8, TL))

masterchain <- left_join(chain, chain2, by=c("lambdaM", "pi0")) %>%
    gather(model, TL, c("TL.x", "TL.y")) %>%
    mutate(model=ifelse(model=="TL.x", "iterated", "approximate")) %>%
    mutate(model=factor(model, ordered=TRUE, levels=c("iterated", "approximate")),
           TL=as.character(TL)) %>%
    mutate(TL=ifelse(TL=="8", ">7", TL))

masterchain %>%
    ggplot +
    geom_raster(aes(x=pi0, y=lambdaM, fill=TL), interpolate=TRUE) +
    scale_x_continuous(name=expression(paste(pi)), expand=c(0,0)) +
    scale_y_continuous(name=expression(paste(lambda[M])), expand=c(0,0)) +
    scale_colour_brewer(palette="PuBu") +
    scale_fill_brewer(palette="PuBu", name="Trophic levels") +
    facet_wrap(~model) +
    theme(panel.grid=element_blank())

## ggsave("../SI/TL_chain_approx.pdf", width=6, height=2.7)

