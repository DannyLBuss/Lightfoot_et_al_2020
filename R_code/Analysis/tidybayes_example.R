library(emmeans)
library(tidybayes)
vignette("models", package = "emmeans")
a<-MCMCglmm(d15Ncoll~TimePeriod-1,
         rcov=~ units, 
         family="gaussian",
         data=paired,
         nitt=nitt,
         burnin=burnin,
         thin=thin,
         prior=ps_prior)

a %>%
  emmeans(~TimePeriod, data=paired) %>%
  contrast(method="pairwise") %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x=.value, y= contrast)) +
  stat_eyeh(.width=c(.95,.8)) + labs(x="")

