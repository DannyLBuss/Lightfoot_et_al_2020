#post.samp<-mcmc chain, N<-number of simulations e.g 5000
sim.pred = function(N, post.samp){
  index <- sample(dim(post.samp)[1],N,replace = T) 
  samp <- post.samp[index,]
  y.pred <- matrix(0, nrow = N, ncol = length(data))
  for(i in 1:N){
    y.pred[i,] <-   rnorm(1,samp[i,1],samp[i,2])
  }
  return(y.pred)
}

post.samp<-test
N<-5000
t<-sim.pred(N,post.samp)


sim.pred = function(N, post.samp){
  index <- sample(dim(post.samp)[1],N,replace = T) 
  samp <- post.samp[index,]
  y.pred <- matrix(0, nrow = N, ncol = length(post.samp))
  for(i in 1:N){
    y.pred[i,] <- rnorm(1,samp[i,1],samp[i,2])
  }
  return(y.pred)
}

Paired_site_NColl<-mclapply(1:3, function(i){
  MCMCglmm(d15Ncoll~TimePeriod-1,
           rcov=~ units, 
           family="gaussian",
           data=paired,
           nitt=nitt,
           burnin=burnin,
           thin=thin,
           prior=ps_prior)
}, mc.cores=3)

c<-list(Paired_site_NColl[[1]]$Sol,Paired_site_NColl[[2]]$Sol,Paired_site_NColl[[3]]$Sol)
d<-list(Paired_site_NColl[[1]]$VCV,Paired_site_NColl[[2]]$VCV,Paired_site_NColl[[3]]$VCV)

t<-chainsPlot(c, densityplot=TRUE,legend.location = 'topleft', cex = 0.7, col=c("grey20","grey80","white"))
chainsPlot(d, densityplot=FALSE,legend.location = 'topleft', cex = 0.7)

samp[5,1]
samp[8,2]
