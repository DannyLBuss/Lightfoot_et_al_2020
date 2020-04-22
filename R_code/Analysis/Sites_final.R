##Sites - nitrogen

prior_final_sites_b<-list(B = list(mu= diag(18)*0, V=diag(18)*1e+03),
                          R = list(nu=0.002, V=1),
                          G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

SiteM1<-MCMCglmm(d15Ncoll~name*Sex,
               ~ idh(TimePeriod),
               rcov=~ units, 
               family="gaussian",
               data=Nitrogen,
               nitt=nitt,
               burnin=burnin,
               thin=thin,
               prior=prior_final_sites_b)

prior_final_sites_c<-list(B = list(mu= diag(10)*0, V=diag(10)*1e+03),
                          R = list(nu=0.002, V=1),
                          G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

SiteM2<-MCMCglmm(d15Ncoll~name + Sex,
               ~ idh(TimePeriod),
               rcov=~ units, 
               family="gaussian",
               data=Nitrogen,
               nitt=nitt,
               burnin=burnin,
               thin=thin,
               prior=prior_final_sites_c)

prior_final_sites_d<-list(B = list(mu= diag(9)*0, V=diag(9)*1e+03),
                          R = list(nu=0.002, V=1),
                          G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

SiteM3<-MCMCglmm(d15Ncoll~name,
                 ~ idh(TimePeriod),
                 rcov=~ units, 
                 family="gaussian",
                 data=Nitrogen,
                 nitt=nitt,
                 burnin=burnin,
                 thin=thin,
                 prior=prior_final_sites_d)

prior_final_sites_e<-list(B = list(mu= diag(2)*0, V=diag(2)*1e+03),
                          R = list(nu=0.002, V=1),
                          G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

SiteM4<-MCMCglmm(d15Ncoll~Sex,
                 ~ idh(TimePeriod),
                 rcov=~ units, 
                 family="gaussian",
                 data=Nitrogen,
                 nitt=nitt,
                 burnin=burnin,
                 thin=thin,
                 prior=prior_final_sites_e)

prior_final_sites_f<-list(B = list(mu= diag(1)*0, V=diag(1)*1e+03),
                          R = list(nu=0.002, V=1),
                          G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

SiteM5<-MCMCglmm(d15Ncoll~1,
                 ~ idh(TimePeriod),
                 rcov=~ units, 
                 family="gaussian",
                 data=Nitrogen,
                 nitt=nitt,
                 burnin=burnin,
                 thin=thin,
                 prior=prior_final_sites_f)

##Sites - carbon collagen

prior_final_sites_b<-list(B = list(mu= diag(18)*0, V=diag(18)*1e+03),
                          R = list(nu=0.002, V=1),
                          G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

SiteM1_CC<-MCMCglmm(d13Ccoll~name*Sex,
                 ~ idh(TimePeriod),
                 rcov=~ units, 
                 family="gaussian",
                 data=CarbCol,
                 nitt=nitt,
                 burnin=burnin,
                 thin=thin,
                 prior=prior_final_sites_b)

prior_final_sites_c<-list(B = list(mu= diag(10)*0, V=diag(10)*1e+03),
                          R = list(nu=0.002, V=1),
                          G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

SiteM2_CC<-MCMCglmm(d13Ccoll~name + Sex,
                 ~ idh(TimePeriod),
                 rcov=~ units, 
                 family="gaussian",
                 data=CarbCol,
                 nitt=nitt,
                 burnin=burnin,
                 thin=thin,
                 prior=prior_final_sites_c)

prior_final_sites_d<-list(B = list(mu= diag(9)*0, V=diag(9)*1e+03),
                          R = list(nu=0.002, V=1),
                          G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

SiteM3_CC<-MCMCglmm(d13Ccoll~name,
                    ~ idh(TimePeriod),
                    rcov=~ units, 
                    family="gaussian",
                    data=CarbCol,
                    nitt=nitt,
                    burnin=burnin,
                    thin=thin,
                    prior=prior_final_sites_d)

prior_final_sites_e<-list(B = list(mu= diag(2)*0, V=diag(2)*1e+03),
                          R = list(nu=0.002, V=1),
                          G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

SiteM4_CC<-MCMCglmm(d13Ccoll~Sex,
                 ~ idh(TimePeriod),
                 rcov=~ units, 
                 family="gaussian",
                 data=CarbCol,
                 nitt=nitt,
                 burnin=burnin,
                 thin=thin,
                 prior=prior_final_sites_e)

prior_final_sites_f<-list(B = list(mu= diag(1)*0, V=diag(1)*1e+03),
                          R = list(nu=0.002, V=1),
                          G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

SiteM5_CC<-MCMCglmm(d13Ccoll~1,
                 ~ idh(TimePeriod),
                 rcov=~ units, 
                 family="gaussian",
                 data=CarbCol,
                 nitt=nitt,
                 burnin=burnin,
                 thin=thin,
                 prior=prior_final_sites_f)

###Plot results
SiteM3_CC<-mclapply(1:3, function(i){
  MCMCglmm(d13Ccoll~name,
                    ~ idh(TimePeriod),
                    rcov=~ units, 
                    family="gaussian",
                    data=CarbCol,
                    nitt=nitt,
                    burnin=burnin,
                    thin=thin,
                    prior=prior_final_sites_d)
}, mc.cores=3)
Sites_CarbCol <- lapply(SiteM3_CC, function(m) m$Sol)
Sites_CarbCol <- do.call(mcmc.list, Sites_CarbCol)
plot.estimates(Sites_CarbCol)

SiteM2_N<-mclapply(1:3, function(i){
  MCMCglmm(d15Ncoll~name + Sex,
                 ~ idh(TimePeriod),
                 rcov=~ units, 
                 family="gaussian",
                 data=Nitrogen,
                 nitt=nitt,
                 burnin=burnin,
                 thin=thin,
                 prior=prior_final_sites_c)
}, mc.cores=3)
SitesM2_Na <- lapply(SiteM2_N, function(m) m$Sol)
SitesM2_Na  <- do.call(mcmc.list, SitesM2_Na)
plot.estimates(SitesM2_Na)

N<-MCMCglmm(d15Ncoll~(name + Sex)-1,
              ~ idh(TimePeriod),
              rcov=~ units, 
              family="gaussian",
              data=Nitrogen,
              nitt=nitt,
              burnin=burnin,
              thin=thin,
              prior=prior_final_sites_c)
summary(N)
plot(N)
Ccol<-MCMCglmm(d13Ccoll~(name + Sex)-1,
                 ~ idh(TimePeriod),
                 rcov=~ units, 
                 family="gaussian",
                 data=CarbCol,
                 nitt=nitt,
                 burnin=burnin,
                 thin=thin,
                 prior=prior_final_sites_c)
summary(Ccol)
plot(Ccol)

prior_final_sites_f<-list(B = list(mu= diag(6)*0, V=diag(6)*1e+03),
                          R = list(nu=0.002, V=1),
                          G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

Ccarb<-MCMCglmm(d13Ccarb~name,
                 rcov=~ units, 
                 family="gaussian",
                 data=CarC_DF,
                 nitt=nitt,
                 burnin=burnin,
                 thin=thin,
                 prior=prior_final_sites_f)

prior_final_sites_f<-list(B = list(mu= diag(6)*0, V=diag(6)*1e+03),
                          R = list(nu=0.002, V=diag(6)),
                          G = list(G1 = list(V = diag(11)*1, nu = 0.02, alpha.mu = rep(0,11), alpha.V = diag(11)*1000),
                                   G2= list(V = diag(11)*1, nu = 0.02, alpha.mu = rep(0,11), alpha.V = diag(11)*1000)))
Ccarb<-MCMCglmm(d13Ccarb~(Sex/TimePeriod) + (name/TimePeriod),
                  rcov=~idh(name):units, 
                  family="gaussian",
                  data=CarC_DF,
                  nitt=nitt,
                  burnin=burnin,
                  thin=thin,
                  prior=prior_final_sites_f)
Ccarb<-MCMCglmm(d13Ccarb~1,
                rcov=~idh(name):units, 
                family="gaussian",
                data=CarC_DF,
                nitt=nitt,
                burnin=burnin,
                thin=thin,
                prior=prior_final_sites_f)
summary(Ccarb)
plot(Ccarb)

###Final_sites_models
Sites_CCarb<- mclapply(1:3, function(i) {
  MCMCglmm(d13Ccarb~(Sex/TimePeriod) + (name/TimePeriod),
                rcov=~idh(name):units, 
                family="gaussian",
                data=CarC_DF,
                nitt=nitt,
                burnin=burnin,
                thin=thin)
  }, mc.cores=3)

Sites_Nit<- mclapply(1:3, function(i) {
  MCMCglmm(d15Ncoll~(Sex/TimePeriod) + (name/TimePeriod),
           rcov=~idh(name):units, 
           family="gaussian",
           data=Nitrogen,
           nitt=nitt,
           burnin=burnin,
           thin=thin)
}, mc.cores=3)

Sites_CColl<- mclapply(1:3, function(i) {
  MCMCglmm(d13Ccoll~(Sex/TimePeriod) + (name/TimePeriod),
           rcov=~idh(name):units, 
           family="gaussian",
           data=CarbCol,
           nitt=nitt,
           burnin=burnin,
           thin=thin)
}, mc.cores=3)

Sites_Nit <- lapply(Sites_Nit, function(m) m$Sol)
Sites_Nit <- do.call(mcmc.list, Sites_Nit)
Sites_CColl <- lapply(Sites_CColl, function(m) m$Sol)
Sites_CColl <- do.call(mcmc.list, Sites_CColl)
Sites_CCarb <- lapply(Sites_CCarb, function(m) m$Sol)
Sites_CCarb <- do.call(mcmc.list, Sites_CCarb)

par(mfrow=c(1,3))
plot.estimates(Sites_Nit)
plot.estimates(Sites_CColl)
plot.estimates(Sites_CCarb)

###Final_sites_models
Sites_CCarb<-MCMCglmm(d13Ccarb~(Sex/TimePeriod) + (name/TimePeriod)-1,
           rcov=~idh(name):units, 
           family="gaussian",
           data=CarC_DF,
           nitt=nitt,
           burnin=burnin,
           thin=thin)

Sites_Nit<-MCMCglmm(d15Ncoll~(Sex/TimePeriod) + (name/TimePeriod)-1,
           rcov=~idh(name):units, 
           family="gaussian",
           data=Nitrogen,
           nitt=nitt,
           burnin=burnin,
           thin=thin)

Sites_CColl<-MCMCglmm(d13Ccoll~(Sex/TimePeriod) + (name/TimePeriod)-1,
           rcov=~idh(name):units, 
           family="gaussian",
           data=CarbCol,
           nitt=nitt,
           burnin=burnin,
           thin=thin)


Sites_CCarb<-MCMCglmm(d13Ccarb~(Sex/name)-1,
                      rcov=~idh(name):units, 
                      family="gaussian",
                      data=CarC_DF,
                      nitt=nitt,
                      burnin=burnin,
                      thin=thin)
summary(Sites_CCarb)


Sites_Nit<-MCMCglmm(d15Ncoll~(Sex/TimePeriod)-1,
                    random =~name,
                    rcov=~idh(name):units, 
                    family="gaussian",
                    data=Nitrogen,
                    nitt=nitt,
                    burnin=burnin,
                    thin=thin)
summary(Sites_Nit)
