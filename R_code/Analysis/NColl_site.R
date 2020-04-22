###Sites model:

#Sites model:
#collagen
yfix <- diag(19) * 1e+03
sites_prior<- list(B = list(mu = rep(0, 19),V = yfix), 
                   R = list(V = 1, nu = 0.002))
#carbonate
yfix_b <- diag(13) * 1e+03
sites_prior2<- list(B = list(mu = rep(0, 13),V = yfix_b), 
                    R = list(V = 1, nu = 0.002))
##Subset data
Nitrogen<-Nit_DF[,c(2,5,7,10)]
Nitrogen<-Nitrogen[!Nitrogen$name=="Josipovo",]
Nitrogen$name<-as.character(Nitrogen$name)
Nitrogen$name<-as.factor(Nitrogen$name)
CarbCol<-Nit_DF[,c(2,5,7,9)]
CarbCol<-CarbCol[!CarbCol$name=="Josipovo",]
CarbCol$name<-as.character(CarbCol$name)
CarbCol$name<-as.factor(CarbCol$name)
CCarb<-Nit_DF[,c(2,5,7,11)]
CCarb<-CCarb[!CCarb$name=="Josipovo",]
CCarb<-CCarb[!is.na(CCarb$d13Ccarb),]
CCarb$name<-as.character(CCarb$name)
CCarb$name<-as.factor(CCarb$name)

##Run sites models in triplicate using parallel processing:
#DIC - Deviance model selection
##Site*Sex - 556
##Site*Sex + TimePeriod - 556
##Site + Sex + TimePeriod - 546
##Site + Sex + (1|TimePeriod)
##Site + Sex - 546
##Site + TimePeriod - 548
##Site - 548
##Site + (1|TimePeriod + Sex) - 546
##Site + (Sex*TimePeriod) - 870
##Site + (1|TimePeriod) - 548
#Null - 870

#collagen
yfix <- diag(1) * 1e+03
sites_prior<- list(B = list(mu = rep(0, 1),V = yfix), 
                   R = list(V = 1, nu = 0.002))
tmp<-MCMCglmm(d15Ncoll~name,
              family="gaussian",
              data=Nitrogen,
              nitt=nitt,
              burnin=burnin,
              thin=thin,
              prior=sites_prior)

yfix <- diag(9) * 1e+03
sites_prior_b<- list(B = list(mu = rep(0, 5),V = yfix), 
                   R = list(V = 1, nu = 0.002))
prior_final_sites<-list(B = list(mu= diag(10)*0, V=diag(10)*1e+03),
                  R = list(nu=0.002, V=1),
                  G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

prior_final_sites_b<-list(B = list(mu= diag(18)*0, V=diag(18)*1e+03),
                        R = list(nu=0.002, V=1),
                        G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

tmp2<-MCMCglmm(d15Ncoll~name*Sex,
               ~ idh(TimePeriod),
              rcov=~ units, 
              family="gaussian",
              data=Nitrogen,
              nitt=nitt,
              burnin=burnin,
              thin=thin,
              prior=prior_final_sites_b)
plot.estimates(tmp2)

Sites_Nitrogen<- mclapply(1:3, function(i) {
  MCMCglmm(d15Ncoll~name*Sex + TimePeriod,
           family="gaussian",
           data=Nitrogen,
           nitt=nitt,
           burnin=burnin,
           thin=thin,
           prior=sites_prior)
}, mc.cores=3)

Sites_CarbCol<- mclapply(1:3, function(i) {
  MCMCglmm(d13Ccoll~name*Sex + TimePeriod,
           family="gaussian",
           data=CarbCol,
           nitt=nitt,
           burnin=burnin,
           thin=thin,
           prior=sites_prior)
}, mc.cores=3)

Sites_CCarb<- mclapply(1:3, function(i) {
  MCMCglmm(d13Ccarb~name*Sex + TimePeriod,
           family="gaussian",
           data=CCarb,
           nitt=nitt,
           burnin=burnin,
           thin=thin,
           prior=sites_prior2)
}, mc.cores=3)

N<-MCMCglmm(d15Ncoll~name*Sex + TimePeriod,
            family="gaussian",
            data=Nitrogen,
            nitt=nitt,
            burnin=burnin,
            thin=thin,
            prior=sites_prior)
