###R_Code for Lightfoot et al. 2019 MCMCglmm runs
#Set working directory and load in data files
setwd("~/Documents/Home_stuff/Arch/Lightfoot_et_al_stats_2019/")
Nit_DF<-read.csv("collagen_nitrogen.csv",header=T)
Car_DF<-read.csv("collagen_carbon.csv",header=T)
CarC_DF<-read.csv("carbonate_carbon.csv",header=T)

#Load packages
library(parallel)
library(coda)
library(MCMCglmm)
library(ggplot2)
library(lattice)
library(latticeExtra)
library(dplyr)
library(reshape2)

###Priors and MCMC settings:
##Set MCMC parameters
#no. of interations
nitt <- c(500000)
#length of burnin
burnin <- c(10000)
#amount of thinning
thin <- c(245)

##Priors
###Non=informative parameter expansion priors (often used in ecology and with non-normal data)

#Paired-sites model:                                          
yfixed <- diag(4) * 1e+03
ps_prior<- list(B = list(mu = rep(0, 4),V = yfixed), 
                R = list(V = 1, nu = 0.002))
#Max_model:
#TP*Sex = 2 intercepts, Alt = 1 slope, DS = 1 slope, DC = 1 slope
prior_final<-list(B = list(mu= diag(5)*0, V=diag(5)*1e+10),
                  R = list(nu=0.002, V=1),
                  G = list(G1 = list(V = diag(10)*1, nu = rep(0.02,10), alpha.mu = rep(0,10), alpha.V = diag(10)*1000)))

#Sites model:
#collagen
yfix <- diag(19) * 1e+03
sites_prior<- list(B = list(mu = rep(0, 19),V = yfix), 
                   R = list(V = 1, nu = 0.002))
#carbonate
yfix_b <- diag(13) * 1e+03
sites_prior2<- list(B = list(mu = rep(0, 13),V = yfix_b), 
                    R = list(V = 1, nu = 0.002))

###Paired-site model:
#subset data
paired<-Nit_DF[Nit_DF$name=="Koprivno",]
paired<-paired[,c(2,5,7,9:11)]

#run 3 chains of mcmc model in parallel for each response variable
Single_Nitrogen<-mclapply(1:3, function(i){
  MCMCglmm(d15Ncoll~TimePeriod*Sex-1,
           rcov=~ units, 
           family="gaussian",
           data=paired,
           nitt=nitt,
           burnin=burnin,
           thin=thin,
           prior=ps_prior)
}, mc.cores=3)

Single_CarbCol<-mclapply(1:3, function(i){
  MCMCglmm(d13Ccoll~TimePeriod*Sex-1,
           rcov=~ units, 
           family="gaussian",
           data=paired,
           nitt=nitt,
           burnin=burnin,
           thin=thin,
           prior=ps_prior)
}, mc.cores=3)

Single_CCarb<-mclapply(1:3, function(i){
  MCMCglmm(d13Ccarb~TimePeriod*Sex-1,
           rcov=~ units, 
           family="gaussian",
           data=paired,
           nitt=nitt,
           burnin=burnin,
           thin=thin,
           prior=ps_prior)
}, mc.cores=3)

#Group 3 chains back together as an mcmcglmm output (Sol - fixed effects)
Single_Nitrogen <- lapply(Single_Nitrogen, function(m) m$Sol)
Single_Nitrogen <- do.call(mcmc.list, Single_Nitrogen)
Single_CarbCol <- lapply(Single_CarbCol, function(m) m$Sol)
Single_CarbCol <- do.call(mcmc.list, Single_CarbCol)
Single_CCarb <- lapply(Single_CCarb, function(m) m$Sol)
Single_CCarb <- do.call(mcmc.list, Single_CCarb)

#Check model convergance - scale reduction factor should be ~ 1
convergence_checks<-list(gelman.diag(Single_Nitrogen),
                         gelman.diag(Single_CarbCol),
                         gelman.diag(Single_CCarb))
###Add some ifelse statement that prints to screen - models passed checkplots
convergence_checks

#Save_alltraceplots_to_file
jpg_names<-c(rep("Nitrogen_Collagen",3),rep("Carbon_Collagen",3),rep("Carbon_Carbonate",3))

#Check tracplots & -> change to - save traceplots to file
plot(Single_Nitrogen[1][1])
plot(Single_Nitrogen[2][1])
plot(Single_Nitrogen[3][1])
plot(Single_CarbCol[1][1])
plot(Single_CarbCol[2][1])
plot(Single_CarbCol[3][1])
plot(Single_CCarb[1][1])
plot(Single_CCarb[2][1])
plot(Single_CCarb[3][1])
dev.off()

#Check model convergance
plot.estimates(Single_Nitrogen)
plot.estimates(Single_CarbCol)
plot.estimates(Single_CCarb)

###Run MCMC model with 3 chains in parallel using 3 cores
Max_prior<-list(B = list(mu= diag(5)*0, V=diag(5)*1e+10),
                  R = list(nu=0.002, V=1),
                  G = list(G1 = list(V = diag(10)*1, nu = 0.02, alpha.mu = rep(0,10), alpha.V = diag(10)*1000)))

Max_Nitrogen<- mclapply(1:3, function(i) {
  MCMCglmm(d15Ncoll~TimePeriod + Alt + Sex + dist_from_Coast,
           ~ idh(name),
           rcov=~ units, 
           family="gaussian",
           data=Nit_DF,
           nitt=nitt,
           burnin=burnin,
           thin=thin,
           prior=Max_prior)
}, mc.cores=3)

Max_CarbCol<- mclapply(1:3, function(i) {
  MCMCglmm(d13Ccoll~TimePeriod + Alt + Sex + dist_from_Coast,
           ~ idh(name),
           rcov=~ units, 
           family="gaussian",
           data=Nit_DF,
           nitt=nitt,
           burnin=burnin,
           thin=thin,
           prior=Max_prior)
}, mc.cores=3)

Max_CCarb<- mclapply(1:3, function(i) {
  MCMCglmm(d13Ccarb~TimePeriod + Alt + Sex + dist_from_Coast,
           ~ idh(name),
           rcov=~ units, 
           family="gaussian",
           data=Nit_DF,
           nitt=nitt,
           burnin=burnin,
           thin=thin,
           prior=Max_prior)
}, mc.cores=3)

Max_Nitrogen <- lapply(Max_Nitrogen, function(m) m$Sol)
Max_Nitrogen <- do.call(mcmc.list, Max_Nitrogen)
Max_CarbCol <- lapply(Max_CarbCol, function(m) m$Sol)
Max_CarbCol <- do.call(mcmc.list, Max_CarbCol)
Max_CCarb <- lapply(Max_CCarb, function(m) m$Sol)
Max_CCarb <- do.call(mcmc.list, Max_CCarb)

#Check model convergance - create table of all 3 
convergence_checks_max<-list(gelman.diag(Max_Nitrogen),
                         gelman.diag(Max_CarbCol),
                         gelman.diag(Max_CCarb))

#Check tracplots
plot(Max_Nitrogen[1][1])
plot(Max_Nitrogen[2][1])
plot(Max_Nitrogen[3][1])
plot(Max_CarbCol[1][1])
plot(Max_CarbCol[2][1])
plot(Max_CarbCol[3][1])
plot(Max_CCarb[1][1])
plot(Max_CCarb[2][1])
plot(Max_CCarb[3][1])

#Check model convergance
plot.estimates(Max_Nitrogen)
plot.estimates(Max_CarbCol)
plot.estimates(Max_CCarb)

###Sites model:
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
Sites_Nitrogen<- mclapply(1:3, function(i) {
  MCMCglmm(d15Ncoll~name*Sex + TimePeriod,
                         family="gaussian",
                         data=Nitrogen,
                         nitt=Nit_DF,
                         burnin=burnin,
                         thin=thin,
           prior=sites_prior)
  }, mc.cores=3)

Sites_CarbCol<- mclapply(1:3, function(i) {
  MCMCglmm(d13Ccoll~name*Sex + TimePeriod,
                        family="gaussian",
                        data=Car_DF,
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

###Combine triplicate runs for each isotope model:
Sites_Nitrogen <- lapply(Sites_Nitrogen, function(m) m$Sol)
Sites_Nitrogen <- do.call(mcmc.list, Sites_Nitrogen)
Sites_CarbCol <- lapply(Sites_CarbCol, function(m) m$Sol)
Sites_CarbCol <- do.call(mcmc.list, Sites_CarbCol)
Sites_CCarb <- lapply(Sites_CCarb, function(m) m$Sol)
Sites_CCarb <- do.call(mcmc.list, Sites_CCarb)

#Check model convergance - create table of all 3 
convergence_checks_Sites<-list(gelman.diag(Sites_Nitrogen),
                         gelman.diag(Sites_CarbCol),
                         gelman.diag(Sites_CCarb))
#Check tracplots
plot(Sites_Nitrogen[1][1])
plot(Sites_Nitrogen[2][1])
plot(Sites_Nitrogen[3][1])
plot(Sites_CarbCol[1][1])
plot(Sites_CarbCol[2][1])
plot(Sites_CarbCol[3][1])
plot(Sites_CCarb[1][1])
plot(Sites_CCarb[2][1])
plot(Sites_CCarb[3][1])

#Check model convergance
par(mfrow=c(1,3))
plot.estimates(Sites_Nitrogen)
plot.estimates(Sites_CarbCol)
plot.estimates(Sites_CCarb)

###Variance models
library(lme4)
library(nlme)
library(dplyr)
library(ggplot2)
a3<-Nit_DF %>%
  group_by(name, Sex, TimePeriod) %>%
  summarise(NColl_Var=var(d15Ncoll),CColl_Var=var(d13Ccoll))
a2<-CarC_DF %>%
  group_by(name, Sex, TimePeriod) %>%
  summarise(Ccarb_Var=var(d13Ccarb))
b1<-Nit_DF %>%
  group_by(name, TimePeriod) %>%
  summarise(NColl_Var=var(d15Ncoll),CColl_Var=var(d13Ccoll))
b2<-CarC_DF %>%
  group_by(name, TimePeriod) %>%
  summarise(Ccarb_Var=var(d13Ccarb))

#Does variance (foraging breadth) change across the two time periods - non-parametric test
#M1<-lmer(CC~1+(1|SXTP))
#boxplot using ggplot2
library(reshape)
library(ggplot2)
library(dplyr)
a3.melt<-melt(a3[,c('NColl_Var','CColl_Var','TimePeriod','Sex')])
names(a1.melt)[5]<-"isotope_value"
names(a1.melt)[4]<-"isotope"
a1.melt<-a1.melt[complete.cases(a1.melt$isotope),]
newdf<-a1.melt[!a1.melt$isotope=="Ccarb_Var",]
F1<-newdf %>%
  ggplot(aes(x=TimePeriod, y=isotope_value, fill=TimePeriod)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) + theme_grey() + 
  facet_wrap(~isotope)
F1


###All correlations by site
head(Nit_DF)
library(dplyr)
cor_res <- Nit_DF %>%
  group_by(name) %>%
  do(cormat = cor(select(., -matches("name"))))

cormat_res <- iris %>%
  group_by(Species) %>%
  do(cormat = cor(select(., -matches("Species"))))
