###Creating plots for Emma's paper
###Priors and MCMC settings:
##Set MCMC parameters
#no. of interations
nitt <- c(500000)
#length of burnin
burnin <- c(10000)
#amount of thinning
thin <- c(245)

datanew<-Nit_DF[!Nit_DF$name=="Sisak",]


#Max model:
Max_prior<-list(B = list(mu= diag(2)*0, V=diag(2)*1e+10),
                R = list(nu=0.002, V=1),
                G = list(G1 = list(V = diag(9)*1, nu = 0.02, alpha.mu = rep(0,9), alpha.V = diag(9)*1000)))
###Paired-site model:
#subset data
paired<-Nit_DF[Nit_DF$name=="Koprivno",]
paired<-paired[,c(2,5,7,9:11)]
yfixed <- diag(4) * 1e+03
ps_prior<- list(B = list(mu = rep(0, 4),V = yfixed), 
                R = list(V = 1, nu = 0.002))

#Sites model:
#collagen
yfix <- diag(19) * 1e+03
sites_prior<- list(B = list(mu = rep(0, 19),V = yfix), 
                   R = list(V = 1, nu = 0.002))
#carbonate
yfix_b <- diag(13) * 1e+03
sites_prior2<- list(B = list(mu = rep(0, 13),V = yfix_b), 
                    R = list(V = 1, nu = 0.002))

yfixed <- diag(4) * 1e+03
ps_prior<- list(B = list(mu = rep(0, 4),V = yfixed), 
                R = list(V = 1, nu = 0.002))
Nit_Max<-MCMCglmm(d15Ncoll~TimePeriod,
         ~ idh(name),
         rcov=~ units, 
         family="gaussian",
         data=datanew,
         nitt=nitt,
         burnin=burnin,
         thin=thin,
         prior=Max_prior)

yfixed <- diag(3) * 1e+03
ps_prior<- list(B = list(mu = rep(0, 3),V = yfixed), 
                R = list(V = 1, nu = 0.002))
Nit_Max<-MCMCglmm(d15Ncoll~TimePeriod,
                  ~ idh(name),
                  rcov=~ units, 
                  family="gaussian",
                  data=datanew,
                  nitt=nitt,
                  burnin=burnin,
                  thin=thin,
                  prior=Max_prior)

mean_EM<-mean(mcmc(Nit_Max$Sol[,"(Intercept)"]))
mean_LM<-mean(mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"TimePeriodLate Mediaeval"])
#mean_male_LM<-mean(mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"TimePeriodLate Mediaeval"]

PS_ccarb<-MCMCglmm(d15Ncoll~TimePeriod + Sex,
         rcov=~ units, 
         family="gaussian",
         data=paired,
         nitt=nitt,
         burnin=burnin,
         thin=thin,
         prior=ps_prior)

#collagen
yfix <- diag(12) * 1e+03
sites_prior<- list(B = list(mu = rep(0, 12),V = yfix), 
                   R = list(V = 1, nu = 0.002))
#carbonate
yfix_b <- diag(13) * 1e+03
sites_prior2<- list(B = list(mu = rep(0, 13),V = yfix_b), 
                    R = list(V = 1, nu = 0.002))

PS_ccarb<-MCMCglmm(d15Ncoll~name + Sex + TimePeriod-1,
         family="gaussian",
         data=Nit_DF,
         nitt=nitt,
         burnin=burnin,
         thin=thin,
         prior=sites_prior)
unique(Nit_DF$name)

S_Sisak<-mcmc(PS_ccarb$Sol[,"nameSisak"])
S_Nova-Raca<-mcmc(PS_ccarb$Sol[,"nameNova-Raca"]) + PS_ccarb$Sol[,"SexMale"]
mean_male_LM<-mcmc(PS_ccarb$Sol[,"(Intercept)"]) + PS_ccarb$Sol[,"TimePeriodLate Mediaeval"] + PS_ccarb$Sol[,"SexMale"]
mean_female_LM<-mcmc(PS_ccarb$Sol[,"(Intercept)"]) + PS_ccarb$Sol[,"TimePeriodLate Mediaeval"]
ccarb_koprivno<-data.frame(
  c(rep("mean_female_EM",2000),rep("mean_male_EM",2000),rep("mean_male_LM",2000),rep("mean_female_LM",2000)),
    c(mean_female_EM,mean_male_EM,mean_male_LM,mean_female_LM))
names(ccarb_koprivno)<-c("TP","mcmc")
dev.off()
ggplot(ccarb_koprivno, aes(mcmc, fill=TP, colour=TP)) + geom_density(alpha=0.5) + theme_classic()

mean_female_EM<-mean(mcmc(Nit_Max$Sol[,"(Intercept)"]))
mean_male_EM<-mean(mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"SexMale"])
mean_male_LM<-mean(mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"TimePeriodLate Mediaeval"] + Nit_Max$Sol[,"SexMale"])
mean_female_LM<-mean(mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"TimePeriodLate Mediaeval"])

lower_HPD_female_EM<-HPDinterval(mcmc(Nit_Max$Sol[,"(Intercept)"]))[,1][1]
lower_HPD_male_EM<-HPDinterval(mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"SexMale"])[,1][1]
lower_HPD_male_LM<-HPDinterval(mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"TimePeriodLate Mediaeval"] + Nit_Max$Sol[,"SexMale"])[,1][1]
lower_HPD_female_LM<-HPDinterval(mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"TimePeriodLate Mediaeval"])[,1][1]

upper_HPD_female_EM<-HPDinterval(mcmc(Nit_Max$Sol[,"(Intercept)"]))[,2][1]
upper_HPD_male_EM<-HPDinterval(mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"SexMale"])[,2][1]
upper_HPD_male_LM<-HPDinterval(mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"TimePeriodLate Mediaeval"] + Nit_Max$Sol[,"SexMale"])[,2][1]
upper_HPD_female_LM<-HPDinterval(mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"TimePeriodLate Mediaeval"])[,2][1]

#CColl
Cc_mean_female_EM<-mean(mcmc(Ccoll_Max$Sol[,"(Intercept)"]))
Cc_mean_male_EM<-mean(mcmc(Ccoll_Max$Sol[,"(Intercept)"]) + Ccoll_Max$Sol[,"SexMale"])
Cc_mean_male_LM<-mean(mcmc(Ccoll_Max$Sol[,"(Intercept)"]) + Ccoll_Max$Sol[,"TimePeriodLate Mediaeval"] + Ccoll_Max$Sol[,"SexMale"])
Cc_mean_female_LM<-mean(mcmc(Ccoll_Max$Sol[,"(Intercept)"]) + Ccoll_Max$Sol[,"TimePeriodLate Mediaeval"])

Cc_lower_HPD_female_EM<-HPDinterval(mcmc(Ccoll_Max$Sol[,"(Intercept)"]))[,1][1]
Cc_lower_HPD_male_EM<-HPDinterval(mcmc(Ccoll_Max$Sol[,"(Intercept)"]) + Ccoll_Max$Sol[,"SexMale"])[,1][1]
Cc_lower_HPD_male_LM<-HPDinterval(mcmc(Ccoll_Max$Sol[,"(Intercept)"]) + Ccoll_Max$Sol[,"TimePeriodLate Mediaeval"] + Ccoll_Max$Sol[,"SexMale"])[,1][1]
Cc_lower_HPD_female_LM<-HPDinterval(mcmc(Ccoll_Max$Sol[,"(Intercept)"]) + Ccoll_Max$Sol[,"TimePeriodLate Mediaeval"])[,1][1]

Cc_upper_HPD_female_EM<-HPDinterval(mcmc(Ccoll_Max$Sol[,"(Intercept)"]))[,2][1]
Cc_upper_HPD_male_EM<-HPDinterval(mcmc(Ccoll_Max$Sol[,"(Intercept)"]) + Ccoll_Max$Sol[,"SexMale"])[,2][1]
Cc_upper_HPD_male_LM<-HPDinterval(mcmc(Ccoll_Max$Sol[,"(Intercept)"]) + Ccoll_Max$Sol[,"TimePeriodLate Mediaeval"] + Nit_Max$Sol[,"SexMale"])[,2][1]
Cc_upper_HPD_female_LM<-HPDinterval(mcmc(Ccoll_Max$Sol[,"(Intercept)"]) + Ccoll_Max$Sol[,"TimePeriodLate Mediaeval"])[,2][1]

#CCarb
Ccb_mean_female_EM<-mean(mcmc(Ccarb_Max$Sol[,"(Intercept)"]))
Ccb_mean_male_EM<-mean(mcmc(Ccarb_Max$Sol[,"(Intercept)"]) + Ccarb_Max$Sol[,"SexMale"])
Ccb_mean_male_LM<-mean(mcmc(Ccarb_Max$Sol[,"(Intercept)"]) + Ccarb_Max$Sol[,"TimePeriodLate Mediaeval"] + Ccarb_Max$Sol[,"SexMale"])
Ccb_mean_female_LM<-mean(mcmc(Ccarb_Max$Sol[,"(Intercept)"]) + Ccarb_Max$Sol[,"TimePeriodLate Mediaeval"])

Ccb_lower_HPD_female_EM<-HPDinterval(mcmc(Ccarb_Max$Sol[,"(Intercept)"]))[,1][1]
Ccb_lower_HPD_male_EM<-HPDinterval(mcmc(Ccarb_Max$Sol[,"(Intercept)"]) + Ccarb_Max$Sol[,"SexMale"])[,1][1]
Ccb_lower_HPD_male_LM<-HPDinterval(mcmc(Ccarb_Max$Sol[,"(Intercept)"]) + Ccarb_Max$Sol[,"TimePeriodLate Mediaeval"] + Ccarb_Max$Sol[,"SexMale"])[,1][1]
Ccb_lower_HPD_female_LM<-HPDinterval(mcmc(Ccarb_Max$Sol[,"(Intercept)"]) + Ccarb_Max$Sol[,"TimePeriodLate Mediaeval"])[,1][1]

Ccb_upper_HPD_female_EM<-HPDinterval(mcmc(Ccarb_Max$Sol[,"(Intercept)"]))[,2][1]
Ccb_upper_HPD_male_EM<-HPDinterval(mcmc(Ccarb_Max$Sol[,"(Intercept)"]) + Ccarb_Max$Sol[,"SexMale"])[,2][1]
Ccb_upper_HPD_male_LM<-HPDinterval(mcmc(Ccarb_Max$Sol[,"(Intercept)"]) + Ccarb_Max$Sol[,"TimePeriodLate Mediaeval"] + Nit_Max$Sol[,"SexMale"])[,2][1]
Ccb_upper_HPD_female_LM<-HPDinterval(mcmc(Ccarb_Max$Sol[,"(Intercept)"]) + Ccarb_Max$Sol[,"TimePeriodLate Mediaeval"])[,2][1]



names_NIT<-c("female_EM","male_EM","male_LM","female_LM")

plotNIT_df<-data.frame(names_NIT,
                       c(mean_female_EM,mean_male_EM,mean_male_LM,mean_female_LM,Cc_mean_female_EM,Cc_mean_male_EM,Cc_mean_male_LM,Cc_mean_female_LM,Ccb_mean_female_EM,Ccb_mean_male_EM,Ccb_mean_male_LM,Ccb_mean_female_LM),
                       c(lower_HPD_female_EM,lower_HPD_male_EM,lower_HPD_male_LM,lower_HPD_female_LM,Cc_lower_HPD_female_EM,Cc_lower_HPD_male_EM,Cc_lower_HPD_male_LM,Cc_lower_HPD_female_LM,Ccb_lower_HPD_female_EM,Ccb_lower_HPD_male_EM,Ccb_lower_HPD_male_LM,Ccb_lower_HPD_female_LM),
                       c(upper_HPD_female_EM,upper_HPD_male_EM,upper_HPD_male_LM,upper_HPD_female_LM,Cc_upper_HPD_female_EM,Cc_upper_HPD_male_EM,Cc_upper_HPD_male_LM,Cc_upper_HPD_female_LM,Ccb_upper_HPD_female_EM,Ccb_upper_HPD_male_EM,Ccb_upper_HPD_male_LM,Ccb_upper_HPD_female_LM),
                       rep(c("Early Modern","Early Modern","Late Mediaeval","Late Mediaeval"),3),
                       rep(c("Female","Male","Male","Female"),3,),
                       c(rep("NitrogenCollagen",4),rep("CarbonCollagen",4),rep("CarbonCarbonate",4))
)

names(plotNIT_df)<-c("names_NIT","means","lowers","uppers","TimePeriod","Sex","Isotope")
str(plotNIT_df)
library(ggplot2)
man_col<-c("grey50","grey2")
ggplot(plotNIT_df, aes(x=names_NIT,y=means, fill=TimePeriod)) + geom_bar(stat="identity") + facet_wrap(.~Isotope) + geom_errorbar(aes(ymin= lowers, ymax = uppers)) + scale_colour_grey()
ggplot(plotNIT_df[plotNIT_df$Isotope=="NitrogenCollagen",], aes(x=names_NIT,y=means, fill=TimePeriod)) + geom_bar(stat="identity") + facet_wrap(.~Isotope) + geom_errorbar(aes(ymin= lowers, ymax = uppers)) + scale_colour_grey()

Ccoll_Max<-MCMCglmm(d13Ccoll~TimePeriod,
                  ~ idh(name),
                  rcov=~ units, 
                  family="gaussian",
                  data=datanew,
                  nitt=nitt,
                  burnin=burnin,
                  thin=thin,
                  prior=Max_prior)

Ccarb_Max<-MCMCglmm(d13Ccarb~TimePeriod,
                  ~ idh(name),
                  rcov=~ units, 
                  family="gaussian",
                  data=datanew,
                  nitt=nitt,
                  burnin=burnin,
                  thin=thin,
                  prior=Max_prior)

convergence_checks_max<-list(gelman.diag(M_Nitrogen),
                             gelman.diag(M_CarbCol),
                             gelman.diag(M_CCarb))

summary(M_CCarb)
#Check tracplots
plot(M_Nitrogen[1][1])
plot(M_Nitrogen[2][1])
plot(M_Nitrogen[3][1])
plot(M_CarbCol[1][1])
plot(M_CarbCol[2][1])
plot(M_CarbCol[3][1])
plot(M_CCarb[1][1])
plot(M_CCarb[2][1])
plot(M_CCarb[3][1])

plot.estimates(M_Nitrogen)
plot.estimates(M_CarbCol)
plot.estimates(M_CCarb)
