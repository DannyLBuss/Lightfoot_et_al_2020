###MAM Models - Lightfoot et al. 2019 
##LM - Sibenik, Koprivno, Stenjavec, Zvon
##EM - Novak, Sisak, Zumb, Udbina DG, Koprivno 

#Load axis labels
N_label<-expression(paste(delta^15, "N collagen (\u2030)",sep=""))
C_label<-expression(paste(delta^13, "C collagen (\u2030)",sep=""))  
Cc_label<-expression(paste(delta^13, "C carbonate (\u2030)",sep=""))     
#xlabuc<-expression(paste("TimePeriod \n Uncorrected"))

#set MCMC parameters
#no. of interations
nitt <- c(500000)
#length of burnin
burnin <- c(10000)
#amount of thinning
thin <- c(245)

#Set priors
yfix <- diag(10) * 1e+03
sites_prior<- list(B = list(mu = rep(0, 10),V = yfix), 
                      R = list(V = 1, nu = 0.002))

yfix <- diag(7) * 1e+03
sites_prior_cc<- list(B = list(mu = rep(0, 7),V = yfix), 
                   R = list(V = 1, nu = 0.002))

#Run top models
SITE_DF$name<-as.character(SITE_DF$name)
SITE_DF$name<-as.factor(SITE_DF$name)
#Nitrogen collagen
S_nitt<-MCMCglmm(d15Ncoll~name + TimePeriod,
             rcov=~ units, 
             family="gaussian",
             data=SITE_DF,
             nitt=nitt,
             burnin=burnin,
             thin=thin,
             prior=sites_prior)
#DIC - site and timeperiod 547.9517
#DIC - site, timeperiod and sex 546.2291
#DIC - site 548.2196
#DIC - null 869.8734

DG<-mcmc(S_nitt$Sol[,"(Intercept)"])
Sib<-mcmc(S_nitt$Sol[,"(Intercept)"]) + S_nitt$Sol[,"TimePeriodLate Mediaeval"] + S_nitt$Sol[,"nameSibenik-sv-Lovre"]
Zum<-mcmc(S_nitt$Sol[,"(Intercept)"]) + S_nitt$Sol[,"nameZumberak"]
Kop<-mcmc(S_nitt$Sol[,"(Intercept)"]) + S_nitt$Sol[,"nameKoprivno"]
KopLM<-mcmc(S_nitt$Sol[,"(Intercept)"]) + S_nitt$Sol[,"nameKoprivno"] + S_nitt$Sol[,"TimePeriodLate Mediaeval"]
Sis<-mcmc(S_nitt$Sol[,"(Intercept)"]) + S_nitt$Sol[,"nameSisak"]
Sten<-mcmc(S_nitt$Sol[,"(Intercept)"]) + S_nitt$Sol[,"nameStenjevec"] + S_nitt$Sol[,"TimePeriodLate Mediaeval"]
NR<-mcmc(S_nitt$Sol[,"(Intercept)"]) + S_nitt$Sol[,"nameNova-Raca"]
Udbina<-mcmc(S_nitt$Sol[,"(Intercept)"]) + S_nitt$Sol[,"nameUdbina"]
Zvon<-mcmc(S_nitt$Sol[,"(Intercept)"]) + S_nitt$Sol[,"TimePeriodLate Mediaeval"] + S_nitt$Sol[,"nameZvonimirovo"]
Nitt_sites<-data.frame(
  c(rep("Drinovci-Greblje (EM)",2000),rep("Sibenik (LM)",2000),rep("Koprivno (EM)",2000),rep("Koprivno (LM)",2000),rep("Sisak (EM)",2000),rep("Stenjevec (LM)",2000),rep("Zumberak (EM)",2000),rep("Nova-Raca (EM)",2000),rep("Udbina (LM)",2000),rep("Zvonimirovo (EM)",2000)),
  c(DG,Sib,Kop,KopLM,Sis,Sten,Zum,NR,Udbina,Zvon),
  rep("Nitrogen collagen", 20000))
names(Nitt_sites)<-c("Site","mcmc","isotope")
Sites_Nitt<-ggplot(Nitt_sites[1:16000,], aes(mcmc, fill=Site, colour=Site)) + geom_density(alpha=0.5) + theme_classic() + labs(x=N_label)


#Carbon collagen
S_ccoll<-MCMCglmm(d13Ccoll~name + TimePeriod,
             rcov=~ units, 
             family="gaussian",
             data=SITE_DF,
             nitt=nitt,
             burnin=burnin,
             thin=thin,
             prior=sites_prior)
#DIC - site and timeperiod ???
#DIC - site, timeperiod and sex ???
#DIC - site ???
#DIC - null ???

DG<-mcmc(S_ccoll$Sol[,"(Intercept)"])
Sib<-mcmc(S_ccoll$Sol[,"(Intercept)"]) + S_ccoll$Sol[,"TimePeriodLate Mediaeval"] + S_ccoll$Sol[,"nameSibenik-sv-Lovre"]
Zum<-mcmc(S_ccoll$Sol[,"(Intercept)"]) + S_ccoll$Sol[,"nameZumberak"]
Kop<-mcmc(S_ccoll$Sol[,"(Intercept)"]) + S_ccoll$Sol[,"nameKoprivno"]
KopLM<-mcmc(S_ccoll$Sol[,"(Intercept)"]) + S_ccoll$Sol[,"nameKoprivno"] + S_ccoll$Sol[,"TimePeriodLate Mediaeval"]
Sis<-mcmc(S_ccoll$Sol[,"(Intercept)"]) + S_ccoll$Sol[,"nameSisak"]
Sten<-mcmc(S_ccoll$Sol[,"(Intercept)"]) + S_ccoll$Sol[,"nameStenjevec"] + S_ccoll$Sol[,"TimePeriodLate Mediaeval"]
NR<-mcmc(S_ccoll$Sol[,"(Intercept)"]) + S_ccoll$Sol[,"nameNova-Raca"]
Udbina<-mcmc(S_ccoll$Sol[,"(Intercept)"]) + S_ccoll$Sol[,"nameUdbina"]
Zvon<-mcmc(S_ccoll$Sol[,"(Intercept)"]) + S_ccoll$Sol[,"TimePeriodLate Mediaeval"] + S_ccoll$Sol[,"nameZvonimirovo"]
Ccoll_sites<-data.frame(
  c(rep("Drinovci-Greblje (EM)",2000),rep("Sibenik (LM)",2000),rep("Koprivno (EM)",2000),rep("Koprivno (LM)",2000),rep("Sisak (EM)",2000),rep("Stenjevec (LM)",2000),rep("Zumberak (EM)",2000),rep("Nova-Raca (EM)",2000),rep("Udbina (LM)",2000),rep("Zvonimirovo (EM)",2000)),
  c(DG,Sib,Kop,KopLM,Sis,Sten,Zum,NR,Udbina,Zvon),
  rep("Carbon collagen", 20000))
names(Ccoll_sites)<-c("Site","mcmc","isotope")
Sites_Ccoll<-ggplot(Ccoll_sites[1:16000,], aes(mcmc, fill=Site, colour=Site)) + geom_density(alpha=0.5) + theme_classic() + labs(x=C_label)


#Carbon carbonate
SITE_DFCC<-SITE_DF[complete.cases(SITE_DF$d13Ccarb),]
SITE_DFCC$name<-as.character(SITE_DFCC$name)
SITE_DFCC$name<-as.factor(SITE_DFCC$name)
S_ccarb<-MCMCglmm(d13Ccarb~name + TimePeriod,
             rcov=~ units, 
             family="gaussian",
             data=SITE_DFCC,
             nitt=nitt,
             burnin=burnin,
             thin=thin,
             prior=sites_prior_cc)

DG<-mcmc(S_ccarb$Sol[,"(Intercept)"])
Sib<-mcmc(S_ccarb$Sol[,"(Intercept)"]) + S_ccarb$Sol[,"TimePeriodLate Mediaeval"] + S_ccarb$Sol[,"nameSibenik-sv-Lovre"]
Zum<-mcmc(S_ccarb$Sol[,"(Intercept)"]) + S_ccarb$Sol[,"nameZumberak"]
Kop<-mcmc(S_ccarb$Sol[,"(Intercept)"]) + S_ccarb$Sol[,"nameKoprivno"]
KopLM<-mcmc(S_ccarb$Sol[,"(Intercept)"]) + S_ccarb$Sol[,"nameKoprivno"] + S_ccarb$Sol[,"TimePeriodLate Mediaeval"]
Sis<-mcmc(S_ccarb$Sol[,"(Intercept)"]) + S_ccarb$Sol[,"nameSisak"]
Sten<-mcmc(S_ccarb$Sol[,"(Intercept)"]) + S_ccarb$Sol[,"nameStenjevec"] + S_ccarb$Sol[,"TimePeriodLate Mediaeval"]
ccarb_sites<-data.frame(
  c(rep("Drinovci-Greblje (EM)",2000),rep("Sibenik (LM)",2000),rep("Koprivno (EM)",2000),rep("Koprivno (LM)",2000),rep("Sisak (EM)",2000),rep("Stenjevec (LM)",2000),rep("Zumberak (EM)",2000)),
  c(DG,Sib,Kop,KopLM,Sis,Sten,Zum),
  rep("Carbon carbonate", 14000))
names(ccarb_sites)<-c("Site","mcmc","isotope")
Sites_CCarb<-ggplot(ccarb_sites[1:14000,], aes(mcmc, fill=Site, colour=Site)) + geom_density(alpha=0.5) + theme_classic() + labs(x=Cc_label)

Output_sites<-rbind(Nitt_sites,Ccoll_sites,ccarb_sites)
library(ggthemr)
ggthemr(palette="dust")
cls<-c("#555555","#db735c","#EFA86E","#9A8A76","#F3C57B","#7A6752","#2A91A2","#87F28A","#6EDCEF","#555567","#db737c","#EFA85E","#9A8A46","#F3C56B","#7A6732","#2A91A1","#67F28A","#4EDCEF")
set_swatch(cls)
tiff("all_sites_models.tiff", unit="in", width=11, height=3.5,res=300)
ggplot(Output_sites, aes(mcmc, fill=Site, colour=Site)) + geom_density(alpha=0.5) + theme_classic() + labs(x=Cc_label) + facet_wrap(.~isotope, scales="free")
dev.off()

#no. of interations
nitt <- c(2000000)
#length of burnin
burnin <- c(40000)
#amount of thinning
thin <- c(980)