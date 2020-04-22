library(MCMCglmm)
library(ggplot2)
Site_prior<-list(B = list(mu= diag(11)*0, V=diag(11)*1e+10),
                R = list(nu=0.002, V=1),
                G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

SITE_DF<-Nit_DF[!Nit_DF$name=="Josipovo",]

yfix <- diag(8) * 1e+03
sites_prior<- list(B = list(mu = rep(0, 8),V = yfix), 
                   R = list(V = 1, nu = 0.002))

PS<-MCMCglmm(d13Ccarb~name + Sex,
         random=~TimePeriod,
         rcov=~ units, 
         family="gaussian",
         data=SITE_DF,
         nitt=nitt,
         burnin=burnin,
         thin=thin,
         prior=sites_prior)
plot(PS)

Site_prior<-list(B = list(mu= diag(7)*0, V=diag(7)*1e+10),
                 R = list(nu=0.002, V=1),
                 G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

unique(SITE_DFCC$name)
SITE_DFCC<-SITE_DF[!SITE_DF$name=="Zvonimirovo" & !SITE_DF$name=="Udbina" & !SITE_DF$name=="Nova-Raca",]
SITE_DFCC$name<-as.character(SITE_DFCC$name)
SITE_DFCC$name<-as.factor(SITE_DFCC$name)

yfix <- diag(7) * 1e+03
sites_prior<- list(B = list(mu = rep(0, 7),V = yfix), 
                   R = list(V = 1, nu = 0.002))
PS<-MCMCglmm(d13Ccarb~name + TimePeriod,
             rcov=~ units, 
             family="gaussian",
             data=SITE_DFCC,
             nitt=nitt,
             burnin=burnin,
             thin=thin,
             prior=sites_prior)
PS$DIC

Site_prior<-list(B = list(mu= diag(6)*0, V=diag(6)*1e+10),
                 R = list(nu=0.002, V=1),
                 G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))
#Best carbonate sites model
PS_name<-MCMCglmm(d13Ccarb~name,
             ~ idh(TimePeriod),
             rcov=~ units, 
             family="gaussian",
             data=SITE_DFCC,
             nitt=nitt,
             burnin=burnin,
             thin=thin,
             prior=Site_prior)

Site_prior<-list(B = list(mu= diag(1)*0, V=diag(1)*1e+10),
                 R = list(nu=0.002, V=1),
                 G = list(G1 = list(V = diag(2)*1, nu = 0.02, alpha.mu = rep(0,2), alpha.V = diag(2)*1000)))

PS_null<-MCMCglmm(d13Ccarb~1,
                  ~ idh(TimePeriod),
                  rcov=~ units, 
                  family="gaussian",
                  data=SITE_DFCC,
                  nitt=nitt,
                  burnin=burnin,
                  thin=thin,
                  prior=Site_prior)

#PS$DIC = 596.3109
#PS_name$DIC = 594.5513
#PS_name_TP$DIC = 594.5575
#PS_null = 745.835

##LM - Sibenik, some Koprivno, Stenjavec, Zvon
##EM - Novak, Sisak, Zumb, Udbina DG, Koprivno 
#no. of interations
nitt <- c(2000000)
#length of burnin
burnin <- c(40000)
#amount of thinning
thin <- c(980)
PS<-MCMCglmm(d13Ccarb~name + TimePeriod,
             rcov=~ units, 
             family="gaussian",
             data=SITE_DFCC,
             nitt=nitt,
             burnin=burnin,
             thin=thin,
             prior=sites_prior)

DG<-mcmc(PS$Sol[,"(Intercept)"])
#DG_M<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"SexMale"]
#Zvon<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"TimePeriodLate Mediaeval"] + PS$Sol[,"nameZvonimirovo"]
#Zvon_M<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"SexMale"] + PS$Sol[,"TimePeriodLate Mediaeval"] + PS$Sol[,"nameZvonimirovo"]
Sib<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"TimePeriodLate Mediaeval"] + PS$Sol[,"nameSibenik-sv-Lovre"]
#Sib_M<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"SexMale"] + PS$Sol[,"TimePeriodLate Mediaeval"] + PS$Sol[,"nameSibenik-sv-Lovre"]
Zum<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"nameZumberak"]
#Zum_M<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"SexMale"] + PS$Sol[,"nameZumberak"]
#Udbina<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"nameUdbina"]
#Udbina_M<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"SexMale"] + PS$Sol[,"nameUdbina"]
Kop<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"nameKoprivno"]
#Kop_M<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"SexMale"] + PS$Sol[,"nameKoprivno"]
KopLM<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"nameKoprivno"] + PS$Sol[,"TimePeriodLate Mediaeval"]
#KopLM_M<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"SexMale"] + PS$Sol[,"nameKoprivno"] + PS$Sol[,"TimePeriodLate Mediaeval"]
#NR<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"nameNova-Raca"]
#NR_M<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"SexMale"] + PS$Sol[,"nameNova-Raca"]
Sis<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"nameSisak"]
#Sis_M<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"SexMale"] + PS$Sol[,"nameSisak"]
Sten<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"nameStenjevec"] + PS$Sol[,"TimePeriodLate Mediaeval"]
#Sten_M<-mcmc(PS$Sol[,"(Intercept)"]) + PS$Sol[,"SexMale"] + PS$Sol[,"nameStenjevec"] + PS$Sol[,"TimePeriodLate Mediaeval"]
ccarb_sites<-data.frame(
  c(rep("DG",2000),rep("Sib",2000),rep("Zum",2000),rep("Koprivno_EM",2000),rep("Koprivno_LM",2000),rep("Sisak",2000),rep("Sten",2000)),
  c(DG,Sib,Zum,Kop,KopLM,Sis,Sten))
names(ccarb_sites)<-c("TP","mcmc")
ccarb<-ggplot(ccarb_sites[1:14000,], aes(mcmc, fill=TP, colour=TP)) + geom_density(alpha=0.5) + theme_classic() + labs(x=xlabel)
ylabel<-expression(paste(delta^15, "N (\u2030)",sep=""))
xlabel<-expression(paste(delta^13, "C (\u2030)",sep=""))     
xlabuc<-expression(paste("TimePeriod \n Uncorrected"))

#nit_sites<-data.frame(
  c(rep("DG_F",2000),rep("DG_M",2000),rep("Zvon_F",2000),rep("Zvon_M",2000),rep("Sib_F",2000),rep("Sib_M",2000),
    rep("Zum_F",2000),rep("Zum_M",2000),rep("Udbina_F",2000),rep("Udbina_M",2000),rep("Kop_F",2000),rep("Kop_M",2000),
    rep("KopLM_F",2000),rep("KopLM_M",2000),rep("NR_F",2000),rep("NR_M",2000),rep("Sis_F",2000),rep("Sis_M",2000),
    rep("Sten_F",2000),rep("Sten_M",2000)),
  c(DG_F,DG_M,Zvon_F,Zvon_M,Sib_F,Sib_M,Zum_F,Zum_M,Udbina_F,Udbina_M,Kop_F,Kop_M,KopLM_F,KopLM_M,NR_F,NR_M,Sis_F,Sis_M,Sten_F,Sten_M))
names(nit_sites)<-c("TP","mcmc")

library(ggthemr)
set_swatch()
a<-ggplot(nit_sites[1:14000,], aes(mcmc, fill=TP, colour=TP)) + geom_density(alpha=0.5) + theme_classic() + scale_colour_manual(values = c("grey10","grey15","grey20","grey25","grey30","grey40","grey45","grey50","grey55","grey60","grey65","grey70","grey75","grey80","grey85","grey90","grey95","grey99","red","blue"))
