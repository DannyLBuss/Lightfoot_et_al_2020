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
TP_prior<-list(B = list(mu= diag(3)*0, V=diag(3)*1e+10),
                 R = list(nu=0.002, V=1),
                 G = list(G1 = list(V = diag(9)*1, nu = 0.02, alpha.mu = rep(0,9), alpha.V = diag(9)*1000)))

datanew<-datanew[!datanew$name=="Sisak",]
Nit_Max<-MCMCglmm(d15Ncoll~TimePeriod + Sex,
                  ~ idh(name),
                  rcov=~ units, 
                  family="gaussian",
                  data=datanew,
                  nitt=nitt,
                  burnin=burnin,
                  thin=thin,
                  prior=TP_prior)

EM_F<-mcmc(Nit_Max$Sol[,"(Intercept)"])
EM_M<-mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"SexMale"]
LM_F<-mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"TimePeriodLate Mediaeval"]
LM_M<-mcmc(Nit_Max$Sol[,"(Intercept)"]) + Nit_Max$Sol[,"SexMale"] + Nit_Max$Sol[,"TimePeriodLate Mediaeval"]
Nitt_TP<-data.frame(
  c(rep("Female in Early Modern",2000),rep("Male in Early Modern",2000),rep("Female in Late Mediaeval",2000),rep("Male in Late Mediaeval",2000)),
  c(EM_F,EM_M,LM_F,LM_M),
  rep("Nitrogen collagen", 8000))
names(Nitt_TP)<-c("Sex_TimePeriod","mcmc","isotope")
ggplot(Nitt_TP[1:8000,], aes(mcmc, fill=Sex_TimePeriod, colour=Sex_TimePeriod)) + geom_density(alpha=0.5) + theme_classic() + labs(x=N_label)

##Carbon collagen
TP_prior<-list(B = list(mu= diag(3)*0, V=diag(3)*1e+10),
               R = list(nu=0.002, V=1),
               G = list(G1 = list(V = diag(9)*1, nu = 0.02, alpha.mu = rep(0,9), alpha.V = diag(9)*1000)))

datanew$name<-as.character(datanew$name)
datanew$name<-as.factor(datanew$name)
Carbon_Max<-MCMCglmm(d13Ccoll~TimePeriod + Sex,
                  ~ idh(name),
                  rcov=~ units, 
                  family="gaussian",
                  data=datanew,
                  nitt=nitt,
                  burnin=burnin,
                  thin=thin,
                  prior=TP_prior)
#TP + Sex - 637.5809
#TP - 635.6237
#plot(Carbon_Max)

EM_F<-mcmc(Carbon_Max$Sol[,"(Intercept)"])
EM_M<-mcmc(Carbon_Max$Sol[,"(Intercept)"]) + Carbon_Max$Sol[,"SexMale"]
LM_F<-mcmc(Carbon_Max$Sol[,"(Intercept)"]) + Carbon_Max$Sol[,"TimePeriodLate Mediaeval"]
LM_M<-mcmc(Carbon_Max$Sol[,"(Intercept)"]) + Carbon_Max$Sol[,"SexMale"] + Carbon_Max$Sol[,"TimePeriodLate Mediaeval"]
Carbon_TP<-data.frame(
  c(rep("Female in Early Modern",2000),rep("Male in Early Modern",2000),rep("Female in Late Mediaeval",2000),rep("Male in Late Mediaeval",2000)),
  c(EM_F,EM_M,LM_F,LM_M),
  rep("Carbon collagen", 8000))
names(Carbon_TP)<-c("TimePeriod","mcmc","isotope")
ggplot(Carbon_TP[1:8000,], aes(mcmc, fill=TimePeriod, colour=TimePeriod)) + geom_density(alpha=0.5) + theme_classic() + labs(x=C_label)

##Carbon carbonate
datanew<-datanew[complete.cases(datanew$d13Ccarb),]
datanew$name<-as.character(datanew$name)
datanew$name<-as.factor(datanew$name)
TP_prior_cc<-list(B = list(mu= diag(3)*0, V=diag(3)*1e+10),
               R = list(nu=0.002, V=1),
               G = list(G1 = list(V = diag(5)*1, nu = 0.02, alpha.mu = rep(0,5), alpha.V = diag(5)*1000)))
CarbonC_Max<-MCMCglmm(d13Ccarb~TimePeriod + Sex,
                     ~ idh(name),
                     rcov=~ units, 
                     family="gaussian",
                     data=datanew,
                     nitt=nitt,
                     burnin=burnin,
                     thin=thin,
                     prior=TP_prior_cc)

EM_F<-mcmc(CarbonC_Max$Sol[,"(Intercept)"])
EM_M<-mcmc(CarbonC_Max$Sol[,"(Intercept)"]) + CarbonC_Max$Sol[,"SexMale"]
LM_F<-mcmc(CarbonC_Max$Sol[,"(Intercept)"]) + CarbonC_Max$Sol[,"TimePeriodLate Mediaeval"]
LM_M<-mcmc(CarbonC_Max$Sol[,"(Intercept)"]) + CarbonC_Max$Sol[,"SexMale"] + CarbonC_Max$Sol[,"TimePeriodLate Mediaeval"]
CarbonC_TP<-data.frame(
  c(rep("Female in Early Modern",2000),rep("Male in Early Modern",2000),rep("Female in Late Mediaeval",2000),rep("Male in Late Mediaeval",2000)),
  c(EM_F,EM_M,LM_F,LM_M),
  rep("Carbon carbonate", 8000))
names(CarbonC_TP)<-c("TimePeriod","mcmc","isotope")
ggplot(CarbonC_TP[1:8000,], aes(mcmc, fill=TimePeriod, colour=TimePeriod)) + geom_density(alpha=0.5) + theme_classic() + labs(x=C_label)

Output_TP<-rbind(Nitt_TP,Carbon_TP,CarbonC_TP)
tiff("tp_sex_models.tiff", unit="in", width=11, height=3.5,res=300)
ggplot(Output_TP, aes(mcmc, fill=TimePeriod, colour=TimePeriod)) + geom_density(alpha=0.5) + theme_classic() + labs(x=Cc_label) + facet_wrap(.~isotope, scales="free")
dev.off()
