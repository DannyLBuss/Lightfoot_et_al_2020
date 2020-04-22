###Paired Site model outputs for Appendix 5
N_label<-expression(paste(delta^15, "N collagen (\u2030)",sep=""))
C_label<-expression(paste(delta^13, "C collagen (\u2030)",sep=""))  
Cc_label<-expression(paste(delta^13, "C carbonate (\u2030)",sep="")) 

#set MCMC parameters
#no. of interations
nitt <- c(1000000)
#length of burnin
burnin <- c(20000)
#amount of thinning
thin <- c(490)

#Set priors
P_prior<-list(B = list(mu= diag(3)*0, V=diag(3)*1e+10),
               R = list(nu=0.002, V=1))
#Subset data
datapaired<-datanew[datanew$name=="Koprivno",]

#Nitrogen collagen
Nit_paired<-MCMCglmm(d15Ncoll~TimePeriod + Sex,
                  rcov=~ units, 
                  family="gaussian",
                  data=datanew,
                  nitt=nitt,
                  burnin=burnin,
                  thin=thin,
                  prior=P_prior)
EM_F<-mcmc(Nit_paired$Sol[,"(Intercept)"])
EM_M<-mcmc(Nit_paired$Sol[,"(Intercept)"]) + Nit_paired$Sol[,"SexMale"]
LM_F<-mcmc(Nit_paired$Sol[,"(Intercept)"]) + Nit_paired$Sol[,"TimePeriodLate Mediaeval"]
LM_M<-mcmc(Nit_paired$Sol[,"(Intercept)"]) + Nit_paired$Sol[,"SexMale"] + Nit_paired$Sol[,"TimePeriodLate Mediaeval"]
Nitt_P<-data.frame(
  c(rep("Female in Early Modern",2000),rep("Male in Early Modern",2000),rep("Female in Late Mediaeval",2000),rep("Male in Late Mediaeval",2000)),
  c(EM_F,EM_M,LM_F,LM_M),
  rep("Nitrogen collagen", 8000))
names(Nitt_P)<-c("Koprivno","mcmc","isotope")
ggplot(Nitt_P[1:8000,], aes(mcmc, fill=Sex_TimePeriod, colour=Sex_TimePeriod)) + geom_density(alpha=0.5) + theme_classic() + labs(x=N_label)

#Carbon collagen
C_paired<-MCMCglmm(d13Ccoll~TimePeriod + Sex,
                     rcov=~ units, 
                     family="gaussian",
                     data=datanew,
                     nitt=nitt,
                     burnin=burnin,
                     thin=thin,
                     prior=P_prior)

EM_F<-mcmc(C_paired$Sol[,"(Intercept)"])
EM_M<-mcmc(C_paired$Sol[,"(Intercept)"]) + C_paired$Sol[,"SexMale"]
LM_F<-mcmc(C_paired$Sol[,"(Intercept)"]) + C_paired$Sol[,"TimePeriodLate Mediaeval"]
LM_M<-mcmc(C_paired$Sol[,"(Intercept)"]) + C_paired$Sol[,"SexMale"] + C_paired$Sol[,"TimePeriodLate Mediaeval"]
Carb_P<-data.frame(
  c(rep("Female in Early Modern",2000),rep("Male in Early Modern",2000),rep("Female in Late Mediaeval",2000),rep("Male in Late Mediaeval",2000)),
  c(EM_F,EM_M,LM_F,LM_M),
  rep("Carbon collagen", 8000))
names(Carb_P)<-c("Koprivno","mcmc","isotope")
ggplot(Carb_P[1:8000,], aes(mcmc, fill=Sex_TimePeriod, colour=Sex_TimePeriod)) + geom_density(alpha=0.5) + theme_classic() + labs(x=C_label)

#Carbon carbonate
Cc_paired<-MCMCglmm(d13Ccarb~TimePeriod + Sex,
                   rcov=~ units, 
                   family="gaussian",
                   data=datanew,
                   nitt=nitt,
                   burnin=burnin,
                   thin=thin,
                   prior=P_prior)
traceplot(Cc_paired$Sol)
EM_F<-mcmc(Cc_paired$Sol[,"(Intercept)"])
EM_M<-mcmc(Cc_paired$Sol[,"(Intercept)"]) + Cc_paired$Sol[,"SexMale"]
LM_F<-mcmc(Cc_paired$Sol[,"(Intercept)"]) + Cc_paired$Sol[,"TimePeriodLate Mediaeval"]
LM_M<-mcmc(Cc_paired$Sol[,"(Intercept)"]) + Cc_paired$Sol[,"SexMale"] + Cc_paired$Sol[,"TimePeriodLate Mediaeval"]
CarbC_P<-data.frame(
  c(rep("Female in Early Modern",2000),rep("Male in Early Modern",2000),rep("Female in Late Mediaeval",2000),rep("Male in Late Mediaeval",2000)),
  c(EM_F,EM_M,LM_F,LM_M),
  rep("Carbon carbonate", 8000))
names(CarbC_P)<-c("Koprivno","mcmc","isotope")
ggplot(CarbC_P[1:8000,], aes(mcmc, fill=Sex_TimePeriod, colour=Sex_TimePeriod)) + geom_density(alpha=0.5) + theme_classic() + labs(x=C_label)
ggplot(CarbC_P[1:8000,], aes(x=Koprivno, y=mcmc, fill=Koprivno, colour=Koprivno)) + geom_bar(stat="identity") + theme_classic() + labs(x=C_label)

Output_Paired<-rbind(Nitt_P,Carb_P,CarbC_P)
tiff("paired_sex_models.tiff", unit="in", width=11, height=3.5,res=300)
ggplot(Output_Paired, aes(mcmc, fill=Koprivno, colour=Koprivno)) + geom_density(alpha=0.5) + theme_classic() + labs(x=Cc_label) + facet_wrap(.~isotope, scales="free")
dev.off()


library(coda)
Single_NColl<-mclapply(1:3, function(i){
  MCMCglmm(d13Ccoll~TimePeriod + Sex,
           rcov=~ units, 
           family="gaussian",
           data=datanew,
           nitt=nitt,
           burnin=burnin,
           thin=thin,
           prior=P_prior)
}, mc.cores=3)
gelman.plot(Single_NColl$[[1]]$Sol)
