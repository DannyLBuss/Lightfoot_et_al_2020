###Pearson's correlations for all models
#d13Ccarb
#d15Ncoll
#d13Ccoll

#AllData_DF
DF<-Nit_DF

DF<-Nit_DF
cor.test(DF$d13Ccoll, DF$d15Ncoll, method=c("pearson"),use = "complete.obs")
cor.test(DF$d13Ccarb, DF$d15Ncoll, method=c("pearson"),use = "complete.obs")
cor.test(DF$d13Ccarb, DF$d13Ccoll, method=c("pearson"),use = "complete.obs")

#All_NoSisak
WS<-Nit_DF[!Nit_DF$name=="Sisak",]


cor.test(WS$d13Ccoll, WS$d15Ncoll, method=c("pearson"),use = "complete.obs")
cor.test(WS$d13Ccarb, WS$d15Ncoll, method=c("pearson"),use = "complete.obs")
cor.test(WS$d13Ccarb, WS$d13Ccoll, method=c("pearson"),use = "complete.obs")

#Koprivno_only
KP<-Nit_DF[Nit_DF$name=="Koprivno",]

cor.test(KP$d13Ccoll, KP$d15Ncoll, method=c("pearson"),use = "complete.obs")
cor.test(KP$d13Ccarb, KP$d15Ncoll, method=c("pearson"),use = "complete.obs")
cor.test(KP$d13Ccarb, KP$d13Ccoll, method=c("pearson"),use = "complete.obs")

#Sisak_only
SIS<-Nit_DF[Nit_DF$name=="Sisak",]

cor.test(SIS$d13Ccoll, SIS$d15Ncoll, method=c("pearson"),use = "complete.obs")
cor.test(SIS$d13Ccarb, SIS$d15Ncoll, method=c("pearson"),use = "complete.obs")
cor.test(SIS$d13Ccarb, SIS$d13Ccoll, method=c("pearson"),use = "complete.obs")

DG<-Nit_DF[Nit_DF$name=="Drinovci-Greblje",]
cor.test(DG$d13Ccoll, DG$d15Ncoll, method=c("pearson"),use = "complete.obs")

DG<-Nit_DF[Nit_DF$name=="Stenjevec",]
cor.test(DG$d13Ccoll, DG$d15Ncoll, method=c("pearson"),use = "complete.obs")

###Normality tests
shapiro.test(DF$d15Ncoll)
qqnorm(DF$d15Ncoll)
shapiro.test(DF$d13Ccoll)
shapiro.test(DF$d13Ccarb)
qqnorm(DF$d13Ccoll)
qqnorm(DF$d13Ccarb)

shapiro.test(WS$d15Ncoll)
qqnorm(WS$d15Ncoll)
shapiro.test(WS$d13Ccoll)
shapiro.test(WS$d13Ccarb)
qqnorm(WS$d13Ccoll)
qqnorm(WS$d13Ccarb)

shapiro.test(KP$d15Ncoll)
qqnorm(KP$d15Ncoll)
shapiro.test(KP$d13Ccoll)
shapiro.test(KP$d13Ccarb)
qqnorm(KP$d13Ccoll)
qqnorm(KP$d13Ccarb)

###One outlier at Koprivno
hist(DF$d15Ncoll)
hist(WS$d15Ncoll)
hist(DF$d13Ccoll)
hist(WS$d13Ccoll)
hist(DF$d13Ccarb)
hist(WS$d13Ccarb)
