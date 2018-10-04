source("functions.R")

metrics<-lapply(output9,"[[",5)
pw.metrics<-lapply(output9,"[[",6)

ordination<-lapply(output9,"[[",7)

metrics1<-do.call(rbind,metrics)
pw.metrics1<-do.call(rbind,pw.metrics)

ordination1<-do.call(rbind,ordination)

library(ggplot2)

metric.plot<-ggplot(aes(x=Treatment,y=mean),data=metrics1)+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(Resolution~Index,scales="free_y")

pw.metrics.plot<-ggplot(aes(x=Treatment,y=Mantel_R),data=pw.metrics1)+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(Resolution~Index,scales="free_y")


metric.plot1<-ggplot(aes(x=Resolution,y=mean,fill=Treatment),
                     data=metrics1[metrics1$Treatment%in%c("dec.key","rd","ru"),])+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~Index,scales="free_y")

met.decKey<-ggplot(aes(x=Resolution,y=mean,fill=Treatment),
                   data=metrics1[metrics1$Treatment%in%unique(metrics1$Treatment)[1:10],])+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~Index,scales="free_y")


pw.metric.plot1<-ggplot(aes(x=Resolution,y=Mantel_R,fill=Treatment),
                     data=pw.metrics1[pw.metrics1$Treatment%in%c("dec.key","rd","ru"),])+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~Index,scales="free_y")

pw.met.decKey<-ggplot(aes(x=Resolution,y=Mantel_R,fill=Treatment),
                   data=pw.metrics1[pw.metrics1$Treatment%in%unique(pw.metrics1$Treatment)[2:11],])+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~Index,scales="free_y")


# metric.plot2<-beanplot::beanplot(value~as.factor(Treatment)*as.factor(Resolution),
#                                  data=metrics1[metrics1$Treatment%in%c("dec.key","rd"),],
#                                  ll=0.04,
#                                  side = "both", xlab="Treatment",
#                                  col = list("purple", c("lightblue", "black")),
#                                  axes=F)

PAord.plot<-ggplot(aes(x=Treatment,y=Scale),data=ordination1[ordination1$Index=="Presence/Absence",])+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(Resolution~Index,scales="free_y",nrow=3)

Abundord.plot<-ggplot(aes(x=Treatment,y=Scale),data=ordination1[ordination1$Index=="Abundance",])+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(Resolution~Index,scales="free_y",nrow=3)
