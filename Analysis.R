library(tidyverse)


# Scenarios ---------------------------------------------------------------

res.ratio1<-data.frame(all.ID=c(5,15,30,60,75,100),
                       non.ID=c(70,55,35,10,0,0),
                       some.ID=c(25,30,35,30,25,0),
                       row.names=rev(c("n.p","n.c","n.o","n.f","n.g","n.s"))
)

res.ratio2<-data.frame(all.ID=c(57,65,73,81,89,100),
                       non.ID=c(13,10,7,4,1,0),
                       some.ID=c(30,25,20,15,10,0),
                       row.names=rev(c("n.p","n.c","n.o","n.f","n.g","n.s"))
)

res.ratio3<-data.frame(all.ID=c(10,28,46,64,82,100),
                       non.ID=c(20,17,14,11,8,0),
                       some.ID=c(70,55,40,25,10,0),
                       row.names=rev(c("n.p","n.c","n.o","n.f","n.g","n.s"))
)


res.ratio<-rbind(res.ratio1,res.ratio2,res.ratio3)
res.ratio$label<-rep(rev(c("n.p","n.c","n.o","n.f","n.g","n.s")),3)
res.ratio$class<-rep(c("Few Species IDs","Most Species IDs","Intermediate Species IDs"),each=6)

res.ratio<-reshape2::melt(res.ratio,id=c("label","class"))
res.ratio<-res.ratio[res.ratio$label!="n.p",]

ggplot(data=res.ratio, aes(x=label,y=value,shape=variable,colour=variable, group=variable)) +
  geom_point()+
  geom_line() +
  scale_x_discrete(limits=rev(c("n.c","n.o","n.f","n.g","n.s")))+
  facet_wrap(~class)

# Simulation Summaries ----------------------------------------------------

sim_stats<-do.call(rbind,lapply(output9,"[[",9)) %>% mutate(comm1=gsub('[[:digit:]]+',"",comm))
var_stats<-do.call(rbind,lapply(output9,"[[",10)) %>% mutate(comm1=gsub('[[:digit:]]+',"",comm))

sim_stats[,2:10]<-apply(sim_stats[,2:10],2,as.numeric) 
var_stats[,2:10]<-apply(var_stats[,2:10],2,as.numeric) 


# # Whole community -------------------------------------------------------
ggplot(sim_stats %>% 
         gather("Attribute","Value",-comm,-rep,-Resolution,-comm1) %>%
         filter(comm=="Whole", !Attribute %in% c("total","unique_taxa_names","n.p")),
       aes(x=Attribute,y=Value)) +
  geom_boxplot() +
  scale_x_discrete(limits=c("n.c","n.o","n.f","n.g","n.s")) #+
  scale_y_log10()

ggplot(sim_stats %>% 
         gather("Attribute","Value",-comm,-rep,-Resolution,-comm1) %>%
         filter(comm=="Whole", Attribute %in% c("total")),
       aes(x=Attribute,y=Value)) +
  geom_boxplot() 



# # Loval species pools and Taxonomy Censoring ----------------------------

d1<-sim_stats %>% 
  gather("Attribute","Number",-comm,-rep,-Resolution,-comm1) %>%
  filter(comm!="Whole", !Attribute %in% c("total","unique_taxa_names","n.p"))

ggplot(var_stats %>% 
         gather("Attribute","Value",-comm,-rep,-Resolution,-comm1) %>%
         filter(!Attribute %in% c("total","unique_taxa_names","n.p")) %>% 
         filter(comm1 %in% c("Sub")) %>% 
         rbind(d1),
       aes(x=Attribute,y=Value,fill=comm1)) +
  geom_boxplot() +
  scale_x_discrete(limits=c("n.c","n.o","n.f","n.g","n.s")) +
  facet_wrap(~Resolution) #+
scale_y_log10()

# # Resolving Taxonomy censoring ------------------------------------------

ggplot(var_stats %>% 
         gather("Attribute","Number",-comm,-rep,-Resolution,-comm1) %>%
         filter(!Attribute %in% c("total","unique_taxa_names","n.p")) %>% 
         rbind(d1) %>% 
         filter(comm1 %in% c("Rem","Sub","dec.key","rd","ru")) %>%
         mutate(Group=as.factor(comm1)) %>% 
         mutate(Group = factor(Group, levels=c("ru","rd","dec.key","Rem","Sub"))) %>% 
         mutate(Attribute = replace(Attribute,Attribute=="n.c","Class"),
                Attribute = replace(Attribute,Attribute=="n.o","Order"),
                Attribute = replace(Attribute,Attribute=="n.f","Family"),
                Attribute = replace(Attribute,Attribute=="n.g","Genus"),
                Attribute = replace(Attribute,Attribute=="n.s","Species"))
         ,
       aes(x=Attribute,y=Number,fill=Group)) +
  geom_boxplot(outlier.size=0) +
  scale_x_discrete(limits=c("Class","Order","Family","Genus","Species")) +
  scale_fill_discrete(labels=c("Roll-up","Roll-down","Decision Key","Raw Censored","Species Pool")) +
  facet_wrap(~Resolution,ncol=1)

# Results  ----------------------------------------------------------------

metrics<-lapply(output9,"[[",5)
pw.metrics<-lapply(output9,"[[",6)

ordination<-lapply(output9,"[[",7)

metrics1<-do.call(rbind,metrics)
pw.metrics1<-do.call(rbind,pw.metrics)

ordination1<-do.call(rbind,ordination)

metric.plot1<-ggplot(aes(x=Resolution,y=mean,fill=Treatment),
                     data=metrics1[metrics1$Treatment%in%c("dec.key","rd","ru"),])+
  geom_boxplot() +
  ylab("Mean Difference from\n True Community") +
  xlab("Taxonomic Censoring")+
  #scale_x_discrete(labels=c("","",""),name="")+
  scale_fill_discrete(labels=c("Decision Key","Roll-down","Roll-up"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept=0,linetype="dashed",col="grey")+
  facet_wrap(~Index,scales="free_y",ncol=3)

pw.metric.plot1<-ggplot(aes(x=Resolution,y=Mantel_R,fill=Treatment),
                     data=pw.metrics1[pw.metrics1$Treatment%in%c("dec.key","rd","ru"),])+
  geom_boxplot() +
  ylab("Mantel R") +
  xlab("Taxonomic Censoring")+
  #scale_x_discrete(labels=c("","",""),name="")+
  scale_fill_discrete(labels=c("Decision Key","Roll-down","Roll-up"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Index,scales="fixed",ncol=3)

Ord.plot<-ggplot(aes(x=Resolution,y=SS,fill=Treatment),data=ordination1[ordination1$Treatment%in%c("dec.key","rd","ru"),])+
  geom_boxplot() +
  ylab("Procrustes Sums of Squares") +
  xlab("Taxonomic Censoring")+
  #scale_x_discrete(labels=c("","",""),name="")+
  scale_fill_discrete(labels=c("Decision Key","Roll-down","Roll-up"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Index,scales="fixed",ncol=3)







metric.plot<-ggplot(aes(x=Treatment,y=mean),data=metrics1)+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(Resolution~Index,scales="free_y")

pw.metrics.plot<-ggplot(aes(x=Treatment,y=Mantel_R),data=pw.metrics1)+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(Resolution~Index,scales="free_y")

met.decKey<-ggplot(aes(x=Resolution,y=mean,fill=Treatment),
                   data=metrics1[metrics1$Treatment%in%unique(metrics1$Treatment)[1:10],])+
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
