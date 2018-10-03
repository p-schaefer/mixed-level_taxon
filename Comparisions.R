rm(list=ls())

source("Create_taxonomies.R")
source("functions.R")

library(dplyr)

set.seed(12345)
###################
# Resolution  Ratios

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

###################
# Create plot resolution ratios

# res.ratio11<-res.ratio1/rowSums(res.ratio1)
# res.ratio22<-res.ratio2/rowSums(res.ratio2)
# res.ratio33<-res.ratio3/rowSums(res.ratio3)

res.ratio<-rbind(res.ratio1,res.ratio2,res.ratio3)
res.ratio$label<-rep(rev(c("n.p","n.c","n.o","n.f","n.g","n.s")),3)
res.ratio$class<-rep(c("Few Species IDs","Most Species IDs","Intermediate Species IDs"),each=6)

res.ratio<-reshape2::melt(res.ratio,id=c("label","class"))

res.ratio<-res.ratio[res.ratio$label!="n.p",]

if (F){
  library(ggplot2)
  ggplot(data=res.ratio, aes(x=label,y=value,shape=variable,colour=variable, group=variable)) +
    geom_point()+
    geom_line() +
    scale_x_discrete(limits=rev(c("n.c","n.o","n.f","n.g","n.s")))+
    facet_wrap(~class)
  #library(ggtern)
  # ggtern(data=res.ratio, aes(x=all.ID,y=non.ID, z=some.ID,colour=class,label=label)) +
  #   geom_point() +
  #   geom_line() +
  #   geom_label() +
  #   theme_showarrows() 
  #geom_line()
}


###################
# loop setup

loop.out<-list()

loop.out$total.comm<-list()
loop.out$ref.comm<-list()
loop.out$rem.comm<-list()
loop.out$complet.index<-data.frame()
loop.out$diversity<-data.frame()
loop.out$pwise.diversity<-data.frame()
loop.out$ordination<-data.frame()
loop.out$whole.comm<-data.frame()
loop.out$sim.stats<-data.frame()
loop.out$var.stats<-data.frame()

out<-list()
out$few.sp.ID<-loop.out
out$most.sp.ID<-loop.out
out$inter.sp.ID<-loop.out

pb = txtProgressBar(min = 0, max = 300, initial = 0,style=3) 

for (n in 1:100){
  for (i in 1:3){
    
    res.ratio<-list(res.ratio1,res.ratio2,res.ratio3)[[i]]
    res<-c("Few Species IDs","Most Species IDs","Intermediate Species IDs")[i]
    #print(i)
    setTxtProgressBar(pb,n*i,label=res)
    
    ###################
    # Create master taxonomy
    repeat{
      ref.out<-master.taxonomy(p.n=5, # number of phyla to start with
                               p.s=2, # Class size parameter - nbinom
                               p.m=2, # Class mu parameter - nbinom
                               c.s=2, # Order size parameter - nbinom
                               c.m=2, # Order size parameter - nbinom
                               o.s=2, # Family size parameter - nbinom
                               o.m=2, # Family size parameter - nbinom
                               f.s=2, # Genus size parameter - nbinom
                               f.m=2, # Genus size parameter - nbinom
                               g.s=2, # Species size parameter - nbinom
                               g.m=2, # Species size parameter - nbinom
                               s.mean=5, # log-normal mean for species abundances
                               s.sd=2 # log-normal sd for species abundances
      )
      if (rowSums(ref.out$ref.out)>600) {break}
    }
    
    #out[[i]]$total.comm[[n]]<-ref.out
    ###################
    # Sample 500 taxa from the master taxonomy 30 times as list and DF
    ref.coms<-list()
    for (x in 1:30){
      repeat{
        max.rar<-500
        
        # if (sum(ref.out$ref.out)<500){
        #   max.rar<-300
        # }
        
        temp.com<-ref.out$ref.out
        
        #specnumb<-floor(rnbinom(1,size=ncol(ref.out$ref.out)/2,mu=ncol(ref.out$ref.out)/2))
        specnumb<-floor(runif(1,ncol(ref.out$ref.out)*(1/4),ncol(ref.out$ref.out)*(3/4)))
        specnumb<-max(5,specnumb)
        specnumb<-min(ncol(ref.out$ref.out),specnumb)
        
        speckeep<-sample(colnames(ref.out$ref.out),specnumb)
        
        temp.com[,!colnames(temp.com)%in%speckeep]<-0
        
        temp.com<-vegan::rrarefy(temp.com,max.rar)
        
        if (rowSums(temp.com)>100 & length(which(temp.com>0))>4){
          break
        }
      }
      ref.coms[[x]]<-temp.com
    }
    

    ###################
    # Perform taxonomic shuffle
    tax.rem<-taxonomy.shuffle(ref.coms=ref.coms,
                              res.ratio=res.ratio)
    
    colnames(tax.rem$rem.df)<-gsub("-NA","",colnames(tax.rem$rem.df))
    
    #out[[i]]$rem.comm[[n]]<-tax.rem
    
    ###################
    # Perform taxonomic roll ups/downs
    dk.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                         taxa.names=tax.rem$rem.taxlist,
                         roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                         assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                         Criteria3.percent=0.2, #critical limit for criteria 3
                         Criteria5a.percent=0.5, #critical limit for criteria 5a
                         Criteria5b.numb=2 #critical limit for criteria 5b
    )
    
    ru.df<-benth.rollup(taxa=tax.rem$for.taxroll, 
                        taxa.names=tax.rem$rem.taxlist,
                        roll.thresh=0.2)
    
    rd.df<-benth.rolldown(taxa=tax.rem$for.taxroll, 
                          taxa.names=tax.rem$rem.taxlist)
    
    
    dk.noRoll.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                         taxa.names=tax.rem$rem.taxlist,
                         roll.down=F, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                         assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                         Criteria3.percent=0.2, #critical limit for criteria 3
                         Criteria5a.percent=0.5, #critical limit for criteria 5a
                         Criteria5b.numb=2 #critical limit for criteria 5b
    )
    dk.noAssign.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                         taxa.names=tax.rem$rem.taxlist,
                         roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                         assign.undet=F, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                         Criteria3.percent=0.2, #critical limit for criteria 3
                         Criteria5a.percent=0.5, #critical limit for criteria 5a
                         Criteria5b.numb=2 #critical limit for criteria 5b
    )
    dk.lowerP3.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                         taxa.names=tax.rem$rem.taxlist,
                         roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                         assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                         Criteria3.percent=0.05, #critical limit for criteria 3
                         Criteria5a.percent=0.5, #critical limit for criteria 5a
                         Criteria5b.numb=2 #critical limit for criteria 5b
    )
    dk.higherP3.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                         taxa.names=tax.rem$rem.taxlist,
                         roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                         assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                         Criteria3.percent=0.5, #critical limit for criteria 3
                         Criteria5a.percent=0.5, #critical limit for criteria 5a
                         Criteria5b.numb=2 #critical limit for criteria 5b
    )
    dk.lowerP5.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                         taxa.names=tax.rem$rem.taxlist,
                         roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                         assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                         Criteria3.percent=0.2, #critical limit for criteria 3
                         Criteria5a.percent=0.3, #critical limit for criteria 5a
                         Criteria5b.numb=2 #critical limit for criteria 5b
    )
    dk.higherP5.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                         taxa.names=tax.rem$rem.taxlist,
                         roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                         assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                         Criteria3.percent=0.2, #critical limit for criteria 3
                         Criteria5a.percent=0.7, #critical limit for criteria 5a
                         Criteria5b.numb=2 #critical limit for criteria 5b
    )
    dk.higherP52.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                                  taxa.names=tax.rem$rem.taxlist,
                                  roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                                  assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                                  Criteria3.percent=0.2, #critical limit for criteria 3
                                  Criteria5a.percent=0.9, #critical limit for criteria 5a
                                  Criteria5b.numb=2 #critical limit for criteria 5b
    )
    
    dk.higherP5b.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                                  taxa.names=tax.rem$rem.taxlist,
                                  roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                                  assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                                  Criteria3.percent=0.2, #critical limit for criteria 3
                                  Criteria5a.percent=0.5, #critical limit for criteria 5a
                                  Criteria5b.numb=3 #critical limit for criteria 5b
    )
    dk.higherP5b2.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                                   taxa.names=tax.rem$rem.taxlist,
                                   roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                                   assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                                   Criteria3.percent=0.2, #critical limit for criteria 3
                                   Criteria5a.percent=0.5, #critical limit for criteria 5a
                                   Criteria5b.numb=4 #critical limit for criteria 5b
    )
    
    ###################
    #  compare raw datasets - metrics
    
    variants<-c("raw","dec.key","dec.key.noRoll","dec.key.noAssign",
                "dec.key.lowerP3","dec.key.higherP3","dec.key.lowerP5",
                "dec.key.higherP5","dec.key.higherP52","dec.key.higherP5b","dec.key.higherP5b2","ru","rd")
    
    div.out<-data.frame(matrix(nrow=length(variants)*3,ncol=30+3))
    colnames(div.out)[1:3]<-c("Resolution","Treatment","Index")
    
    div.out$Index<-c(rep(c("Shannon","Richness","Bray-Curtis"),each=length(variants))) 
    div.out$Treatment<-c(rep(variants,times=3)) 
    div.out$Resolution<-res
    
    pwise.diff<-data.frame(matrix(nrow=length(variants)*3,ncol=3+2))
    colnames(pwise.diff)[1:3]<-c("Resolution","Treatment","Index")
    
    pwise.diff$Index<-c(rep(c("Shannon","Richness","Bray-Curtis"),each=length(variants))) 
    pwise.diff$Treatment<-c(rep(variants,times=3)) 
    pwise.diff$Resolution<-res
    colnames(pwise.diff)[4:5]<-c("Mantel_R","pval")
    
    ord.out<-data.frame(matrix(nrow=length(variants)*2,ncol=3+3))
    colnames(ord.out)[1:3]<-c("Resolution","Treatment","Index")
    
    ord.out$Index<-c(rep(c("Abundance","Presence/Absence"),each=length(variants))) 
    ord.out$Treatment<-c(rep(variants,times=2)) 
    ord.out$Resolution<-res
    colnames(ord.out)[4:6]<-c("SS","Scale","pval")
    
    var.stats.out<-data.frame()
    
    for (z in variants) {
      if (z=="raw"){
        test.df<-tax.rem$rem.df
      }
      if (z=="dec.key"){
        test.df<-dk.df
      }
      if (z=="dec.key.noRoll"){
        test.df<-dk.noRoll.df
      }
      if (z=="dec.key.noAssign"){
        test.df<-dk.noAssign.df
      }
      if (z=="dec.key.lowerP3"){
        test.df<-dk.lowerP3.df
      }
      if (z=="dec.key.higherP3"){
        test.df<-dk.higherP3.df
      }
      if (z=="dec.key.lowerP5"){
        test.df<-dk.lowerP5.df
      }
      if (z=="dec.key.higherP5"){
        test.df<-dk.higherP5.df
      }
      if (z=="dec.key.higherP52"){
        test.df<-dk.higherP52.df
      }
      if (z=="dec.key.higherP5b"){
        test.df<-dk.higherP5b.df
      }
      if (z=="dec.key.higherP5b2"){
        test.df<-dk.higherP5b2.df
      }
      if (z=="ru"){
        test.df<-ru.df
      }
      if (z=="rd"){
        test.df<-rd.df
      }
      
      #Pairwise differences
      shan.p<-vegan::mantel(dist(vegan::diversity(test.df)),dist(vegan::diversity(tax.rem$ref.df)))
      pwise.diff$Mantel_R[div.out$Index=="Shannon" & div.out$Treatment==z]<-shan.p$statistic
      pwise.diff$pval[div.out$Index=="Shannon" & div.out$Treatment==z]<-shan.p$signif
      
      rich.p<-vegan::mantel(dist(vegan::specnumber(test.df)),dist(vegan::specnumber(tax.rem$ref.df)))
      pwise.diff$Mantel_R[div.out$Index=="Richness" & div.out$Treatment==z]<-rich.p$statistic
      pwise.diff$pval[div.out$Index=="Richness" & div.out$Treatment==z]<-rich.p$signif
      
      bray.p<-vegan::mantel(vegan::vegdist(test.df),vegan::vegdist(tax.rem$ref.df))
      pwise.diff$Mantel_R[div.out$Index=="Bray-Curtis" & div.out$Treatment==z]<-bray.p$statistic
      pwise.diff$pval[div.out$Index=="Bray-Curtis" & div.out$Treatment==z]<-bray.p$signif

      #Diversity + Richness
      div.out[div.out$Index=="Shannon" & div.out$Treatment==z,-c(1:3)]<-vegan::diversity(test.df)-vegan::diversity(tax.rem$ref.df)
      div.out[div.out$Index=="Richness" & div.out$Treatment==z,-c(1:3)]<-vegan::specnumber(test.df)-vegan::specnumber(tax.rem$ref.df)
      
      #Bray-Curtis
      bray.diff<-NA
      for (x in 1:30){
        t1<-plyr::rbind.fill(test.df[x,],tax.rem$ref.df[x,])
        t1[is.na(t1)]<-0
        bray.diff[x]<-vegan::vegdist(t1)
      }
      
      div.out[div.out$Index=="Bray-Curtis" & div.out$Treatment==z,-c(1:3)]<-bray.diff
      
      #Abundance Ordination
      rem.forOrd<-log10(test.df+1)
      ref.forOrd<-log10(tax.rem$ref.df+1)

      rem.forOrd<-rem.forOrd[,apply(rem.forOrd,2,function(x)length(which(x>0)))>1]
      ref.forOrd<-ref.forOrd[,apply(ref.forOrd,2,function(x)length(which(x>0)))>1]

      ref.ca<-vegan::cca(ref.forOrd)
      rem.ca<-vegan::cca(rem.forOrd)

      rem.proc<-vegan::protest(ref.ca,rem.ca)

      ord.out[ord.out$Index=="Abundance" & ord.out$Treatment==z,-c(1:3)]<-c(rem.proc$ss,rem.proc$scale,rem.proc$signif)
      
      #P/A Ordination
      rem.forOrd<-test.df
      rem.forOrd[rem.forOrd>0]<-1
      ref.forOrd<-tax.rem$ref.df
      ref.forOrd[ref.forOrd>0]<-1

      rem.forOrd<-rem.forOrd[,apply(rem.forOrd,2,function(x)length(which(x>0)))>1]
      ref.forOrd<-ref.forOrd[,apply(ref.forOrd,2,function(x)length(which(x>0)))>1]

      ref.ca<-vegan::cca(ref.forOrd)
      rem.ca<-vegan::cca(rem.forOrd)

      rem.proc<-vegan::protest(ref.ca,rem.ca,symmetric = TRUE)

      ord.out[ord.out$Index=="Presence/Absence" & ord.out$Treatment==z,-c(1:3)]<-c(rem.proc$ss,rem.proc$scale,rem.proc$signif)
      
      
      ############################
      # Stats on variants

      var.stats<-data.frame(plyr::ldply(strsplit(split="-",x=colnames(test.df)), rbind),stringsAsFactors = F)
      
      var.stats<-cbind(do.call(rbind,replicate(30,var.stats,simplify=F)),c(t(test.df)))
      
      colnames(var.stats)<-c("n.p",
                             "n.c",
                             "n.o",
                             "n.f",
                             "n.g",
                             "n.s",
                             "count")
      
      var.stats$comm<-rep(paste0(z,1:30),each=ncol(test.df))
      
      var.stats<-var.stats %>% 
        filter(count>0)%>%
        group_by(comm) %>%
        arrange(desc(n.s)) %>%
        dplyr::summarize(n.p=length(unique(n.p)),
                         n.c=length(unique(n.c)),
                         n.o=length(unique(n.o)),
                         n.f=length(unique(n.f)),
                         n.g=length(unique(n.g)),
                         n.s=length(unique(n.s)),
                         unique_taxa_names=length(count),
                         total=sum(count))
      var.stats$rep<-n
      
      var.stats.out<-rbind(var.stats.out,var.stats)
      
      rm(test.df)
    }
    
    ## Simulation stats
    
    temp1<-data.frame(matrix((unlist(lapply(ref.coms,function(x) strsplit(split="-",x=colnames(x))))),ncol=6,byrow=T),stringsAsFactors = F)
    colnames(temp1)<-c("n.p",
                       "n.c",
                       "n.o",
                       "n.f",
                       "n.g",
                       "n.s")
    
    temp1$count<-unlist(ref.coms)
    temp1$comm<-rep(paste0("Sub",1:30),each=as.numeric(lapply(ref.coms,ncol)[1]))
    
    temp2<-do.call(rbind,tax.rem$rem.coms)
    temp2$comm<-rep(paste0("Rem",1:30),lapply(tax.rem$rem.coms,nrow))
    temp2<-temp2 %>% select(n.p,n.c,n.o,n.f,n.g,n.s,count,comm)
    
    temp3<-rbind(temp1,temp2)
    
    temp3<-temp3 %>% 
      filter(count>0)%>%
      group_by(comm) %>%
      arrange(desc(n.s)) %>%
      dplyr::summarize(n.p=length(unique(n.p)),
                       n.c=length(unique(n.c)),
                       n.o=length(unique(n.o)),
                       n.f=length(unique(n.f)),
                       n.g=length(unique(n.g)),
                       n.s=length(unique(n.s)),
                       unique_taxa_names=length(count),
                       total=sum(count))
    temp3[nrow(temp3)+1,]<-c("Whole",apply(ref.out$taxa[,1:6],2,function(x) length(unique(x))),
                             length(ref.out$taxa$n.s),
                             sum(ref.out$taxa$count))
    temp3$rep<-n
    

    div.out$rep<-n
    ord.out$rep<-n
    
    div.out <- div.out %>% 
      tidyr::gather("sample","difference",-c(1:3)) %>%
      group_by(Treatment,Index) %>%
      summarize(min=min(difference),
                mean=mean(difference),
                max=max(difference),
                sd=sd(difference)) %>%
      arrange(Index) %>%
      as.data.frame()
    
    
    #Paste all values to output
    temp3$Resolution<-res
    var.stats.out$Resolution<-res
    div.out$Resolution<-res
    div.out$rep<-n
    pwise.diff$rep<-n
    
    
    out[[i]]$sim.stats<-rbind(out[[i]]$sim.stats,temp3)
    
    out[[i]]$var.stats<-rbind(out[[i]]$var.stats,var.stats.out)
    
    out[[i]]$pwise.diversity<-rbind(out[[i]]$pwise.diversity,pwise.diff)
    out[[i]]$diversity<-rbind(out[[i]]$diversity,div.out)
    out[[i]]$ordination<-rbind(out[[i]]$ordination,ord.out)
    
    #print(n)
  }
  
}



