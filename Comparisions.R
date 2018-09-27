rm(list=ls())

source("Create_taxonomies.R")
source("functions.R")

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
loop.out$ordination<-data.frame()
loop.out$whole.comm<-data.frame()

out<-list()
out$few.sp.ID<-loop.out
out$most.sp.ID<-loop.out
out$inter.sp.ID<-loop.out

for (i in 1:3){
  res.ratio<-list(res.ratio1,res.ratio2,res.ratio3)[[i]]
  res<-c("Few Species IDs","Most Species IDs","Intermediate Species IDs")[i]
  print(i)
  for (n in 1:100){
    
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
        
        if (sum(ref.out$ref.out)<500){
          max.rar<-300
        }
        
        temp.com<-ref.out$ref.out
        
        specnumb<-floor(rnorm(1,mean(ncol(ref.out$ref.out)/2),ncol(ref.out$ref.out)/6))
        specnumb<-max(1,specnumb)
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
    
    #out[[i]]$ref.comm[[n]]<-ref.coms
    
    ###################
    # Ratios for identification rollups
    #res.ratio<-data.frame(all.ID=c(2,4,8,16,32,64),
    #                      non.ID=c(16,4,2,1,1,0),
    #                      some.ID=c(5,4,3,2,1,0),
    #                      row.names=rev(c("n.p","n.c","n.o","n.f","n.g","n.s"))
    #                      )
    
    #res.ratio1<-res.ratio/rowSums(res.ratio)
    
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
                         Criteria3.percent=0.1, #critical limit for criteria 3
                         Criteria5a.percent=0.5, #critical limit for criteria 5a
                         Criteria5b.numb=2 #critical limit for criteria 5b
    )
    dk.higherP3.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                         taxa.names=tax.rem$rem.taxlist,
                         roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                         assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                         Criteria3.percent=0.3, #critical limit for criteria 3
                         Criteria5a.percent=0.5, #critical limit for criteria 5a
                         Criteria5b.numb=2 #critical limit for criteria 5b
    )
    dk.lowerP5.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                         taxa.names=tax.rem$rem.taxlist,
                         roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                         assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                         Criteria3.percent=0.2, #critical limit for criteria 3
                         Criteria5a.percent=0.4, #critical limit for criteria 5a
                         Criteria5b.numb=2 #critical limit for criteria 5b
    )
    dk.higherP5.df<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                         taxa.names=tax.rem$rem.taxlist,
                         roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                         assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                         Criteria3.percent=0.2, #critical limit for criteria 3
                         Criteria5a.percent=0.6, #critical limit for criteria 5a
                         Criteria5b.numb=2 #critical limit for criteria 5b
    )
    
    ###################
    #  Calculate index of dataset completeness
    
    # tax.lowest.rank <-
    #   unlist(lapply(lapply(tax.rem$rem.taxlist, "[[", 2), tail, n = 1), use.names = T)
    # 
    # t1<-colSums(tax.rem$rem.df)
    # 
    # totala<-sum(t1)
    # totala.s<-sum(t1[tax.lowest.rank=="species"])/totala
    # totala.g<-sum(t1[tax.lowest.rank=="genus"])/totala
    # totala.f<-sum(t1[tax.lowest.rank=="family"])/totala
    # totala.o<-sum(t1[tax.lowest.rank=="order"])/totala
    # totala.c<-sum(t1[tax.lowest.rank=="class"])/totala
    # totala.p<-sum(t1[tax.lowest.rank=="phylum"])/totala
    # 
    # ca.score<-((totala.p*1) +
    #              (totala.c*2) +
    #              (totala.o*4) +
    #              (totala.f*8) +
    #              (totala.g*16) +
    #              (totala.s*32))/
    #   (32+16+8+4+2+1)
    
    
    ###################
    #  compare raw datasets - metrics
    
    variants<-c("raw","dec.key","dec.key.noRoll","dec.key.noAssign",
                "dec.key.lowerP3","dec.key.higherP3","dec.key.lowerP5",
                "dec.key.higherP5","ru","rd")
    
    div.out<-data.frame(matrix(nrow=30,ncol=30+3))
    colnames(div.out)[1:3]<-c("Resolution","Treatment","Index")
    
    div.out$Index<-c(rep(c("Shannon","Richness","Bray-Curtis"),each=10)) 
    div.out$Treatment<-c(rep(variants,times=3)) 
    div.out$Resolution<-res
    
    ord.out<-data.frame(matrix(nrow=20,ncol=3+3))
    colnames(ord.out)[1:3]<-c("Resolution","Treatment","Index")
    
    ord.out$Index<-c(rep(c("Abundance","Presence/Absence"),each=10)) 
    ord.out$Treatment<-c(rep(variants,times=2)) 
    ord.out$Resolution<-res
    colnames(ord.out)[4:6]<-c("SS","Scale","pval")
    
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
      if (z=="ru"){
        test.df<-ru.df
      }
      if (z=="rd"){
        test.df<-rd.df
      }
      
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
      
    }
    
    out[[i]]$diversity<-rbind(out[[i]]$diversity,div.out)
    out[[i]]$ordination<-rbind(out[[i]]$ordination,ord.out)
    
    
    ###################
    #  compare raw datasets - matrix procrustes
    
    if (F) {
      rem.forOrd<-log10(tax.rem$rem.df+1)
      ref.forOrd<-log10(tax.rem$ref.df+1)
      dk.forOrd<-log10(dk.df+1)
      ru.forOrd<-log10(ru.df+1)
      rd.forOrd<-log10(rd.df+1)
      
      rem.forOrd<-rem.forOrd[,apply(rem.forOrd,2,function(x)length(which(x>0)))>1]
      ref.forOrd<-ref.forOrd[,apply(ref.forOrd,2,function(x)length(which(x>0)))>1]
      dk.forOrd<-dk.forOrd[,apply(dk.forOrd,2,function(x)length(which(x>0)))>1]
      ru.forOrd<-ru.forOrd[,apply(ru.forOrd,2,function(x)length(which(x>0)))>1]
      rd.forOrd<-rd.forOrd[,apply(rd.forOrd,2,function(x)length(which(x>0)))>1]
      
      rem.forOrd<-plyr::rbind.fill(rem.forOrd,ref.forOrd)
      rem.forOrd[is.na(rem.forOrd)]<-0
      
      dk.forOrd<-plyr::rbind.fill(dk.forOrd,ref.forOrd)
      dk.forOrd[is.na(dk.forOrd)]<-0
      
      ru.forOrd<-plyr::rbind.fill(ru.forOrd,ref.forOrd)
      ru.forOrd[is.na(ru.forOrd)]<-0
      
      rd.forOrd<-plyr::rbind.fill(rd.forOrd,ref.forOrd)
      rd.forOrd[is.na(rd.forOrd)]<-0
      
      rem.proc.m<-vegan::procrustes(rem.forOrd[31:60,],rem.forOrd[1:30,],symmetric = TRUE)
      dk.proc.m<-vegan::procrustes(dk.forOrd[31:60,],dk.forOrd[1:30,],symmetric = TRUE)
      ru.proc.m<-vegan::procrustes(ru.forOrd[31:60,],ru.forOrd[1:30,],symmetric = TRUE)
      rd.proc.m<-vegan::procrustes(rd.forOrd[31:60,],rd.forOrd[1:30,],symmetric = TRUE)
    }
    print(n)
  }
}



