source("Create_taxonomies.R")
source("functions.R")

###################
# Resolution  Ratios

res.ratio1<-data.frame(all.ID=c(2,4,8,16,32,64),
                      non.ID=c(16,4,2,1,1,0),
                      some.ID=c(5,4,3,2,1,0),
                      row.names=rev(c("n.p","n.c","n.o","n.f","n.g","n.s"))
)

res.ratio2<-data.frame(all.ID=c(2,4,8,16,32,64),
                       non.ID=c(5,4,3,2,1,0),
                       some.ID=c(16,4,2,1,1,0),
                       row.names=rev(c("n.p","n.c","n.o","n.f","n.g","n.s"))
)

res.ratio3<-data.frame(all.ID=c(16,4,2,1,1,0),
                       non.ID=c(5,4,3,2,1,0),
                       some.ID=c(2,4,8,16,32,64),
                       row.names=rev(c("n.p","n.c","n.o","n.f","n.g","n.s"))
)

res.ratio11<-res.ratio1/rowSums(res.ratio1)
res.ratio22<-res.ratio2/rowSums(res.ratio2)
res.ratio33<-res.ratio3/rowSums(res.ratio3)

res.ratio<-rbind(res.ratio11,res.ratio22,res.ratio33)
res.ratio$label<-rep(rev(c("n.p","n.c","n.o","n.f","n.g","n.s")),3)
res.ratio$class<-rep(c("non.ID","some.ID","all.ID"),each=6)

library(ggplot2)
library(ggtern)
ggtern(data=res.ratio, aes(x=all.ID,y=non.ID, z=some.ID,colour=class,label=label)) +
  geom_point() +
  geom_line() +
  geom_label() +
  theme_showarrows() 
  #geom_line()







loop.out<-list()

loop.out$total.comm<-data.frame()
loop.out$ref.comm<-data.frame()
loop.out$rem.comm<-data.frame()
loop.out$complet.index<-data.frame()
loop.out$diversity<-data.frame()
loop.out$ordination<-data.frame()
loop.out$whole.comm<-data.frame()

###################
# Create master taxonomy
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

###################
# Sample 500 taxa from the master taxonomy 30 times as list and DF
ref.coms<-list()
for (x in 1:30){
  temp.com<-ref.out$ref.out
  
  specnumb<-floor(rnorm(1,mean(ncol(ref.out$ref.out)/2),ncol(ref.out$ref.out)/6))
  specnumb<-max(1,specnumb)
  specnumb<-min(ncol(ref.out$ref.out),specnumb)
  
  speckeep<-sample(colnames(ref.out$ref.out),specnumb)
  
  temp.com[,!colnames(temp.com)%in%speckeep]<-0
  
  ref.coms[[x]]<-vegan::rrarefy(temp.com,500)
}

###################
# Ratios for identification rollups
res.ratio<-data.frame(all.ID=c(2,4,8,16,32,64),
                      non.ID=c(16,4,2,1,1,0),
                      some.ID=c(5,4,3,2,1,0),
                      row.names=rev(c("n.p","n.c","n.o","n.f","n.g","n.s"))
                      )

res.ratio1<-res.ratio/rowSums(res.ratio)

library(ggplot2)
library(ggtern)
ggtern(data=res.ratio1, aes(x=all.ID,y=non.ID, z=some.ID)) +
  geom_line() +
  geom_label(label=rownames(res.ratio1)) +
  theme_showarrows()

###################
# Perform taxonomic shuffle
tax.rem<-taxonomy.shuffle(ref.coms=ref.coms,
                           res.ratio=res.ratio)

colnames(tax.rem$rem.df)<-gsub("-NA","",colnames(tax.rem$rem.df))
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



###################
#  Calculate index of dataset completeness

tax.lowest.rank <-
  unlist(lapply(lapply(tax.rem$rem.taxlist, "[[", 2), tail, n = 1), use.names = T)

t1<-colSums(tax.rem$rem.df)

totala<-sum(t1)
totala.s<-sum(t1[tax.lowest.rank=="species"])/totala
totala.g<-sum(t1[tax.lowest.rank=="genus"])/totala
totala.f<-sum(t1[tax.lowest.rank=="family"])/totala
totala.o<-sum(t1[tax.lowest.rank=="order"])/totala
totala.c<-sum(t1[tax.lowest.rank=="class"])/totala
totala.p<-sum(t1[tax.lowest.rank=="phylum"])/totala


ca.score<-((totala.p*1) +
             (totala.c*2) +
             (totala.o*4) +
             (totala.f*8) +
             (totala.g*16) +
             (totala.s*32))/
  (32+16+8+4+2+1)

###################
#  compare raw datasets - metrics

rem.div.sh.dif<-vegan::diversity(tax.rem$rem.df)-vegan::diversity(tax.rem$ref.df)
dk.div.sh.dif<-vegan::diversity(dk.df)-vegan::diversity(tax.rem$ref.df)
ru.div.sh.dif<-vegan::diversity(ru.df)-vegan::diversity(tax.rem$ref.df)
rd.div.sh.dif<-vegan::diversity(rd.df)-vegan::diversity(tax.rem$ref.df)

rem.ric.dif<-vegan::specnumber(tax.rem$rem.df)-vegan::specnumber(tax.rem$ref.df)
dk.ric.dif<-vegan::specnumber(dk.df)-vegan::specnumber(tax.rem$ref.df)
ru.ric.dif<-vegan::specnumber(ru.df)-vegan::specnumber(tax.rem$ref.df)
rd.ric.dif<-vegan::specnumber(rd.df)-vegan::specnumber(tax.rem$ref.df)

rem.bray.diff<-NA
dk.bray.diff<-NA
ru.bray.diff<-NA
rd.bray.diff<-NA

for (x in 1:30){
  t1<-plyr::rbind.fill(tax.rem$rem.df[x,],tax.rem$ref.df[x,])
  t1[is.na(t1)]<-0
  rem.bray.diff[x]<-vegan::vegdist(t1)
  
  t1<-plyr::rbind.fill(dk.df[x,],tax.rem$ref.df[x,])
  t1[is.na(t1)]<-0
  dk.bray.diff[x]<-vegan::vegdist(t1)
  
  t1<-plyr::rbind.fill(ru.df[x,],tax.rem$ref.df[x,])
  t1[is.na(t1)]<-0
  ru.bray.diff[x]<-vegan::vegdist(t1)
  
  t1<-plyr::rbind.fill(rd.df[x,],tax.rem$ref.df[x,])
  t1[is.na(t1)]<-0
  rd.bray.diff[x]<-vegan::vegdist(t1)
}

###################
#  compare raw datasets - ordinations

rem.forOrd<-log10(tax.rem$rem.df+1)
ref.forOrd<-log10(tax.rem$ref.df+1)
dk.forOrd<-log10(dk.df+1)
ru.forOrd<-log10(ru.df+1)
rd.forOrd<-log10(rd.df+1)

rem.forOrd<-rem.forOrd[,apply(rem.forOrd,2,function(x)length(which(x>0)))>2]
ref.forOrd<-ref.forOrd[,apply(ref.forOrd,2,function(x)length(which(x>0)))>2]
dk.forOrd<-dk.forOrd[,apply(dk.forOrd,2,function(x)length(which(x>0)))>2]
ru.forOrd<-ru.forOrd[,apply(ru.forOrd,2,function(x)length(which(x>0)))>2]
rd.forOrd<-rd.forOrd[,apply(rd.forOrd,2,function(x)length(which(x>0)))>2]

ref.ca<-vegan::cca(ref.forOrd)
rem.ca<-vegan::cca(rem.forOrd)
dk.ca<-vegan::cca(dk.forOrd)
ru.ca<-vegan::cca(ru.forOrd)
rd.ca<-vegan::cca(rd.forOrd)

rem.proc<-vegan::procrustes(ref.ca,rem.ca,symmetric = TRUE)
dk.proc<-vegan::procrustes(ref.ca,dk.ca,symmetric = TRUE)
ru.proc<-vegan::procrustes(ref.ca,ru.ca,symmetric = TRUE)
rd.proc<-vegan::procrustes(ref.ca,rd.ca,symmetric = TRUE)

###################
#  compare raw datasets - matrix procrustes

if (F) {
  rem.forOrd<-log10(tax.rem$rem.df+1)
  ref.forOrd<-log10(tax.rem$ref.df+1)
  dk.forOrd<-log10(dk.df+1)
  ru.forOrd<-log10(ru.df+1)
  rd.forOrd<-log10(rd.df+1)
  
  rem.forOrd<-rem.forOrd[,apply(rem.forOrd,2,function(x)length(which(x>0)))>2]
  ref.forOrd<-ref.forOrd[,apply(ref.forOrd,2,function(x)length(which(x>0)))>2]
  dk.forOrd<-dk.forOrd[,apply(dk.forOrd,2,function(x)length(which(x>0)))>2]
  ru.forOrd<-ru.forOrd[,apply(ru.forOrd,2,function(x)length(which(x>0)))>2]
  rd.forOrd<-rd.forOrd[,apply(rd.forOrd,2,function(x)length(which(x>0)))>2]
  
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


