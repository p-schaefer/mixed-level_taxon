source("Create_taxonomies.R")
source("functions.R")

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
  ref.coms[[x]]<-vegan::rrarefy(ref.out,500)
}

ref.df<-data.frame(matrix(unlist(ref.coms), nrow=ncol(ref.out), byrow=F))
rownames(ref.df)<-colnames(ref.out)
ref.df<-data.frame(t(ref.df))

###################
# Ratios for identification rollups
res.ratio<-data.frame(all.ID=c(2,4,8,16,32,64),
                      non.ID=c(16,4,2,1,0,0),
                      some.ID=c(5,4,3,2,1,0),
                      row.names=rev(c("n.p","n.c","n.o","n.f","n.g","n.s"))
                      )

res.ratio1<-res.ratio/rowSums(res.ratio)

library(ggplot2)
library(ggtern)
ggtern(data=res.ratio1, aes(x=all.ID,y=non.ID, z=some.ID)) +
  geom_point() +
  geom_label(label=rownames(res.ratio1))

###################
# Perform taxonomic shuffle
tax.rem<-taxonomy.shuffle(ref.coms=ref.coms,
                           res.ratio=res.ratio)

###################
# Perform taxonomic roll ups/downs
taxroll.dk<-benth.taxroll(taxa=tax.rem$for.taxroll, 
                          taxa.names=tax.rem$rem.taxlist,
                          roll.down=F, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                          assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                          Criteria3.percent=0.2, #critical limit for criteria 3
                          Criteria5a.percent=0.5, #critical limit for criteria 5a
                          Criteria5b.numb=2 #critical limit for criteria 5b
)

taxroll.up<-benth.rollup(taxa=tax.rem$for.taxroll, 
                         taxa.names=tax.rem$rem.taxlist,
                         roll.thresh=0.2)

taxroll.down<-benth.rolldown(taxa=tax.rem$for.taxroll, 
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
#  compare raw datasets

whole.cca<-vegan::cca(log10(ref.df+1))
rem.cca<-vegan::cca(log10(rem.df+1))
proc<-vegan::procrustes(whole.cca,rem.cca)
vegan::protest(whole.cca,rem.cca)
