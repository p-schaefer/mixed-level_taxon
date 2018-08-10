#set.seed(12345)

loop.mem.sum<-function(loop.mem,loop.list){
  start<-1
  out<-vector()
  for (i in loop.list){
    end<-start+i-1
    out<-c(out,sum(loop.mem[start:end]))
    start<-end+1
  }
  return(out)
}

n.p<-rnbinom(5,size=2,mu=2)
n.p[n.p==0]<-1
n.c<-rnbinom(sum(n.p),size=2,mu=2)
n.c[n.c==0]<-1
n.o<-rnbinom(sum(n.c),size=2,mu=2)
n.o[n.o==0]<-1
n.f<-rnbinom(sum(n.o),size=2,mu=2)
n.f[n.f==0]<-1
n.g<-rnbinom(sum(n.f),size=2,mu=2)
n.g[n.g==0]<-1
n.s<-rep(1,sum(n.g))

taxa<-data.frame(matrix(nrow=sum(n.s),ncol=6))
colnames(taxa)<-c("n.p","n.c","n.o","n.f","n.g","n.s")

taxa.list<-list(n.p=n.p,
                n.c=n.c,
                n.o=n.o,
                n.f=n.f,
                n.g=n.g,
                n.s=n.s
)

for (i in names(rev(taxa.list))){
  if (i=="n.s"){
    taxa[,i]<-paste0("s",1:sum(n.s))
    next
  }
  if (i=="n.g"){
    loop.list<-taxa.list[[i]]
    taxa[,i]<-paste0(gsub("n.","",i),rep(1:length(loop.list),times=loop.list))
    loop.mem<-loop.list
    next
  }
  loop.list<-taxa.list[[i]]
  rep.times<-loop.mem.sum(loop.mem,loop.list)
  taxa[,i]<-paste0(gsub("n.","",i),rep(1:length(loop.list),times=rep.times))
  loop.mem<-rep.times
}

taxa$count<-ceiling(rlnorm(nrow(taxa),5,2))

####################
#  remove taxa

res.ratio<-data.frame(all.ID=c(1,2,3,4,5,5),
                      non.ID=c(20,10,5,4,3,0),
                      some.ID=c(1,1,1,1,1,0)
)
rownames(res.ratio)<-rev(c("n.p","n.c","n.o","n.f","n.g","n.s"))

checked<-c(F,F,F,F,F,F)
names(checked)<-c("n.p","n.c","n.o","n.f","n.g","n.s")

res.ratio<-t(apply(res.ratio,1,cumsum))

taxa.whole<-taxa

for (i in names(rev(taxa.list))) {
  if (i=="n.s"){
    for (n in unique(taxa[,i])) {
      decision<-round(runif(1,min=1,max=sum(res.ratio[i,])))
      if (decision<=res.ratio[i,"all.ID"]){
        next
      } else if (decision<=res.ratio[i,"non.ID"]){
        new.tax<-taxa[taxa[,i]==n,]
        new.tax[,i]<-"NA"
        taxa<-rbind(taxa,new.tax)
        taxa<-taxa[taxa[,i]!=n,]
      } else { #some.IF
        numb.kept<-ceiling(runif(1,0.1,0.99)*taxa$count[taxa[,i]==n])
        numb.upped<-taxa$count[taxa[,i]==n]-numb.kept
        new.tax<-taxa[taxa[,i]==n,]
        new.tax[,i]<-"NA"
        new.tax$count<-numb.upped
        taxa$count[taxa[,i]==n]<-numb.kept
        taxa<-rbind(taxa,new.tax)
      }
    }
    taxa$full.taxon<-apply(taxa[,colnames(taxa)%in%c("n.p","n.c","n.o","n.f","n.g","n.s")],
                           1,paste0,collapse="-")
    taxa1<-aggregate(count~full.taxon,data=taxa,sum)
    taxa<-cbind(taxa1,do.call(rbind, strsplit(taxa1$full.taxon,"-")))
    colnames(taxa)[3:8]<-c("n.p","n.c","n.o","n.f","n.g","n.s")
    taxa[3:8]<-apply(taxa[3:8],2,as.character)
    checked[i]<-T
  } else {
    for (n in unique(
      taxa$full.taxon[taxa[,names(checked[which.max(checked)])]=="NA"]
    )
    ) {
      decision<-round(runif(1,min=1,max=res.ratio[i,"some.ID"]))
      if (decision<=res.ratio[i,"all.ID"]){
        next
      } else if (decision<=res.ratio[i,"non.ID"]){
        new.tax<-taxa[taxa$full.taxon==n,]
        new.tax[,i]<-"NA"
        new.tax$full.taxon<-paste0(new.tax[,colnames(new.tax)%in%c("n.p","n.c","n.o","n.f","n.g","n.s")], collapse="-")
        taxa<-rbind(taxa,new.tax)
        taxa<-taxa[taxa$full.taxon!=n,]
        next
      } else { #some.IF
        numb.kept<-ceiling(runif(1,0.1,0.99)*taxa$count[taxa$full.taxon==n])
        numb.upped<-taxa$count[taxa$full.taxon==n]-numb.kept
        new.tax<-taxa[taxa$full.taxon==n,]
        new.tax[,i]<-"NA"
        new.tax$count<-numb.upped
        new.tax$full.taxon<-paste0(new.tax[,colnames(new.tax)%in%c("n.p","n.c","n.o","n.f","n.g","n.s")], collapse="-")
        taxa$count[taxa$full.taxon==n]<-numb.kept
        taxa<-rbind(taxa,new.tax)
        next
      }
    }
    taxa$full.taxon<-apply(taxa[,colnames(taxa)%in%c("n.p","n.c","n.o","n.f","n.g","n.s")],
                           1,paste0,collapse="-")
    taxa1<-aggregate(count~full.taxon,data=taxa,sum)
    taxa<-cbind(taxa1,do.call(rbind, strsplit(taxa1$full.taxon,"-")))
    colnames(taxa)[3:8]<-c("n.p","n.c","n.o","n.f","n.g","n.s")
    taxa[3:8]<-apply(taxa[3:8],2,as.character)
    checked[i]<-T
  }
}

#sum(taxa$count)
#sum(taxa.whole$count)

taxa.mod<-taxa

###################
#  Calculate index of dataset completeness

total<-nrow(taxa)
total.s<-length(taxa$n.s[taxa$n.s!="NA"])/total
total.g<-length(taxa$n.g[taxa$n.g!="NA"])/total
total.f<-length(taxa$n.f[taxa$n.f!="NA"])/total
total.o<-length(taxa$n.o[taxa$n.o!="NA"])/total
total.c<-length(taxa$n.c[taxa$n.c!="NA"])/total
total.p<-length(taxa$n.p[taxa$n.p!="NA"])/total

c.score<-((total.p*1) +
            (total.c*2) +
            (total.o*4) +
            (total.f*8) +
            (total.g*16) +
            (total.s*32))/
  (32+16+8+4+2+1)

totala<-sum(taxa$count)
totala.s<-sum(taxa$count[taxa$n.g!="NA" & taxa$n.s!="NA"])/totala
totala.g<-(sum(taxa$count[taxa$n.f!="NA" & taxa$n.g!="NA"])/totala)-totala.s
totala.f<-(sum(taxa$count[taxa$n.o!="NA" & taxa$n.f!="NA"])/totala)-(totala.s+totala.g)
totala.o<-(sum(taxa$count[taxa$n.c!="NA" & taxa$n.o!="NA"])/totala)-(totala.s+totala.g+totala.f)
totala.c<-(sum(taxa$count[taxa$n.p!="NA" & taxa$n.c!="NA"])/totala)-(totala.s+totala.g+totala.f+totala.o)
totala.p<-(sum(taxa$count[taxa$n.p!="NA"])/totala)-(totala.s+totala.g+totala.f+totala.o+totala.c)


ca.score<-((totala.p*1) +
             (totala.c*2) +
             (totala.o*4) +
             (totala.f*8) +
             (totala.g*16) +
             (totala.s*32))/
  (32+16+8+4+2+1)

###################
#  compare raw datasets

whole.pca<-vegan::rda()

