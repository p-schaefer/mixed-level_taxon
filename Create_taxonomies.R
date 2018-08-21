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
taxa$full.taxon<-apply(taxa[,colnames(taxa)%in%c("n.p","n.c","n.o","n.f","n.g","n.s")],
                       1,paste0,collapse="-")

ref.out<-t(as.numeric(taxa$count))
colnames(ref.out)<-taxa$full.taxon

ref.coms<-list()
for (x in 1:30){
  ref.coms[[x]]<-vegan::rrarefy(ref.out,500)
}

ref.df<-data.frame(matrix(unlist(ref.coms), nrow=length(taxa$full.taxon), byrow=F))

rownames(ref.df)<-taxa$full.taxon
ref.df<-data.frame(t(ref.df))

#lapply(ref.coms,t)
####################
#  remove taxa

res.ratio<-data.frame(all.ID=c(1,2,3,4,5,5),
                      non.ID=c(20,10,5,4,3,0),
                      some.ID=c(1,1,1,1,1,0)
)
rownames(res.ratio)<-rev(c("n.p","n.c","n.o","n.f","n.g","n.s"))

res.ratio<-t(apply(res.ratio,1,cumsum))

rem.coms<-list()

for (x in 1:30){
  checked<-c(F,F,F,F,F,F)
  names(checked)<-c("n.p","n.c","n.o","n.f","n.g","n.s")
  
  taxa<-t(ref.coms[[x]])
  taxa<-cbind(rownames(taxa),taxa)
  colnames(taxa)<-c("full.taxon","count")
  taxa<-data.frame(taxa)
  taxa<-cbind(taxa,do.call(rbind, strsplit(as.character(taxa$full.taxon),"-")))
  colnames(taxa)[3:8]<-c("n.p","n.c","n.o","n.f","n.g","n.s")
  taxa[c(1,3:8)]<-apply(taxa[c(1,3:8)],2,as.character)
  taxa$count<-as.numeric(as.character(taxa$count))
  
  for (i in names(rev(taxa.list))) {
    if (i=="n.s"){
      for (n in unique(taxa[,i])) {
        if (taxa$count[taxa$n.s==n]==0) next
        
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
        if (taxa$count[taxa$full.taxon==n]==0) next
        
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
  
  rem.coms[[x]]<-taxa
}

rem.out<-lapply(rem.coms,"[[",'count')
rem.out<-lapply(rem.out,data.frame)
rem.out.names<-lapply(rem.coms,"[[",'full.taxon')
for (x in 1:30){
  row.names(rem.out[[x]])<-rem.out.names[[x]]
}
rem.out<-lapply(rem.out,t)
#rem.out<-lapply(rem.out,data.frame)

unique.taxa<-unique(unlist(lapply(rem.out,colnames)))

rem.df<-data.frame(matrix(nrow=30,ncol=length(unique.taxa)))
colnames(rem.df)<-unique.taxa

unique.taxa.all<-unique(c(unlist(lapply(rem.out,colnames)),unlist(lapply(rem.coms,"[[",1))))

all.df<-data.frame(matrix(nrow=30,ncol=length(unique.taxa.all)))
colnames(all.df)<-unique.taxa.all

for (x in 1:30) {
  rem.df[x,]<-rem.out[[x]][,match(unique.taxa,colnames(rem.out[[x]]),nomatch =NA)]
}

all.df<-data.frame(matrix(nrow=60,ncol=length(unique.taxa.all)))
colnames(all.df)<-unique.taxa.all
rownames(all.df)<-c(paste0("rem",1:30),paste0("ref",1:30))

for (x in 1:30) {
  all.df[x,]<-rem.out[[x]][,match(unique.taxa.all,colnames(rem.out[[x]]),nomatch =NA)]
  all.df[x+30,]<-ref.coms[[x]][,match(unique.taxa.all,colnames(ref.coms[[x]]),nomatch =NA)]
}

all.df[is.na(all.df)]<-0
all.df<-all.df[,colSums(all.df)>0]

rem.df[is.na(rem.df)]<-0
rem.df<-rem.df[,colSums(rem.df)>0]

ref.df[is.na(ref.df)]<-0
ref.df<-ref.df[,colSums(ref.df)>0]

#heatmap(log(as.matrix(all.df)+1))

ord1<-vegan::cca(log10(rem.df+1))
vegan::ordiplot(vegan::cca(log10(all.df+1),scale=T),type="t",scaling="symmetric",display = "sites")


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

whole.cca<-vegan::cca(log10(ref.df+1))
rem.cca<-vegan::cca(log10(rem.df+1))
proc<-vegan::procrustes(whole.cca,rem.cca)
vegan::protest(whole.cca,rem.cca)
