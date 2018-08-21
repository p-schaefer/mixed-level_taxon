########################################################
#Create a single master taxonomy reflecting all taxa abundances at a site
#
########################################################

master.taxonomy<-function(p.n=5,
                          p.s=2,
                          p.m=2,
                          c.s=2,
                          c.m=2,
                          o.s=2,
                          o.m=2,
                          f.s=2,
                          f.m=2,
                          g.s=2,
                          g.m=2,
                          s.mean=5,
                          s.sd=2){
  source("functions.R")
  
  n.p<-rnbinom(p.n,size=p.s,mu=p.m)
  n.p[n.p==0]<-1
  n.c<-rnbinom(sum(n.p),size=c.s,mu=c.m)
  n.c[n.c==0]<-1
  n.o<-rnbinom(sum(n.c),size=o.s,mu=o.m)
  n.o[n.o==0]<-1
  n.f<-rnbinom(sum(n.o),size=f.s,mu=f.m)
  n.f[n.f==0]<-1
  n.g<-rnbinom(sum(n.f),size=g.s,mu=g.m)
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
  
  taxa$count<-ceiling(rlnorm(nrow(taxa),s.mean,s.sd))
  taxa$full.taxon<-apply(taxa[,colnames(taxa)%in%c("n.p","n.c","n.o","n.f","n.g","n.s")],
                         1,paste0,collapse="-")
  
  ref.out<-t(as.numeric(taxa$count))
  colnames(ref.out)<-taxa$full.taxon
  
  return(ref.out)
}


########################################################
#randomly move individuals up the taxonomic heirarchy
#
########################################################

taxonomy.shuffle<-function(ref.coms,res.ratio){
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
  
  rem.taxlist<-list()
  
  for (i in 1:length(unique.taxa)) {
    x<-strsplit(unique.taxa[i],"-")[[1]]
    t1<-data.frame(name=as.character(x),rank=c("phylum","class","order","family","genus","species"))
    t1<-t1[t1$name!="NA",]
    t1$name<-as.character(t1$name)
    t1$rank<-as.character(t1$rank)
    rem.taxlist[[i]]<-t1
  }
  
  names(rem.taxlist)<-unlist(lapply(lapply(rem.taxlist, "[[", 1), tail, n = 1), use.names = T)
  
  for.taxroll<-data.frame(t(rem.df))
  for.taxroll<-cbind(unlist(lapply(lapply(rem.taxlist, "[[", 1), tail, n = 1), use.names = T),for.taxroll)
  colnames(for.taxroll)<-c("taxa",paste0("rem",1:30))
  rownames(for.taxroll) <-
    unlist(lapply(lapply(rem.taxlist, "[[", 1), tail, n = 1), use.names = T)
  
  out<-list()
  
  out$rem.coms<-rem.coms
  out$ref.coms<-ref.coms
  
  out$rem.df<-rem.df
  out$ref.df<-ref.df
  out$all.df<-all.df

  out$rem.df[is.na(out$rem.df)]<-0
  out$ref.df[is.na(out$ref.df)]<-0
  out$all.df[is.na(out$all.df)]<-0
  
  out$rem.taxlist<-rem.taxlist
  out$for.taxroll<-for.taxroll
  
  return(out)
}



