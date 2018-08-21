
########################################################
#Utility functions
#
########################################################

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

paste.rep<-function(x,y){ # Utility function needed for describing decision rules
  if (length(x)==1){
    paste0(x,"; ",y)
  } else {
    paste(x,rep(y,times=length(x)),sep="; ")
    
  }
}

print.benth.taxroll<-function(x,...){
  print(x$decisions[,grepl("Decision",colnames(x$decisions))])
}

########################################################
#Taxa roll decision key function
#
########################################################

benth.taxroll<-function(taxa, #taxa by site matrix
                        taxa.names=NA, #output of benth.taxnames
                        roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                        assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                        Criteria3.percent=0.2, #critical limit for criteria 3
                        Criteria5a.percent=0.5, #critical limit for criteria 5a
                        Criteria5b.numb=2 #critical limit for criteria 5b
) {

  if (!is.data.frame(taxa)) stop("taxa must be a data frame")
  if (!all(apply(taxa[,-c(1)],2,is.numeric)))  stop("taxa can only contain numberic values past row 1")
  
  taxa[,1]<-as.character(taxa[,1])
  
  taxa[is.na(taxa)]<-0
  
  tax<-taxa.names[unique(names(taxa.names))]
  
  tax.full.ranks <- lapply(tax, "[[", 2)
  tax.lowest.rank <-
    unlist(lapply(lapply(tax, "[[", 2), tail, n = 1), use.names = T)
  
  tax.full.ranks.name <- lapply(tax, "[[", 1)
  tax.lowest.rank.name <-
    unlist(lapply(lapply(tax, "[[", 1), tail, n = 1), use.names = T)
  
  dat.out<-aggregate(taxa[,-c(1)],by=list(names(taxa.names)),sum)
  dat.out<-dat.out[match(unique(names(taxa.names)),dat.out$Group.1),]
  dat.out<-dat.out[dat.out$Group.1%in%unique(paste0(as.character(tax.lowest.rank.name))),-c(1)]
  row.names(dat.out) <- unique(paste0(as.character(tax.lowest.rank.name)))
  dat.orig<-dat.out
  
  for (i in c("phylum","class","order","family","genus","species")) {

    at.i <-
      tax.lowest.rank[names(tax.lowest.rank) %in% rownames(dat.out)][which(tax.lowest.rank[names(tax.lowest.rank) %in% rownames(dat.out)] == i)]
    for (n in names(at.i)) {
      over.n <- tax[[n]]
      
      if (any(is.na(over.n))& !any(grepl("artificial",over.n))) {
        stop("NA error in taxonomy.")
      }
      
      #criteria 1 - taxon at lowest level
      other.at.n <-
        which(sapply(tax.full.ranks.name[!names(tax.full.ranks.name) %in% n &
                                           names(tax.full.ranks.name) %in% rownames(dat.out)], function(x)
                                             any(as.character(x) == as.character(over.n$name[nrow(over.n)]))))
      if (length(other.at.n) == 0) {

        suppressWarnings(
          rm(
            other.at.n,
            other.below.n,
            at.and.below.n,
            abund.at.n,
            abund.below.n,
            prop.below,
            which.max.below.n
          )
        )
        next
      }
      #criteria 2 - taxon contains only 1 lower-level identification
      other.below.n <-other.at.n
      if (length(other.below.n) == 1) {
        dat.out[over.n$name[nrow(over.n)], ] <-
          colSums(dat.out[c(over.n$name[nrow(over.n)], names(other.below.n)), ])
        dat.out <- dat.out[!row.names(dat.out) %in% names(other.below.n), ]
        
        suppressWarnings(
          rm(
            over.n,
            other.at.n,
            other.below.n,
            at.and.below.n,
            abund.at.n,
            abund.below.n,
            prop.below,
            which.max.below.n
          )
        )
        
        next
      }
      #criteria 3 - Abundance at higher level is <20% of
      #             the total abundance within that taxon and lower
      at.and.below.n <-
        c(over.n$name[nrow(over.n)], names(other.below.n))
      abund.at.n <-
        sum(dat.out[rownames(dat.out) %in% at.and.below.n[1], ])
      abund.below.n <-
        sum(rowSums(dat.out[rownames(dat.out) %in% at.and.below.n[-c(1)], ]))
      which.max.below.n<-names(which.max(rowSums(dat.out[rownames(dat.out) %in% at.and.below.n[-c(1)], ])))
      if (abund.at.n / sum(abund.at.n, abund.below.n) <= Criteria3.percent) {
        if (roll.down) {
          for (x in 1:ncol(dat.out)) {
            prop.below <-
              dat.out[at.and.below.n[-c(1)], x] / sum(dat.out[at.and.below.n[-c(1)], x])
            if (all(is.nan(prop.below))) {
              if (!assign.undet) {
                next
              } else {
                dat.out[rownames(dat.out) %in% which.max.below.n, x]<-dat.out[rownames(dat.out) %in% at.and.below.n[1], x]
                next
                #assign undetermined taxa to most abundant taxon
              }
            } 
            dat.out[rownames(dat.out) %in% at.and.below.n[-c(1)], x] <-
              rowSums(cbind(prop.below * dat.out[at.and.below.n[1], x], dat.out[rownames(dat.out) %in%
                                                                                  at.and.below.n[-c(1)], x]))
          }
        }
        dat.out <- dat.out[!row.names(dat.out) %in% at.and.below.n[1], ]
        suppressWarnings(
          rm(
            over.n,
            other.at.n,
            other.below.n,
            at.and.below.n,
            abund.at.n,
            abund.below.n,
            prop.below,
            which.max.below.n
          )
        )
        
        next
      }
      
      #criteria 5a - Higher-level taxon contains 2 lower-level identifications and >50%
      #             of the total abundance within that taxon and lower
      if (abund.at.n / sum(abund.at.n, abund.below.n) >= Criteria5a.percent) {
        if (!length(at.and.below.n[-c(1)]) > Criteria5b.numb) {
          dat.out[row.names(dat.out) %in% at.and.below.n[1], ] <-
            colSums(dat.out[rownames(dat.out) %in% at.and.below.n, ])
          dat.out <-
            dat.out[!rownames(dat.out) %in% at.and.below.n[-c(1)], ]
          
          suppressWarnings(
            rm(
              over.n,
              other.at.n,
              other.below.n,
              at.and.below.n,
              abund.at.n,
              abund.below.n,
              prop.below,
              which.max.below.n
            )
          )
          next
        } else {
          next
        } 
      } else {
        if (roll.down) {
          for (x in 1:ncol(dat.out)) {
            prop.below <-
              dat.out[at.and.below.n[-c(1)], x] / sum(dat.out[at.and.below.n[-c(1)], x])
            if (all(is.nan(prop.below))) {
              if (!assign.undet) {
                next
              } else {
                dat.out[rownames(dat.out) %in% which.max.below.n, x]<-dat.out[rownames(dat.out) %in% at.and.below.n[1], x]
                next
                #assign undetermined taxa to most abundant taxon
              }
            } 
            dat.out[rownames(dat.out) %in% at.and.below.n[-c(1)], x] <-
              rowSums(cbind(prop.below * dat.out[at.and.below.n[1], x], dat.out[rownames(dat.out) %in%
                                                                                  at.and.below.n[-c(1)], x]))
          }
        }
        dat.out <- dat.out[!row.names(dat.out) %in% at.and.below.n[1], ]
        
        suppressWarnings(
          rm(
            over.n,
            other.at.n,
            other.below.n,
            at.and.below.n,
            abund.at.n,
            abund.below.n,
            prop.below,
            which.max.below.n
          )
        )
        next
      }
    }
  }
  
  #browser()
  
  out<-dat.out

  return(out)
}

########################################################
#Taxa roll down function
#
########################################################

benth.rolldown<-function(taxa, #taxa by site matrix
                         taxa.names=NA #output of benth.taxnames
) {
  if (!is.data.frame(taxa)) stop("taxa must be a data frame")
  if (!all(apply(taxa[,-c(1)],2,is.numeric)))  stop("taxa can only contain numberic values past row 1")
  
  taxa[,1]<-as.character(taxa[,1])
  
  taxa[is.na(taxa)]<-0
  
  tax<-taxa.names[unique(names(taxa.names))]
  
  tax.full.ranks <- lapply(tax, "[[", 2)
  tax.lowest.rank <-
    unlist(lapply(lapply(tax, "[[", 2), tail, n = 1), use.names = T)
  
  tax.full.ranks.name <- lapply(tax, "[[", 1)
  tax.lowest.rank.name <-
    unlist(lapply(lapply(tax, "[[", 1), tail, n = 1), use.names = T)
  
  dat.out<-aggregate(taxa[,-c(1)],by=list(names(taxa.names)),sum)
  dat.out<-dat.out[match(unique(names(taxa.names)),dat.out$Group.1),]
  dat.out<-dat.out[dat.out$Group.1%in%unique(paste0(as.character(tax.lowest.rank.name))),-c(1)]
  row.names(dat.out) <- unique(paste0(as.character(tax.lowest.rank.name)))
  dat.orig<-dat.out
  
  for (i in c("phylum","class","order","family","genus","species")) {

    at.i <-
      tax.lowest.rank[names(tax.lowest.rank) %in% rownames(dat.out)][which(tax.lowest.rank[names(tax.lowest.rank) %in% rownames(dat.out)] == i)]
    for (n in names(at.i)) {
      over.n <- tax[[n]]
      
      if (any(is.na(over.n))& !any(grepl("artificial",over.n))) {
        stop("NA error in taxonomy.")
      }
      
      other.at.n <-
        which(sapply(tax.full.ranks.name[!names(tax.full.ranks.name) %in% n &
                                           names(tax.full.ranks.name) %in% rownames(dat.out)], function(x)
                                             any(as.character(x) == as.character(over.n$name[nrow(over.n)]))))
      if (length(other.at.n) == 0) {
        suppressWarnings(
          rm(
            other.at.n
          )
        )
        next
      }
      if (length(other.at.n) == 1) {
        dat.out[rownames(dat.out) %in% names(other.at.n), ] <-
          dat.out[rownames(dat.out) %in% names(other.at.n), ] +
          dat.out[rownames(dat.out) %in% over.n$name[nrow(over.n)], ]
        dat.out[rownames(dat.out) %in% over.n$name[nrow(over.n)], ]<-0
        suppressWarnings(
          rm(
            other.at.n
          )
        )
        next
      }
      
      at.and.below.n <-
        c(over.n$name[nrow(over.n)], names(other.at.n))
      abund.at.n <-
        dat.out[rownames(dat.out) %in% at.and.below.n[1], ]
      abund.below.n <-
        dat.out[rownames(dat.out) %in% at.and.below.n[-c(1)], ]
      prop.below.n<-apply(abund.below.n,2,function(x) x/sum(x))
      prop.below.n[is.nan(prop.below.n)]<-0
      for (x in 1:ncol(prop.below.n)){
        dat.out[rownames(dat.out) %in% at.and.below.n[-c(1)], x] <-
          abund.below.n[,x]+(prop.below.n[,x]*as.numeric(abund.at.n[x]))
        if (sum(prop.below.n[,x])!=0){
          dat.out[rownames(dat.out) %in% over.n$name[nrow(over.n)], x]<-0
        }
      }
      suppressWarnings(
        rm(
          other.at.n,
          other.below.n,
          at.and.below.n,
          abund.at.n,
          abund.below.n,
          prop.below
        )
      )
    }
  }
  #browser()
  
  out<-dat.out
  
  return(out)
}

########################################################
#Taxa roll up function
#
########################################################

benth.rollup<-function(taxa, #taxa by site matrix
                       taxa.names=NA, #output of benth.taxnames
                       roll.thresh=0.2 #threshold for rolling up taxa
) {
  if (!is.data.frame(taxa)) stop("taxa must be a data frame")
  if (!all(apply(taxa[,-c(1)],2,is.numeric)))  stop("taxa can only contain numberic values past row 1")
  
  taxa[,1]<-as.character(taxa[,1])
  
  taxa[is.na(taxa)]<-0
  
  tax<-taxa.names[unique(names(taxa.names))]
  
  tax.full.ranks <- lapply(tax, "[[", 2)
  tax.lowest.rank <-
    unlist(lapply(lapply(tax, "[[", 2), tail, n = 1), use.names = T)
  
  tax.full.ranks.name <- lapply(tax, "[[", 1)
  tax.lowest.rank.name <-
    unlist(lapply(lapply(tax, "[[", 1), tail, n = 1), use.names = T)
  
  dat.out<-aggregate(taxa[,-c(1)],by=list(names(taxa.names)),sum)
  dat.out<-dat.out[match(unique(names(taxa.names)),dat.out$Group.1),]
  dat.out<-dat.out[dat.out$Group.1%in%unique(paste0(as.character(tax.lowest.rank.name))),-c(1)]
  row.names(dat.out) <- unique(paste0(as.character(tax.lowest.rank.name)))
  dat.orig<-dat.out
  
  for (i in c("phylum","class","order","family","genus","species")) {
    at.i <-
      tax.lowest.rank[names(tax.lowest.rank) %in% rownames(dat.out)][which(tax.lowest.rank[names(tax.lowest.rank) %in% rownames(dat.out)] == i)]
    for (n in names(at.i)) {
      over.n <- tax[[n]]
      
      if (any(is.na(over.n))& !any(grepl("artificial",over.n))) {
        stop("NA error in taxonomy.")
      }
      
      other.at.n <-
        which(sapply(tax.full.ranks.name[!names(tax.full.ranks.name) %in% n &
                                           names(tax.full.ranks.name) %in% rownames(dat.out)], function(x)
                                             any(as.character(x) == as.character(over.n$name[nrow(over.n)]))))
      if (length(other.at.n) == 0) {
        suppressWarnings(
          rm(
            other.at.n
          )
        )
        next
      }
      
      at.and.below.n <-
        c(over.n$name[nrow(over.n)], names(other.at.n))
      abund.at.n <-
        dat.out[rownames(dat.out) %in% at.and.below.n[1], ]
      abund.below.n <-
        dat.out[rownames(dat.out) %in% at.and.below.n[-c(1)], ]
      thresh<-sum(abund.at.n)/(sum(rowSums(abund.below.n))+sum(abund.at.n))
      thresh[is.nan(thresh)]<-0
      
      if (thresh<(roll.thresh)){
        dat.out[rownames(dat.out) %in% over.n$name[nrow(over.n)], ]<-0
      } else {
        dat.out[rownames(dat.out) %in% over.n$name[nrow(over.n)], ]<-
          dat.out[rownames(dat.out) %in% over.n$name[nrow(over.n)], ]+
          colSums(dat.out[rownames(dat.out) %in% at.and.below.n[-c(1)], ])
        dat.out[rownames(dat.out) %in% at.and.below.n[-c(1)], ]<-0
      }
      suppressWarnings(
        rm(
          other.at.n,
          other.below.n,
          at.and.below.n,
          abund.at.n,
          abund.below.n,
          prop.below
        )
      )
    }
  }
  
  out<-dat.out
  return(out)
}
