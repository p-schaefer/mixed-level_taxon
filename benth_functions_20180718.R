########################################################
#Utility functions
#
########################################################
paste.rep<-function(x,y){ # Utility function needed for describing decision rules
  if (length(x)==1){
    paste0(x,"; ",y)
  } else {
    paste(x,rep(y,times=length(x)),sep="; ")
    
  }
}

print.benth.taxnames<-function(x,...){
  print(as.data.frame(attr(x,"Name_match")[,c(1,2,6)],row.names=1:nrow(attr(x,"Name_match"))))
}

print.benth.taxroll<-function(x,...){
  print(x$decisions[,grepl("Decision",colnames(x$decisions))])
}

print.twogroup_comparision<-function(x,...){
  print(x$table)
}


########################################################
#Function to clean taxonomic names, match with ITIS (secondary matched to NCBI and EOL). 
#
#
########################################################

benth.taxnames<- function(x,
                          #use.invalid=F, # use taxa not currenly valid in ITIS? Not currently implimented
                          gnrDB=c(3,4), #order of prefed DB for looking up names - from gnr_datasources()
                          oligoHairs=T # Treat oligochaete taxa with and without hairs as different taxa
) 
{
  if (!require(taxize,quietly = T)) install.packages('taxize')
  require(taxize,quietly = T)
  
  #x<-c("Simuliidae","Prosimulium mixum","Helodon","Simulium sp.", "Prosimulium/Helodon", "Trichoptera","Trichoptera/Plecoptera")
  
  x<-as.character(x)
  
  if (any(!gnrDB %in% gnr_datasources()[,1])) stop ("gnrDB must be one or more numeric values 
                                               corresponding to database IDs in taxize::gnr_datasources()")
  
  if(!is.vector(x)) stop("x must be a vector of taxon names")
  
  if (any(grepl("indeterminate|immature",x))) stop("Taxon names cannot contain: indeterminate or immature")
  
  # if (any(grepl("Phylum|Subphylum|Class|Order|Family|Subfamily|Tribe| sp| spp
  #               |group| complex| grp| grp|larvae|Larvae|(larvae)|(Larvae)
  #               |adult|Adult|(adult)|(Adult)
  #               |immature|Immature|(immature)|(Immature)|
  #               indet|Indet|pupae|Pupae",x))) stop("Invalid names found (i.e. group, complex, immature, pupae, etc.")
  
  tax.names<-x
  tax.names.orig<-x
  
  #  Section 1: remove useless text from taxa names
  
  
  tax.names <-
    gsub("Phylum|Subphylum|Class|Order|Family|Subfamily|Tribe",
         "",
         tax.names)
  
  tax.names <-
    gsub(" spp| group| complex| grp|larvae|Larvae|(larvae)|(Larvae)
         |adult|Adult|(adult)|(Adult)
         |immature|Immature|(immature)|(Immature)|
         indet|Indet|pupae|Pupae",
         "",
         tax.names,
         fixed = F)
  
  taxa.groups<-tax.names.orig[grepl(" group| complex| grp| spp",tax.names.orig)]
  taxa.slash<-tax.names.orig[grepl("/",tax.names.orig)]
  
  tax.names <- gsub(" sp.", "", tax.names, fixed = T)
  tax.names<-gsub("\\s*\\([^\\)]+\\)","",tax.names)
  
  tax.names<-trimws(tax.names)
  
  #  Section 2: Check to see if any duplicates remain
  if (any(duplicated(tax.names))) stop(paste0("Duplicated taxa not permitted. Please consolidate: ",
                                              paste0(tax.names[duplicated(tax.names)], collapse=", ")))
  
  #  Section 3: Validate each taxon name
  tax.names1<-as.data.frame(resolve("Simuliidae",with_canonical_ranks=T,with_context = T,
                                    best_match_only=T,preferred_data_sources=gnrDB[1],fields="all")$gnr[F,])
  tax.names1<-tax.names1[1:length(tax.names),]
  rownames(tax.names1)<-tax.names
  message("Retrieving taxonomic information...")
  pb <- txtProgressBar(min=1,max=length(tax.names),style = 3)
  
  for (i in 1:length(tax.names)){
    #for debugging
    
    #if (tax.names[i]=="Rhyacophila brunnea/vemna"){
    #  browser()
    #}
    
    setTxtProgressBar(pb, i)
    #browser()
    # a. Check if the name UNIQUELY exists in ITIS
    
    for (n in gnrDB){
      if (!is.na(tax.names1[i,1])) next
      temp1<-resolve(tax.names[i],with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=n,fields="all")
      if(nrow(temp1$gnr)==1) { #If 1 match
        if(!grepl("Animalia|Bilateria",temp1$gnr$classification_path)){  #make sure its an animal (some plants or fungi have identical names as some animals)
          temp1<-resolve(tax.names[i],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=3,fields="all")
          if (any(grepl("Animalia|Bilateria",temp1$gnr$classification_path))){ #This will select the one that is an animal
            tax.names1[i,]<-temp1$gnr[which(grepl("Animalia|Bilateria",temp1$gnr$classification_path))[1],]
          } 
        } else {
          tax.names1[i,]<-temp1$gnr[1,]
        }
      }
      rm(temp1)
    }
  }
  
  #colnames(tax.names1)<-colnames(temp1$gnr)
  #browser()
  #tax.names1$matched_name2[tax.names1$matched_name2!=tax.names1$submitted_name]<-NA #IS THIS CORRECT?
  
  # b. Check any taxa that could not be matched UNIQUELY in ITIS
  if (any(is.na(tax.names1$matched_name2))){
    no.match<-which(is.na(tax.names1$matched_name2))
    error<-0
    resolved<-F
    for (i in no.match) { #For each that couldnt be found
      
      #if (tax.names[i]=="Bezzia/ Palpomyia"){
      #  browser()
      #}
      
      #browser()
      if (grepl("/",tax.names[i],fixed=T)){ # Check if taxa contain a "/"
        #If it does, separate the taxa and evaluate them individually
        #browser()
        error<-error+1
        if (grepl(" .*/.*",tax.names[i])){
          two.tax<-tax.names[i]
          genus<-strsplit(two.tax,"\\s.*")[[1]]
          species<-strsplit(two.tax,"*.\\s")[[1]][2]
          species<-strsplit(species,"/")[[1]]
          species<-trimws(species)
          two.tax<-paste(genus,species,sep=" ")
        } else {
          two.tax<-strsplit(tax.names[i],"/")[[1]]
          two.tax<-trimws(two.tax)
          
        }
        if (any(grepl(paste0(two.tax,collapse="|"),tax.names))){ #If any other taxa are in the dataset, treat undetermiend fraction as own taxa
          
          for (n in gnrDB){
            if (resolved) next
            temp2<-resolve(two.tax,with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=n,fields="all")
            if(nrow(temp2$gnr)==2) {
              if(!any(grepl("Animalia|Bilateria",temp2$gnr$classification_path))){
                temp2.1<-resolve(two.tax[1],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=n,fields="all")
                temp2.2<-resolve(two.tax[2],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=n,fields="all")
                
                if (any(grepl("Animalia|Bilateria",c(temp2.1$gnr$classification_path,temp2.2$gnr$classification_path)))){
                  temp2$gnr<-rbind(temp2.1$gnr[which(grepl("Animalia|Bilateria",temp2$gnr$classification_path))[1],],
                                   temp2.2$gnr[which(grepl("Animalia|Bilateria",temp2$gnr$classification_path))[1],])
                }
              }
              higher.two.tax1<-t(plyr::rbind.fill(lapply(strsplit(temp2$gnr$classification_path,"|",fixed=T),function(x) as.data.frame(t(as.data.frame(x))))))
              higher.two.rank1<-t(plyr::rbind.fill(lapply(strsplit(temp2$gnr$classification_path_ranks,"|",fixed=T),function(x) as.data.frame(t(as.data.frame(x))))))
              higher.shared.rank<-max(na.omit(match(higher.two.rank1[,1], higher.two.rank1[,2])))-1
              
              tax.names1[i,]<-resolve(as.character(higher.two.tax1[higher.shared.rank,1]),with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=n,fields="all")$gnr
              
              tax.names1[i,c(1,2,6,7,8,9,10)]<-c(tax.names[i],NA,tax.names[i],
                                                 paste0(tax.names1$classification_path[i],"|",tax.names[i]),
                                                 paste0(tax.names1$classification_path_ranks[i],"|artificial_",error),
                                                 NA,
                                                 NA)
              other.i<-grep(paste0(two.tax,collapse="|"),tax.names)
              other.i<-other.i[which(!tax.names[i]==rownames(tax.names1[other.i,]))]
              for(n in other.i){
                template<-unlist(strsplit(tax.names1$classification_path[i],"|",fixed=T))
                class.path<-unlist(strsplit(tax.names1$classification_path[n],"|",fixed=T))
                new<-which(is.na(match(class.path,template)))
                new.class.path<-paste0(c(class.path[1:(new[1]-1)],tax.names[i],class.path[new]),collapse="|")
                
                template1<-unlist(strsplit(tax.names1$classification_path_ranks[i],"|",fixed=T))
                class.path1<-unlist(strsplit(tax.names1$classification_path_ranks[n],"|",fixed=T))
                new1<-which(is.na(match(class.path,template)))
                new.class.path1<-paste0(c(class.path1[1:(new1[1]-1)],template1[length(template1)],class.path1[new1]),collapse="|")
                
                tax.names1[n,c(7,8,9)]<-c(new.class.path,
                                          new.class.path1,
                                          NA)
              }
              resolved<-T
              suppressWarnings(rm(temp2,temp2.1,temp2.2,higher.two.tax1,higher.shared.rank,higher.two.rank1,other.i,template,class.path,new,
                                  new.class.path,template1,class.path1,new1,new.class.path1))
            } else {
              #Only 1 taxa found in supplied database, or more than 2 matches found
            }
          }
        } else if (!all(tax.names==two.tax[1],tax.names==two.tax[2])) {#If there are no other taxa in the pair in the dataset, just pick one and use it
          for (n in gnrDB){
            if (!is.na(tax.names1[i,1])) next
            temp1<-resolve(tax.names[i],with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=n,fields="all")
            if(nrow(temp1$gnr)>0) { #If 1 match
              if(!grepl("Animalia|Bilateria",temp1$gnr$classification_path)){  #make sure its an animal (some plants or fungi have identical names as some animals)
                temp1<-resolve(tax.names[i],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=3,fields="all")
                if (any(grepl("Animalia|Bilateria",temp1$gnr$classification_path))){ #This will select the one that is an animal
                  tax.names1[i,]<-temp1$gnr[which(grepl("Animalia|Bilateria",temp1$gnr$classification_path))[1],]
                } 
              } else {
                tax.names1[i,]<-temp1$gnr[1,]
              }
            }
            rm(temp1)
          }
        } 
      }
    }
  }
  
  if (oligoHairs){
    #browser()
    hair<-grep("hair|chaetae|Hair|Chaetae",tax.names1$user_supplied_name)
    if (length(hair)>0){
      hair.dat<-tax.names1[hair,]
      for (i in 1:nrow(hair.dat)) {
        tax.names1[hair[i],c(6,7,8,9,10)]<-c(hair.dat$submitted_name[i],
                                             paste0(hair.dat$classification_path[i],"|",hair.dat$submitted_name[i]),
                                             paste0(hair.dat$classification_path_ranks[i],"|artificial_hair"),
                                             NA,NA)
      }
    }
  }
  
  #browser()
  l.out<-list()
  for (i in 1:nrow(tax.names1)){
    tax<-as.character(strsplit(tax.names1$classification_path[i],"|",fixed=T)[[1]])
    rank<-tolower(as.character(strsplit(tax.names1$classification_path_ranks[i],"|",fixed=T)[[1]]))
    id<-as.character(strsplit(tax.names1$classification_path_ids[i],"|",fixed=T)[[1]])
    if(is.na(tax[1])) next
    #browser()
    if (tax[length(tax)]!=tax.names1$matched_name2[i]){
      #browser()
      tax.names1$matched_name2[i]<-tax[length(tax)]
    }
    
    template<-data.frame(name=tax,
                         rank=rank,
                         id=id,
                         stringsAsFactors=F
    )
    template[template=="Not assigned"]<-NA
    l.out[[i]]<-template
    names(l.out)[i]<-tax.names1$matched_name2[i]
  }
  
  if (any(is.na(tax.names1$user_supplied_name))){
    tax.names1$user_supplied_name[is.na(tax.names1$user_supplied_name)]<-row.names(tax.names1)[is.na(tax.names1$user_supplied_name)]
  }
  
  attr(l.out, 'Name_match')<-tax.names1
  class(l.out)<-"benth.taxnames"
  
  mis.match<-tax.names1[,colnames(tax.names1)%in%c("user_supplied_name","matched_name2")]
  mis.match[is.na(mis.match)]<-"NA"
  mis.match<-mis.match[mis.match[,1]!=mis.match[,2],]
  no.match<-mis.match[mis.match$matched_name2=="NA",]
  if (nrow(no.match)>0)  rownames(no.match)<-1:nrow(no.match)
  
  mis.match<-mis.match[mis.match$matched_name2!="NA",]
  if (nrow(mis.match)>0)   rownames(mis.match)<-1:nrow(mis.match)
  
  if (length(gnrDB)>1){
    multi.source<-tax.names1[,colnames(tax.names1)%in%c("user_supplied_name","matched_name2","data_source_title")]
    multi.source<-multi.source[tax.names1$data_source_title!=gnr_datasources()$title[gnrDB[1]]&!is.na(tax.names1$data_source_title),]
  } else {
    multi.source<-matrix(ncol=0,nrow=0)
  }
  
  if (nrow(mis.match)>0){
    message("")
    message("")
    message("The following taxonomic name changes were made:")
    print(mis.match)
    message("Please review original data matrix, especially if species and genera names are now different, ")
    message("i.e. - Hydropsyche morosa became Ceratopsyche morosa, but Hydropsyche stayed Hydropsyche")
    message("")
  }
  if (nrow(no.match)>0){
    message("")
    message("The following taxa could not be matched:")
    print(data.frame(no.match))
    message("Please manually supply the information with mod.taxnames(), or")
    message("try using additional datasources from gnr_datasources()")
    message("")
  }
  if (length(gnrDB)>1 & nrow(multi.source)>0){
    message("")
    message("The following taxa could not be matched in ITIS:")
    print(data.frame(multi.source))
    message("Please verify higher level taxonomy matches other taxa")
    message("")
  }
  if (length(taxa.slash)>0){
    message("")
    message("Input data contains '/' taxa:")
    print(data.frame(taxa.slash))
    message("Please review these taxa for potential redundancy")
    message("")
  }
  if (length(taxa.groups)>0){
    message("")
    message("The following taxa groups were in the original dataset:")
    print(data.frame(taxa.groups))
    message("Please confirm there are no redundant taxa in these groups")
    message("")
  }
  
  # if (write.output){
  #   saveRDS(l.out,file=paste0(Sys.Date(),"_tax_names(IMPORTANT).rds"))
  # }
  
  return(l.out)
}

########################################################
#Function to modify benth.taxnames(). 
#
#
########################################################

mod.taxnames<-function(taxon,replacement,benth.taxnames) {
  if (class(benth.taxnames)!="benth.taxnames") stop("benth.taxnames must be an output of benth.taxnames()")
  if (length(taxon)!=length(replacement)) stop("taxon and replacement must be the same length")
  if (ncol(replacement[[1]])!=3) stop("replacement must have same table structure as benth.taxnames() objects")
  #browser()
  
  benth.taxnames[match(taxon,names(benth.taxnames))]<-replacement
  names(benth.taxnames[match(taxon,names(benth.taxnames))])<-unlist(lapply(replacement,function(x) x[nrow(x),1]))
  
  Names_match<-attr(benth.taxnames,"Name_match")
  Names_match$classification_path[match(taxon,names(benth.taxnames))]<-unlist(lapply(replacement,function(x)paste0(x[,1],collapse = "|")))
  Names_match$classification_path_ranks[match(taxon,names(benth.taxnames))]<-unlist(lapply(replacement,function(x)paste0(x[,2],collapse = "|")))
  Names_match$classification_path_ids[match(taxon,names(benth.taxnames))]<-unlist(lapply(replacement,function(x)paste0(x[,3],collapse = "|")))
  Names_match$taxon_id[match(taxon,names(benth.taxnames))]<-unlist(lapply(replacement,function(x) x[nrow(x),3]))
  Names_match$matched_name2[match(taxon,names(benth.taxnames))]<-unlist(lapply(replacement,function(x) x[nrow(x),1]))
  Names_match$data_source_title[match(taxon,names(benth.taxnames))]<-"User Entered"
  Names_match$data_source_id[match(taxon,names(benth.taxnames))]<-"NA"
  rownames(Names_match)[match(taxon,names(benth.taxnames))]<-unlist(lapply(replacement,function(x) x[nrow(x),1]))
  
  attr(benth.taxnames,"Name_match")<-Names_match
  
  return(benth.taxnames)
}

########################################################
#Function to get higher level taxonomy in benth.taxnames(). 
#
#
########################################################

higher.taxa<-function(x,taxonomy=NA){
  l.out<-x
  tax.names1<-attr(x,"Name_match")
  tax.final <- l.out
  if (!is.na(taxonomy)) {
    tax.full1 <- taxonomy
  } else {
    tax.full1 <-
      c(
        "kingdom",
        "subkingdom",
        "infrakingdom",
        "superphylum",
        "phylum",
        "subphylum",
        "class",
        "subclass",
        "infraclass",
        "superorder",
        "order",
        "suborder",
        "infraorder",
        "superfamily",
        "family",
        "subfamily",
        "tribe",
        "genus",
        "subgenus",
        "species",
        "subspecies"
      )
  }
  
  if (any(!grepl(paste0(tax.full1,collapse="|"),unique(unlist(lapply(l.out, "[[", 2)))))){
    numb.miss<-unique(unlist(lapply(l.out, "[[", 2)))[which(!grepl(paste0(tax.full1,collapse="|"),unique(unlist(lapply(l.out, "[[", 2)))))]
    for(i in numb.miss){
      if (i=="") next
      has.artif<-which(grepl(i,lapply(l.out, "[[", 2)))
      has.artif<-l.out[[has.artif[which.max(lapply(l.out[has.artif],nrow))]]][,2]
      has.artif<-has.artif[(grep(i,has.artif)-1):(grep(i,has.artif)+1)]
      tax.full1<-append(tax.full1,has.artif[2],which(match(tax.full1,has.artif,nomatch=0)==1))
    }
  }
  
  tax.full1 <-
    tax.full1[tax.full1 %in% as.character(unique(unlist(lapply(tax.final, "[[", 2), use.names =
                                                          F)))]
  tax.full1<-tax.full1[!tax.full1%in%c(
    "kingdom",
    "subkingdom",
    "infrakingdom",
    "superphylum"
  )]
  #browser()
  
  mat.out <- matrix(nrow = length(l.out), ncol = length(tax.full1))
  colnames(mat.out) <- tax.full1
  rownames(mat.out) <- tax.names1$matched_name2
  for (i in l.out) {
    if (any(rownames(mat.out) %in% i$name[nrow(i)])){
      mat.out[rownames(mat.out) %in% i$name[nrow(i)], colnames(mat.out) %in% i$rank] <-
        as.character(i$name[i$rank %in% colnames(mat.out)])
    } else {
      t1<-which(tax.names1$matched_name2==i$name[nrow(i)])
      for (n in 1:length(t1)){
        mat.out[rownames(mat.out) %in% rownames(tax.names1)[t1[n]], colnames(mat.out) %in% i$rank] <-
          as.character(i$name[i$rank %in% colnames(mat.out)])
      }
    }
  }
  return(mat.out)
}



########################################################
#Function to roll up benthic taxonomy. 
#
#
########################################################

benth.taxroll<-function(taxa, #taxa by site matrix
                        taxa.names=NA, #output of benth.taxnames
                        roll.down=T, #apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                        fam.roll.down=T, #apply roll downs in cases of unresolved taxonomy for family level output
                        assign.undet=T, #if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                        fam.assign.undet=T, #same as assign.undet but for family level output
                        Criteria3.percent=0.2, #critical limit for criteria 3
                        Criteria5a.percent=0.5, #critical limit for criteria 5a
                        Criteria5b.numb=2, #critical limit for criteria 5b
                        output.taxonomy=T #add full taxonomy to output matricesS
) {
  if (!require(plyr,quietly = T)) install.packages('plyr')
  require(plyr,quietly = T)
  
  if (!require(taxize,quietly = T)) install.packages('taxize')
  require(taxize,quietly = T)
  
  if (!is.data.frame(taxa)) stop("taxa must be a data frame")
  if (!all(apply(taxa[,-c(1)],2,is.numeric)))  stop("taxa can only contain numberic values past row 1")

  taxa[,1]<-as.character(taxa[,1])
  
  if(class(taxa.names)!="benth.taxnames"){
    message("Starting name search")
    taxa.names<-benth.taxnames(taxa[,1])
    message("Finished name search")
  } else {
    taxa.names<-taxa.names
  }
  
  taxa[is.na(taxa)]<-0
  
  tax<-taxa.names[unique(names(taxa.names))]
  
  tax.full <-
    c("kingdom",
      "subkingdom",
      "infrakingdom",
      "superphylum",
      "phylum",
      "subphylum",
      "class",
      "subclass",
      "infraclass",
      "superorder",
      "order",
      "suborder",
      "infraorder",
      "section",
      "superfamily",
      "family",
      "subfamily",
      "tribe",
      "genus",
      "subgenus",
      "species",
      "subspecies"
    )
  if (any(!grepl(paste0(tax.full,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))){
    numb.miss<-unique(unlist(lapply(tax, "[[", 2)))[which(!grepl(paste0(tax.full,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))]
    for(i in numb.miss){
      if (i=="") next
      has.artif<-which(grepl(i,lapply(tax, "[[", 2)))
      has.artif<-tax[[has.artif[which.max(lapply(tax[has.artif],nrow))]]][,2]
      has.artif<-has.artif[(grep(i,has.artif)-1):(grep(i,has.artif)+1)]
      tax.full<-append(tax.full,has.artif[2],which(match(tax.full,has.artif,nomatch=0)==1))
    }
  }
  
  tax.full<-tax.full[!tax.full%in%c("kingdom",
                                    "subkingdom",
                                    "infrakingdom",
                                    "superphylum")]
  
  tax.full <-
    tax.full[tax.full %in% unique(unlist(lapply(lapply(tax, "[[", 2), tail, n =
                                                  1), use.names = F))]
  
  tax.full.ranks <- lapply(tax, "[[", 2)
  tax.lowest.rank <-
    unlist(lapply(lapply(tax, "[[", 2), tail, n = 1), use.names = T)
  
  tax.full.ranks.ID <- lapply(tax, "[[", 3)
  tax.lowest.rank.ID <-
    unlist(lapply(lapply(tax, "[[", 3), tail, n = 1), use.names = T)
  
  tax.full.ranks.name <- lapply(tax, "[[", 1)
  tax.lowest.rank.name <-
    unlist(lapply(lapply(tax, "[[", 1), tail, n = 1), use.names = T)
  
  #browser()
  
  dat.out<-aggregate(taxa[,-c(1)],by=list(names(taxa.names)),sum)
  dat.out<-dat.out[match(unique(names(taxa.names)),dat.out$Group.1),]
  dat.out<-dat.out[dat.out$Group.1%in%unique(paste0(as.character(tax.lowest.rank.name))),-c(1)]
  row.names(dat.out) <- unique(paste0(as.character(tax.lowest.rank.name)))
  dat.orig<-dat.out
  
  tax.final <- tax
  tax.full1 <-
    c(
      "kingdom",
      "subkingdom",
      "infrakingdom",
      "superphylum",
      "phylum",
      "subphylum",
      "class",
      "subclass",
      "infraclass",
      "superorder",
      "order",
      "suborder",
      "infraorder",
      "superfamily",
      "family",
      "subfamily",
      "tribe",
      "genus",
      "subgenus",
      "species",
      "subspecies"
    )
  
  
  if (any(!grepl(paste0(tax.full1,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))){
    numb.miss<-unique(unlist(lapply(tax, "[[", 2)))[which(!grepl(paste0(tax.full1,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))]
    for(i in numb.miss){
      if (i=="") next
      has.artif<-which(grepl(i,lapply(tax, "[[", 2)))
      has.artif<-tax[[has.artif[which.max(lapply(tax[has.artif],nrow))]]][,2]
      has.artif<-has.artif[(grep(i,has.artif)-1):(grep(i,has.artif)+1)]
      tax.full1<-append(tax.full1,has.artif[2],which(match(tax.full1,has.artif,nomatch=0)==1))
    }
  }
  
  tax.full1 <-
    tax.full1[tax.full1 %in% as.character(unique(unlist(lapply(tax.final, "[[", 2), use.names =
                                                          F)))]
  mat.out <- matrix(nrow = nrow(taxa), ncol = length(tax.full))
  colnames(mat.out) <- tax.full
  rownames(mat.out) <- rownames(attr(taxa.names,"Name_match"))
  for (i in tax.final) {
    if (any(rownames(mat.out) %in% i$name[nrow(i)])){
      mat.out[rownames(mat.out) %in% i$name[nrow(i)], colnames(mat.out) %in% i$rank] <-
        as.character(i$name[i$rank %in% colnames(mat.out)])
    } else {
      t1<-which(attr(taxa.names,"Name_match")$matched_name2==i$name[nrow(i)])
      for (n in 1:length(t1)){
        mat.out[rownames(mat.out) %in% rownames(attr(taxa.names,"Name_match"))[t1[n]], colnames(mat.out) %in% i$rank] <-
          as.character(i$name[i$rank %in% colnames(mat.out)])
      }
    }
  }
  
  decisions.out<-data.frame(Orig.Taxa=rownames(attr(taxa.names,"Name_match")),Matched.Taxa=attr(taxa.names,"Name_match")$matched_name2,Decision="",stringsAsFactors = F)
  decisions.out <- suppressWarnings(cbind(mat.out[match(rownames(attr(taxa.names,"Name_match")),rownames(mat.out)),],decisions.out))
  decisions.out <- cbind(attr(taxa.names,"Name_match")[,c("data_source_title","taxon_id")],decisions.out)
  
  decisions.out$Decision<-as.character(decisions.out$Decision)
  if (any(rowSums(is.na(mat.out))==ncol(mat.out))){
    decisions.out$Decision[which(rowSums(is.na(mat.out))==ncol(mat.out))]<-paste0("Could not match with valid taxa - see matched name column")
  }
  decisions.out$Decision[duplicated(decisions.out$Matched.Taxa)]<-paste0("Duplicated taxa summed")
  
  decisions.out$Decision[is.na(decisions.out$Matched.Taxa)]<-paste0("Could not match with valid taxa - see matched name column")
  decisions.out$Matched.Taxa[is.na(decisions.out$Matched.Taxa)]<-"no match"
  
  message("Starting taxa redundnacy screen:")
  
  for (i in tax.full) {
    #if(i=="genus") stop()
    message(paste0("Screening ",i))
    
    #if(i=="artificial_2") browser()
    
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
        decisions.out$Decision[decisions.out$Matched.Taxa == n] <-
          paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa == n], paste0("0A - Taxon is at lowest level. From ",n))
        
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
        
        decisions.out$Decision[decisions.out$Matched.Taxa == n] <-
          paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa == n],paste0(
            "1A - Only 1 lower level taxa. Recieved abundance from: ", names(other.below.n)))
        
        decisions.out$Decision[decisions.out$Matched.Taxa %in% names(other.below.n)] <-
          paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa %in% names(other.below.n)],paste0(
            "1A - Only 1 lower level taxa. Taxa was rolled up to: ",n))
        
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
        
        if (roll.down) {
          if (assign.undet){
            decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]] <-
              paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]],paste(
                "3A - Abundance less than ",
                Criteria3.percent*100,
                "% of total abundance within taxon and lower. Proportionally assigned abundance to: ",
                paste(at.and.below.n[-c(1)],collapse=", "),
                ". If no lower level taxa identified at the site, abundances assigned to: ", paste0(which.max.below.n)
              ))
            decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[-c(1)]] <-
              paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[-c(1)]],paste(
                "3A - Abundance of higher level taxa less than ",
                Criteria3.percent*100,
                "% of total abundance within ",n," and lower. Recieved abundance from ",n
              ))
            
          } else {
            decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]] <-
              paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]],paste(
                "3A - Abundance less than ",
                Criteria3.percent*100,
                "% of total abundance within taxon and lower. Proportionally assigned abundance to: ", paste(at.and.below.n[-c(1)],collapse=", ")
              ))
            decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[-c(1)]] <-
              paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[-c(1)]],paste(
                "3A - Abundance of higher level taxa less than ",
                Criteria3.percent*100,
                "% of total abundance within ",n," and lower. Recieved abundance from ",n
              ))
          }
          
        } else {
          decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]] <-
            paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]],paste(
              "3A - Abundance less than ",
              Criteria3.percent*100,
              "% of total abundance within taxon and lower. Taxon deleted."
            ))
        }
        
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
          
          decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[-c(1)]] <-
            paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[-c(1)]],paste(
              "5B - Dataset contains only ",Criteria5b.numb," taxa below ", n," and abundance of ",n," is equal to or greater than ",
              Criteria5a.percent*100,
              "% of total abundance within the higher taxon and lower. Abundance was rolled up to: ",paste(at.and.below.n[1],collapse=", ")
            ))
          
          decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[c(1)]] <-
            paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[c(1)]],paste(
              "5B - Taxon contains only ",Criteria5b.numb," lower level taxa and abundance of higher taxa is equal to or greater than ",
              Criteria5a.percent*100,
              "% of total abundance within the higher taxon and lower. Recieved abundance from: ",paste(at.and.below.n[-c(1)],collapse=", ")
            ))
          
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
          
          decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n] <-
            paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n],paste(
              "5C - Taxon contains more than ",Criteria5b.numb," lower level taxa and abundance of higher taxa is equal to or greater than ",
              Criteria5a.percent*100,
              "% of total abundance within the higher taxon and lower. Kept higher and lower level taxa. Applied rule from: ",n
            ))
          
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
        
        if (roll.down) {
          if (assign.undet){
            decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]] <-
              paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]],paste(
                "5A - Taxon contains ",Criteria5b.numb," or more lower level taxa and abundance of higher taxa is less than ",
                Criteria5a.percent*100,
                "% of total abundance within taxon and lower. Proportionally assigned abundance to: ",
                paste(at.and.below.n[-c(1)],collapse=", "),
                ". If no lower level taxa identified at the site, abundances assigned to: ", paste0(which.max.below.n)
              ))
            decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[-c(1)]] <-
              paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[-c(1)]],paste(
                "5A - Recieved abundance from: ",n
              ))
            
          } else {
            decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]] <-
              paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]],paste(
                "5A - Taxon contains ",Criteria5b.numb," or more lower level taxa and abundance of higher taxa is less than ",
                Criteria5a.percent*100,
                "% of total abundance within taxon and lower. Proportionally assigned abundance to: ", paste(at.and.below.n[-c(1)],collapse=", ")
              ))
            decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[-c(1)]] <-
              paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa %in% at.and.below.n[-c(1)]],paste(
                "5A - Recieved abundance from: ",n
              ))
            
          }
          
        } else {
          decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]] <-
            paste.rep(decisions.out$Decision[decisions.out$Matched.Taxa == at.and.below.n[1]],paste(
              "5A - Taxon contains ",Criteria5b.numb," or more lower level taxa and abundance of higher taxa is less than ",
              Criteria5a.percent*100,
              "% of total abundance within taxon and lower. Taxon deleted."
            ))
        }
        
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
  
  d.temp<-do.call(plyr::rbind.fill,lapply(strsplit(decisions.out$Decision,";"),function(x) data.frame(t(data.frame(x)),stringsAsFactors = F)))
  d.temp[is.na(d.temp)]<-""
  colnames(d.temp)<-paste0("Decision ",1:ncol(d.temp))
  decisions.out<-decisions.out[,!colnames(decisions.out)%in%"Decision"]
  decisions.out<-cbind(decisions.out,d.temp)
  
  #browser()
  
  fams<-unique(unlist(lapply(taxa.names, function(x) x[x[2]=="family", 1])))
  
  fams<-fams[!is.na(fams)]
  
  tax.full.ranks.name1<-tax.full.ranks.name[!is.na(names(tax.full.ranks.name))]
  tax.full.ranks1<-tax.full.ranks[!is.na(names(tax.full.ranks.name))]
  
  dat.fam <- data.frame(matrix(nrow=length(fams),ncol=ncol(dat.orig)))
  rownames(dat.fam)<-fams
  colnames(dat.fam)<-colnames(dat.orig)
  for (i in fams) {
    dat.fam[row.names(dat.fam) == i, ] <-
      colSums(dat.orig[grepl(i,tax.full.ranks.name1), ])
  }
  
  dat.orig.mod<-dat.orig
  higher.tax<-tax.lowest.rank[!tax.lowest.rank%in%tax.full1[which(tax.full1=="family"):length(tax.full1)]]
  for (i in names(higher.tax)){
    n.higher<-length(which(lapply(lapply(tax.full.ranks.name1,function(x) which(x==i)),sum)>0))
    if (n.higher<3){
      dat.fam<-rbind(dat.orig.mod[rownames(dat.orig.mod)%in%i,],dat.fam)
    } else if (fam.roll.down) {
      #browser()
      at.i<-tax.full.ranks.name1[grep(i,tax.full.ranks.name1)]
      at.i<-at.i[!names(at.i)%in%i]
      for (x in 1:ncol(dat.fam)){
        if (dat.orig.mod[i,x]==0) next
        if (fam.assign.undet & sum(dat.orig.mod[names(at.i),x])==0){
          #browser()
          highest.prop<-names(at.i)[which.max(rowSums(dat.orig.mod[names(at.i),]))]
          dat.orig.mod[highest.prop,x]<-dat.orig.mod[i,x]
        } else {
          lower.prop<-dat.orig.mod[names(at.i),x]/sum(dat.orig.mod[names(at.i),x])
          dat.orig.mod[names(at.i),x]<-rowSums(cbind(dat.orig.mod[i,x]*lower.prop,dat.orig.mod[names(at.i),x]))
        }
      }
      
      fams.new<-unique(unlist(lapply(at.i, function(x) x[grepl("idae",x)])))
      for (u in fams.new) {
        dat.fam[row.names(dat.fam) == u, ] <-
          colSums(dat.orig.mod[grepl(u,tax.full.ranks.name1), ])
      }
    }
  }
  
  
  lpl.temp<-merge(dat.out,decisions.out[!apply(decisions.out[,colnames(decisions.out)%in%colnames(mat.out)],1,function(x)all(is.na(x))),],all.x=T,all.y=F,by.x = 0,by.y="Matched.Taxa")
  lpl.temp<-lpl.temp[match(rownames(dat.out),lpl.temp$Row.names),]
  rownames(lpl.temp)<-lpl.temp$Row.names
  lpl.temp<-lpl.temp[,match(c(colnames(mat.out),colnames(dat.out)),colnames(lpl.temp))]
  lpl.tax<-cbind(lpl.temp[,colnames(mat.out)])#,rownames(lpl.temp))
  lpl.tax<-apply(lpl.tax,2,as.character)
  lpl.tax.out<-apply(lpl.tax,1,paste0,collapse=";")
  lpl.tax.out<-gsub("NA","",lpl.tax.out)
  
  fam.tax<-resolve(rownames(dat.fam),with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=3,fields="all")$gnr
  l.out<-list()
  for (i in 1:nrow(fam.tax)){
    tax<-as.character(strsplit(fam.tax$classification_path[i],"|",fixed=T)[[1]])
    rank<-tolower(as.character(strsplit(fam.tax$classification_path_ranks[i],"|",fixed=T)[[1]]))
    template<-data.frame(name=tax,
                         rank=rank,
                         stringsAsFactors=F
    )
    l.out[[i]]<-template
    names(l.out)[i]<-fam.tax$matched_name2[i]
  }
  
  fam.mat.out <- matrix(nrow = nrow(dat.fam), ncol = length(tax.full1[-c((which(tax.full1=="family")+1):length(tax.full1))]))
  colnames(fam.mat.out) <- tax.full1[-c((which(tax.full1=="family")+1):length(tax.full1))]
  fam.mat.out<-fam.mat.out[,!colnames(fam.mat.out)%in%c("kingdom","subkingdom","infrakingdom","superphylum","infraclass","superorder","infraorder","superfamily")]
  rownames(fam.mat.out) <- rownames(dat.fam)
  for (i in l.out) {
    if (any(rownames(fam.mat.out) %in% i$name[nrow(i)])){
      fam.mat.out[rownames(fam.mat.out) %in% i$name[nrow(i)], colnames(fam.mat.out) %in% i$rank] <-
        as.character(i$name[i$rank %in% colnames(fam.mat.out)])
    } else {
      t1<-which(attr(taxa.names,"Name_match")$matched_name2==i$name[nrow(i)])
      for (n in 1:length(t1)){
        fam.mat.out[rownames(fam.mat.out) %in% rownames(attr(taxa.names,"Name_match"))[t1[n]], colnames(fam.mat.out) %in% i$rank] <-
          as.character(i$name[i$rank %in% colnames(fam.mat.out)])
      }
    }
  }
  
  fam.tax<-apply(fam.mat.out,2,as.character)
  
  fam.tax.out<-apply(fam.tax,1,paste0,collapse=";")
  fam.tax.out<-gsub("NA","",fam.tax.out)
  
  out<-list()
  out$lpl.matrix<-dat.out
  out$fam.matrix<-dat.fam
  
  out$decisions<-decisions.out
  out$taxa.names<-taxa.names
  
  if (output.taxonomy){
    out$lpl.taxonomy<-lpl.tax.out
    out$fam.taxonomy<-fam.tax.out
    
  }
  
  class(out)<-"benth.taxroll"
  return(out)
}

########################################################
#Function to roll down benthic taxonomy - consistent with old Minnow format. 
#
#
########################################################

benth.rolldown<-function(taxa, #taxa by site matrix
                         taxa.names=NA #output of benth.taxnames
) {
  if (!require(plyr,quietly = T)) install.packages('plyr')
  #require(plyr,quietly = T)
  
  if (!require(taxize,quietly = T)) install.packages('taxize')
  require(taxize,quietly = T)
  
  if (!is.data.frame(taxa)) stop("taxa must be a data frame")
  if (!all(apply(taxa[,-c(1)],2,is.numeric)))  stop("taxa can only contain numberic values past row 1")
  
  CABIN.omit.taxa<-c("Cladocera","Rotifera","Copepoda","Ostracoda","Nematoda","Porifera","Platyhelminthes")
  
  taxa[,1]<-as.character(taxa[,1])
  
  if(class(taxa.names)!="benth.taxnames"){
    message("Starting name search")
    taxa.names<-benth.taxnames(taxa[,1])
    message("Finished name search")
  } else {
    taxa.names<-taxa.names
  }
  
  taxa[is.na(taxa)]<-0
  
  tax<-taxa.names[unique(names(taxa.names))]
  
  tax.full <-
    c("kingdom",
      "subkingdom",
      "infrakingdom",
      "superphylum",
      "phylum",
      "subphylum",
      "class",
      "subclass",
      "infraclass",
      "superorder",
      "order",
      "suborder",
      "infraorder",
      "section",
      "superfamily",
      "family",
      "subfamily",
      "tribe",
      "genus",
      "subgenus",
      "species",
      "subspecies"
    )
  if (any(!grepl(paste0(tax.full,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))){
    numb.miss<-unique(unlist(lapply(tax, "[[", 2)))[which(!grepl(paste0(tax.full,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))]
    for(i in numb.miss){
      if (i=="") next
      has.artif<-which(grepl(i,lapply(tax, "[[", 2)))
      has.artif<-tax[[has.artif[which.max(lapply(tax[has.artif],nrow))]]][,2]
      has.artif<-has.artif[(grep(i,has.artif)-1):(grep(i,has.artif)+1)]
      tax.full<-append(tax.full,has.artif[2],which(match(tax.full,has.artif,nomatch=0)==1))
    }
  }
  
  tax.full<-tax.full[!tax.full%in%c("kingdom",
                                    "subkingdom",
                                    "infrakingdom",
                                    "superphylum")]
  
  tax.full <-
    tax.full[tax.full %in% unique(unlist(lapply(lapply(tax, "[[", 2), tail, n =
                                                  1), use.names = F))]
  
  tax.full.ranks <- lapply(tax, "[[", 2)
  tax.lowest.rank <-
    unlist(lapply(lapply(tax, "[[", 2), tail, n = 1), use.names = T)
  
  tax.full.ranks.ID <- lapply(tax, "[[", 3)
  tax.lowest.rank.ID <-
    unlist(lapply(lapply(tax, "[[", 3), tail, n = 1), use.names = T)
  
  tax.full.ranks.name <- lapply(tax, "[[", 1)
  tax.lowest.rank.name <-
    unlist(lapply(lapply(tax, "[[", 1), tail, n = 1), use.names = T)
  
  #browser()
  
  dat.out<-aggregate(taxa[,-c(1)],by=list(names(taxa.names)),sum)
  dat.out<-dat.out[match(unique(names(taxa.names)),dat.out$Group.1),]
  dat.out<-dat.out[dat.out$Group.1%in%unique(paste0(as.character(tax.lowest.rank.name))),-c(1)]
  row.names(dat.out) <- unique(paste0(as.character(tax.lowest.rank.name)))
  dat.orig<-dat.out
  
  tax.final <- tax
  tax.full1 <-
    c(
      "kingdom",
      "subkingdom",
      "infrakingdom",
      "superphylum",
      "phylum",
      "subphylum",
      "class",
      "subclass",
      "infraclass",
      "superorder",
      "order",
      "suborder",
      "infraorder",
      "superfamily",
      "family",
      "subfamily",
      "tribe",
      "genus",
      "subgenus",
      "species",
      "subspecies"
    )
  
  if (any(!grepl(paste0(tax.full1,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))){
    numb.miss<-unique(unlist(lapply(tax, "[[", 2)))[which(!grepl(paste0(tax.full1,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))]
    for(i in numb.miss){
      if (i=="") next
      has.artif<-which(grepl(i,lapply(tax, "[[", 2)))
      has.artif<-tax[[has.artif[which.max(lapply(tax[has.artif],nrow))]]][,2]
      has.artif<-has.artif[(grep(i,has.artif)-1):(grep(i,has.artif)+1)]
      tax.full1<-append(tax.full1,has.artif[2],which(match(tax.full1,has.artif,nomatch=0)==1))
    }
  }
  
  tax.full1 <-
    tax.full1[tax.full1 %in% as.character(unique(unlist(lapply(tax.final, "[[", 2), use.names =
                                                          F)))]
  message("Starting taxa redundnacy screen:")
  
  for (i in tax.full) {
    #browser()
    
    message(paste0("Screening ",i))
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
  
  fams<-unique(unlist(lapply(taxa.names, function(x) x[x[2]=="family", 1])))
  
  fams<-fams[!is.na(fams)]
  
  tax.full.ranks.name1<-tax.full.ranks.name[!is.na(names(tax.full.ranks.name))]
  tax.full.ranks1<-tax.full.ranks[!is.na(names(tax.full.ranks.name))]
  
  dat.fam <- data.frame(matrix(nrow=length(fams),ncol=ncol(dat.orig)))
  rownames(dat.fam)<-fams
  colnames(dat.fam)<-colnames(dat.orig)
  for (i in fams) {
    dat.fam[row.names(dat.fam) == i, ] <-
      colSums(dat.orig[grepl(i,tax.full.ranks.name1), ])
  }
  
  dat.orig.mod<-dat.orig
  higher.tax<-tax.lowest.rank[!tax.lowest.rank%in%tax.full1[which(tax.full1=="family"):length(tax.full1)]]
  for (i in names(higher.tax)){
    n.higher<-length(which(lapply(lapply(tax.full.ranks.name1,function(x) which(x==i)),sum)>0))
    if (sum(dat.out[rownames(dat.out)%in%i,])>0){
      dat.fam<-rbind(dat.out[rownames(dat.out)%in%i,],dat.fam)
    }
  }
  
  out<-list()
  out$lpl.matrix<-dat.out
  out$fam.matrix<-dat.fam
  
  #browser()
  
  lpl.tax.out<-data.frame(matrix(nrow=length(tax),ncol=length(tax.full)))
  colnames(lpl.tax.out)<-tax.full
  rownames(lpl.tax.out)<-names(tax)
  for (x in 1:length(tax)){
    lpl.tax.out[x,colnames(lpl.tax.out)%in%tax[[x]][,'rank']]<-tax[[x]][tax[[x]][,'rank']%in%colnames(lpl.tax.out),'name']
  }
  
  #browser()
  
  fam.tax<-resolve(rownames(dat.fam),with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=c(3,4),fields="all")$gnr
  
  for (i in unique(fam.tax$matched_name2)){
    if (length(which(fam.tax$matched_name2==i))==2){
      fam.tax<-fam.tax[-c(which(fam.tax$matched_name2==i)[2]),]
    }
  }
  
  l.out<-list()
  for (i in 1:nrow(fam.tax)){
    tax<-as.character(strsplit(fam.tax$classification_path[i],"|",fixed=T)[[1]])
    rank<-tolower(as.character(strsplit(fam.tax$classification_path_ranks[i],"|",fixed=T)[[1]]))
    template<-data.frame(name=tax,
                         rank=rank,
                         stringsAsFactors=F
    )
    l.out[[i]]<-template
    names(l.out)[i]<-fam.tax$matched_name2[i]
  }
  tax.full.fam<-tax.full[1:which(tax.full=="family")]
  fam.tax.out<-data.frame(matrix(nrow=length(l.out),ncol=length(tax.full.fam)))
  colnames(fam.tax.out)<-tax.full.fam
  rownames(fam.tax.out)<-names(l.out)
  for (x in 1:length(l.out)){
    fam.tax.out[x,colnames(fam.tax.out)%in%l.out[[x]][,'rank']]<-l.out[[x]][l.out[[x]][,'rank']%in%colnames(fam.tax.out),'name']
  }
  
  out$lpl.taxonomy<-apply(lpl.tax.out,1,paste0,collapse=";")
  out$lpl.taxonomy<-gsub("NA","",out$lpl.taxonomy)
  out$fam.taxonomy<-apply(fam.tax.out,1,paste0,collapse=";")
  out$fam.taxonomy<-gsub("NA","",out$fam.taxonomy)
  
  class(out)<-"benth.rolldown"
  return(out)
}

########################################################
#Function to roll up benthic taxonomy 
#
#
########################################################

benth.rollup<-function(taxa, #taxa by site matrix
                       taxa.names=NA, #output of benth.taxnames
                       roll.thresh=0.2 #threshold for rolling up taxa
) {
  if (!require(plyr,quietly = T)) install.packages('plyr')
  #require(plyr,quietly = T)
  
  if (!require(taxize,quietly = T)) install.packages('taxize')
  require(taxize,quietly = T)
  
  if (!is.data.frame(taxa)) stop("taxa must be a data frame")
  if (!all(apply(taxa[,-c(1)],2,is.numeric)))  stop("taxa can only contain numberic values past row 1")
  
  CABIN.omit.taxa<-c("Cladocera","Rotifera","Copepoda","Ostracoda","Nematoda","Porifera","Platyhelminthes")
  
  taxa[,1]<-as.character(taxa[,1])
  
  if(class(taxa.names)!="benth.taxnames"){
    message("Starting name search")
    taxa.names<-benth.taxnames(taxa[,1])
    message("Finished name search")
  } else {
    taxa.names<-taxa.names
  }
  
  taxa[is.na(taxa)]<-0
  
  tax<-taxa.names[unique(names(taxa.names))]
  
  tax.full <-
    c("kingdom",
      "subkingdom",
      "infrakingdom",
      "superphylum",
      "phylum",
      "subphylum",
      "class",
      "subclass",
      "infraclass",
      "superorder",
      "order",
      "suborder",
      "infraorder",
      "section",
      "superfamily",
      "family",
      "subfamily",
      "tribe",
      "genus",
      "subgenus",
      "species",
      "subspecies"
    )
  if (any(!grepl(paste0(tax.full,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))){
    numb.miss<-unique(unlist(lapply(tax, "[[", 2)))[which(!grepl(paste0(tax.full,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))]
    for(i in numb.miss){
      if (i=="") next
      has.artif<-which(grepl(i,lapply(tax, "[[", 2)))
      has.artif<-tax[[has.artif[which.max(lapply(tax[has.artif],nrow))]]][,2]
      has.artif<-has.artif[(grep(i,has.artif)-1):(grep(i,has.artif)+1)]
      tax.full<-append(tax.full,has.artif[2],which(match(tax.full,has.artif,nomatch=0)==1))
    }
  }
  
  tax.full<-tax.full[!tax.full%in%c("kingdom",
                                    "subkingdom",
                                    "infrakingdom",
                                    "superphylum")]
  
  tax.full <-
    tax.full[tax.full %in% unique(unlist(lapply(lapply(tax, "[[", 2), tail, n =
                                                  1), use.names = F))]
  
  tax.full.ranks <- lapply(tax, "[[", 2)
  tax.lowest.rank <-
    unlist(lapply(lapply(tax, "[[", 2), tail, n = 1), use.names = T)
  
  tax.full.ranks.ID <- lapply(tax, "[[", 3)
  tax.lowest.rank.ID <-
    unlist(lapply(lapply(tax, "[[", 3), tail, n = 1), use.names = T)
  
  tax.full.ranks.name <- lapply(tax, "[[", 1)
  tax.lowest.rank.name <-
    unlist(lapply(lapply(tax, "[[", 1), tail, n = 1), use.names = T)
  
  #browser()
  
  dat.out<-aggregate(taxa[,-c(1)],by=list(names(taxa.names)),sum)
  dat.out<-dat.out[match(unique(names(taxa.names)),dat.out$Group.1),]
  dat.out<-dat.out[dat.out$Group.1%in%unique(paste0(as.character(tax.lowest.rank.name))),-c(1)]
  row.names(dat.out) <- unique(paste0(as.character(tax.lowest.rank.name)))
  dat.orig<-dat.out
  
  tax.final <- tax
  tax.full1 <-
    c(
      "kingdom",
      "subkingdom",
      "infrakingdom",
      "superphylum",
      "phylum",
      "subphylum",
      "class",
      "subclass",
      "infraclass",
      "superorder",
      "order",
      "suborder",
      "infraorder",
      "superfamily",
      "family",
      "subfamily",
      "tribe",
      "genus",
      "subgenus",
      "species",
      "subspecies"
    )
  
  if (any(!grepl(paste0(tax.full1,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))){
    numb.miss<-unique(unlist(lapply(tax, "[[", 2)))[which(!grepl(paste0(tax.full1,collapse="|"),unique(unlist(lapply(tax, "[[", 2)))))]
    for(i in numb.miss){
      if (i=="") next
      has.artif<-which(grepl(i,lapply(tax, "[[", 2)))
      has.artif<-tax[[has.artif[which.max(lapply(tax[has.artif],nrow))]]][,2]
      has.artif<-has.artif[(grep(i,has.artif)-1):(grep(i,has.artif)+1)]
      tax.full1<-append(tax.full1,has.artif[2],which(match(tax.full1,has.artif,nomatch=0)==1))
    }
  }
  
  tax.full1 <-
    tax.full1[tax.full1 %in% as.character(unique(unlist(lapply(tax.final, "[[", 2), use.names =
                                                          F)))]
  message("Starting taxa redundnacy screen:")
  
  for (i in tax.full) {
    #browser()
    
    message(paste0("Screening ",i))
    at.i <-
      tax.lowest.rank[names(tax.lowest.rank) %in% rownames(dat.out)][which(tax.lowest.rank[names(tax.lowest.rank) %in% rownames(dat.out)] == i)]
    for (n in names(at.i)) {
      #message(paste0("Screening ",n))
      
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
  #browser()
  
  fams<-unique(unlist(lapply(taxa.names, function(x) x[x[2]=="family", 1])))
  
  fams<-fams[!is.na(fams)]
  
  tax.full.ranks.name1<-tax.full.ranks.name[!is.na(names(tax.full.ranks.name))]
  tax.full.ranks1<-tax.full.ranks[!is.na(names(tax.full.ranks.name))]
  
  dat.fam <- data.frame(matrix(nrow=length(fams),ncol=ncol(dat.orig)))
  rownames(dat.fam)<-fams
  colnames(dat.fam)<-colnames(dat.orig)
  for (i in fams) {
    dat.fam[row.names(dat.fam) == i, ] <-
      colSums(dat.orig[grepl(i,tax.full.ranks.name1), ])
  }
  
  dat.orig.mod<-dat.orig
  higher.tax<-tax.lowest.rank[!tax.lowest.rank%in%tax.full1[which(tax.full1=="family"):length(tax.full1)]]
  for (i in names(higher.tax)){
    n.higher<-length(which(lapply(lapply(tax.full.ranks.name1,function(x) which(x==i)),sum)>0))
    if (sum(dat.out[rownames(dat.out)%in%i,])>0){
      dat.fam<-rbind(dat.out[rownames(dat.out)%in%i,],dat.fam)
    }
  }
  
  out<-list()
  out$lpl.matrix<-dat.out
  out$fam.matrix<-dat.fam
  
  #browser()
  
  lpl.tax.out<-data.frame(matrix(nrow=length(tax),ncol=length(tax.full)))
  colnames(lpl.tax.out)<-tax.full
  rownames(lpl.tax.out)<-names(tax)
  for (x in 1:length(tax)){
    lpl.tax.out[x,colnames(lpl.tax.out)%in%tax[[x]][,'rank']]<-tax[[x]][tax[[x]][,'rank']%in%colnames(lpl.tax.out),'name']
  }
  
  #browser()
  
  fam.tax<-resolve(rownames(dat.fam),with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=c(3,4),fields="all")$gnr
  
  for (i in unique(fam.tax$matched_name2)){
    if (length(which(fam.tax$matched_name2==i))==2){
      fam.tax<-fam.tax[-c(which(fam.tax$matched_name2==i)[2]),]
    }
  }
  
  l.out<-list()
  for (i in 1:nrow(fam.tax)){
    tax<-as.character(strsplit(fam.tax$classification_path[i],"|",fixed=T)[[1]])
    rank<-tolower(as.character(strsplit(fam.tax$classification_path_ranks[i],"|",fixed=T)[[1]]))
    template<-data.frame(name=tax,
                         rank=rank,
                         stringsAsFactors=F
    )
    l.out[[i]]<-template
    names(l.out)[i]<-fam.tax$matched_name2[i]
  }
  tax.full.fam<-tax.full[1:which(tax.full=="family")]
  fam.tax.out<-data.frame(matrix(nrow=length(l.out),ncol=length(tax.full.fam)))
  colnames(fam.tax.out)<-tax.full.fam
  rownames(fam.tax.out)<-names(l.out)
  for (x in 1:length(l.out)){
    fam.tax.out[x,colnames(fam.tax.out)%in%l.out[[x]][,'rank']]<-l.out[[x]][l.out[[x]][,'rank']%in%colnames(fam.tax.out),'name']
  }
  
  out$lpl.taxonomy<-apply(lpl.tax.out,1,paste0,collapse=";")
  out$lpl.taxonomy<-gsub("NA","",out$lpl.taxonomy)
  out$fam.taxonomy<-apply(fam.tax.out,1,paste0,collapse=";")
  out$fam.taxonomy<-gsub("NA","",out$fam.taxonomy)
  
  class(out)<-"benth.rollup"
  return(out)
}


