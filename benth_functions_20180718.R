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
  
  x<-c("Simuliidae","Prosimulium mixum","Helodon","Simulium sp.", "Prosimulium/Helodon")
  
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
    gsub(" spp| group| complex| grp| grp|larvae|Larvae|(larvae)|(Larvae)
         |adult|Adult|(adult)|(Adult)
         |immature|Immature|(immature)|(Immature)|
         indet|Indet|pupae|Pupae",
         "",
         tax.names,
         fixed = F)
  
  taxa.groups<-tax.names.orig[grepl(" group| complex| grp| grp.| spp.",tax.names.orig)]
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
    
    # if(nrow(temp1$gnr)==1) { #If 1 match
    #   if(!grepl("Animalia|Bilateria",temp1$gnr$classification_path)){  #make sure its an animal (some plants or fungi have identical names as some animals)
    #     temp1<-resolve(tax.names[i],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=3,fields="all")
    #     if (any(grepl("Animalia|Bilateria",temp1$gnr$classification_path))){ #This will select the one that is an animal
    #       tax.names1[i,]<-temp1$gnr[which(grepl("Animalia|Bilateria",temp1$gnr$classification_path))[1],]
    #     } 
    #   } else {
    #     tax.names1[i,]<-temp1$gnr[1,]
    #   }
    # } else { # If no match is found, and use.NCBI is an option, look there for the taxon
    #   if (use.NCBI) { 
    #     temp1<-resolve(tax.names[i],with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=4,fields="all")
    #     if(nrow(temp1$gnr)==1) {
    #       if(!grepl("Bilateria",temp1$gnr$classification_path)){
    #         temp1<-resolve(tax.names[i],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=4,fields="all")
    #         if (any(grepl("Bilateria",temp1$gnr$classification_path))){
    #           tax.names1[i,]<-temp1$gnr[which(grepl("Bilateria",temp1$gnr$classification_path))[1],]
    #         } 
    #       } else {
    #         tax.names1[i,]<-temp1$gnr[1,]
    #       }
    #     }
    #   }
    # }
  }
  
  #colnames(tax.names1)<-colnames(temp1$gnr)
  #browser()
  tax.names1$matched_name2[tax.names1$matched_name2!=tax.names1$submitted_name]<-NA
  
  # b. Check any taxa that could not be matched UNIQUELY in ITIS
  if (any(is.na(tax.names1$matched_name2))){
    no.match<-which(is.na(tax.names1$matched_name2))
    error<-0
    for (i in no.match) { #For each that couldnt be found
      
      #if (tax.names[i]=="Bezzia/ Palpomyia"){
      #  browser()
      #}
      
      #browser()
      if (grepl("/",tax.names[i],fixed=T)){# Check if taxa contain a "/"
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
          #browser()
          temp2<-resolve(two.tax,with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=3,fields="all")
          if(nrow(temp2$gnr)==2) {
            if(!any(grepl("Animalia|Bilateria",temp2$gnr$classification_path))){
              temp2.1<-resolve(two.tax[1],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=3,fields="all")
              temp2.2<-resolve(two.tax[2],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=3,fields="all")
              
              if (any(grepl("Animalia|Bilateria",c(temp2.1$gnr$classification_path,temp2.2$gnr$classification_path)))){
                temp2$gnr<-rbind(temp2.1$gnr[which(grepl("Animalia|Bilateria",temp2$gnr$classification_path))[1],],
                                 temp2.2$gnr[which(grepl("Animalia|Bilateria",temp2$gnr$classification_path))[1],])
                higher.two.tax1<-as.data.frame(strsplit(temp2$gnr$classification_path,"|",fixed=T))
                higher.two.rank1<-as.data.frame(strsplit(temp2$gnr$classification_path_ranks,"|",fixed=T))
                #higher.shared.rank<-match(higher.two.rank1[,1], higher.two.rank1[,2])[length(higher.two.rank1[,1])-1]
                higher.shared.rank<-max(na.omit(match(higher.two.rank1[,1], higher.two.rank1[,2])))-1
                
                tax.names1[i,]<-resolve(as.character(higher.two.tax1[higher.shared.rank,1]),with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=3,fields="all")$gnr
                
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
              }
            } else {
              higher.two.tax1<-t(plyr::rbind.fill(lapply(strsplit(temp2$gnr$classification_path,"|",fixed=T),function(x) as.data.frame(t(as.data.frame(x))))))
              higher.two.rank1<-t(plyr::rbind.fill(lapply(strsplit(temp2$gnr$classification_path_ranks,"|",fixed=T),function(x) as.data.frame(t(as.data.frame(x))))))
              #higher.two.tax1<-as.data.frame(strsplit(temp2$gnr$classification_path,"|",fixed=T))
              #higher.two.rank1<-as.data.frame(strsplit(temp2$gnr$classification_path_ranks,"|",fixed=T))
              #higher.shared.rank<-match(higher.two.rank1[,1], higher.two.rank1[,2])[length(higher.two.rank1[,1])-1]
              higher.shared.rank<-max(na.omit(match(higher.two.rank1[,1], higher.two.rank1[,2])))-1
              
              tax.names1[i,]<-resolve(as.character(higher.two.tax1[higher.shared.rank,1]),with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=3,fields="all")$gnr
              
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
            }
          } else {
            if (use.NCBI) {
              temp2<-resolve(two.tax,with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=4,fields="all")
              if(nrow(temp2$gnr)==2) {
                if(!any(grepl("Bilateria",temp2$gnr$classification_path))){
                  temp2.1<-resolve(two.tax[1],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=4,fields="all")
                  temp2.2<-resolve(two.tax[2],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=4,fields="all")
                  
                  if (any(grepl("Bilateria",c(temp2.1$gnr$classification_path,temp2.2$gnr$classification_path)))){
                    temp2$gnr<-rbind(temp2.1$gnr[which(grepl("Bilateria",temp2$gnr$classification_path))[1],],
                                     temp2.2$gnr[which(grepl("Bilateria",temp2$gnr$classification_path))[1],])
                    higher.two.tax1<-as.data.frame(strsplit(temp2$gnr$classification_path,"|",fixed=T))
                    higher.two.rank1<-as.data.frame(strsplit(temp2$gnr$classification_path_ranks,"|",fixed=T))
                    higher.shared.rank<-max(na.omit(match(higher.two.rank1[,1], higher.two.rank1[,2])))-1
                    
                    #higher.shared.rank<-match(higher.two.rank1[,1], higher.two.rank1[,2])[length(higher.two.rank1[,1])-1]
                    tax.names1[i,]<-resolve(as.character(higher.two.tax1[higher.shared.rank,1]),with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=4,fields="all")$gnr
                    
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
                  }
                } else {
                  #browser()
                  higher.two.tax1<-t(plyr::rbind.fill(lapply(strsplit(temp2$gnr$classification_path,"|",fixed=T),function(x) as.data.frame(t(as.data.frame(x))))))
                  higher.two.rank1<-t(plyr::rbind.fill(lapply(strsplit(temp2$gnr$classification_path_ranks,"|",fixed=T),function(x) as.data.frame(t(as.data.frame(x))))))
                  
                  #higher.shared.rank<-match(higher.two.rank1[,1], higher.two.rank1[,2])[length(higher.two.rank1[,1])-1]
                  
                  higher.shared.rank<-max(na.omit(match(higher.two.rank1[,1], higher.two.rank1[,2])))-1
                  tax.names1[i,]<-resolve(as.character(higher.two.tax1[higher.shared.rank,1]),with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=4,fields="all")$gnr
                  
                  tax.names1[i,c(1,2,6,7,8,9,10)]<-c(tax.names[i],NA,tax.names[i],
                                                     paste0(tax.names1$classification_path[i],"|",tax.names[i]),
                                                     paste0(tax.names1$classification_path_ranks[i],"|artificial_",error),
                                                     NA,
                                                     NA)
                  other.i<-grep(paste0(two.tax,collapse="|"),tax.names)
                  other.i<-other.i[which(!tax.names[i]==rownames(tax.names1[other.i,]))]
                  
                  temp2.1<-resolve(tax.names[other.i[1]],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=4,fields="all")
                  temp2.2<-resolve(tax.names[other.i[2]],with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=4,fields="all")
                  
                  if (any(grepl("Bilateria",temp2.1$gnr$classification_path))){
                    tax.names1[other.i[1],]<-temp2.1$gnr[which(grepl("Bilateria",temp2.1$gnr$classification_path))[1],]
                  }
                  if (any(grepl("Bilateria",temp2.2$gnr$classification_path))){
                    tax.names1[other.i[2],]<-temp2.2$gnr[which(grepl("Bilateria",temp2.2$gnr$classification_path))[1],]
                  }
                  
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
                    
                    
                    #template<-unlist(strsplit(tax.names1$classification_path[i],"|",fixed=T))
                    #class.path<-unlist(strsplit(tax.names1$classification_path[n],"|",fixed=T))
                    #new<-which(is.na(match(class.path,template)))
                    #new.class.path<-c(class.path[1:(new[1]-1)],tax.names[i],class.path[new])
                    
                    #template1<-unlist(strsplit(tax.names1$classification_path_ranks[i],"|",fixed=T))
                    #class.path1<-unlist(strsplit(tax.names1$classification_path_ranks[n],"|",fixed=T))
                    #new1<-which(is.na(match(class.path,template)))
                    #new.class.path1<-c(class.path1[1:(new1[1]-1)],template1[length(template1)],class.path1[new1])
                    
                    #tax.names1[n,c(7,8,9)]<-c(new.class.path,
                    #                          new.class.path1,
                    #                          NA)
                  }
                }
              }
            }
          }
        } else if (!all(tax.names==two.tax[1],tax.names==two.tax[2])) {#If there are no other taxa in the pair in the dataset, just pick one and use it
          temp1<-resolve(two.tax,with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=3,fields="all")
          if(nrow(temp1$gnr)>=1) {
            if(!any(grepl("Animalia|Bilateria",temp1$gnr$classification_path))){
              temp1<-resolve(two.tax,with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=3,fields="all")
              if (any(grepl("Animalia|Bilateria",temp1$gnr$classification_path))){
                tax.names1[i,]<-temp1$gnr[which(grepl("Animalia|Bilateria",temp1$gnr$classification_path))[1],]
              } 
            } else {
              tax.names1[i,]<-temp1$gnr[1,]
            }
          } else {
            if (use.NCBI) {
              temp1<-resolve(two.tax,with_canonical_ranks=T,with_context = T,best_match_only=T,preferred_data_sources=4,fields="all")
              if(nrow(temp1$gnr)>=1) {
                if(!any(grepl("Bilateria",temp1$gnr$classification_path))){
                  temp1<-resolve(two.tax,with_canonical_ranks=T,with_context = T,best_match_only=F,preferred_data_sources=4,fields="all")
                  if (any(grepl("Bilateria",temp1$gnr$classification_path))){
                    tax.names1[i,]<-temp1$gnr[which(grepl("Bilateria",temp1$gnr$classification_path))[1],]
                  } 
                } else {
                  tax.names1[i,]<-temp1$gnr[1,]
                }
              }
            }
          }
        } else {#If only 1 other taxa in the pair is present in the dataset
          
        }
        
      }
    }
  }
  
  if (oligo.hairs){
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
  
  if (use.NCBI){
    multi.source<-tax.names1[,colnames(tax.names1)%in%c("user_supplied_name","matched_name2","data_source_title")]
    multi.source<-multi.source[tax.names1$data_source_title!="ITIS"&!is.na(tax.names1$data_source_title),]
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
    message("try setting use.NCBI=TRUE in benth.taxnames()")
    message("")
  }
  if (use.NCBI & nrow(multi.source)>0){
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
    message("Please review these taxa")
    message("Functional traits (feeding and habitat) may not properly be calculated for these taxa")
    message("")
  }
  if (length(taxa.groups)>0){
    message("")
    message("The following taxa groups were in the original dataset:")
    print(data.frame(taxa.groups))
    message("Please confirm there are no redundant taxa in these groups")
    message("")
  }
  
  if (write.output){
    saveRDS(l.out,file=paste0(Sys.Date(),"_tax_names(IMPORTANT).rds"))
  }
  
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
  Names_match$data_source_title[match(taxon,names(benth.taxnames))]<-"Manual override"
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

higher.taxa<-function(x){
  l.out<-x
  tax.names1<-attr(x,"Name_match")
  tax.final <- l.out
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
                        output.taxonomy=T, #add full taxonomy to output matricesS
                        CABIN.taxa=F
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
  
  if (CABIN.taxa){
    lpl.omits<-grepl(paste0(CABIN.omit.taxa,collapse="|"),lpl.tax.out)
    fam.omits<-grepl(paste0(CABIN.omit.taxa,collapse="|"),fam.tax.out)
    
    lpl.tax.out<-lpl.tax.out[!lpl.omits]
    fam.tax.out<-fam.tax.out[!fam.omits]
    
    dat.out<-dat.out[!lpl.omits,]
    dat.fam<-dat.fam[!fam.omits,]
    
    decision.omits<-sapply(attr(taxa.names,"Name_match")$classification_path,function(x) any(grepl(paste0(CABIN.omit.taxa,collapse="|"),x)))
    decisions.out[decision.omits,grepl("Decision",colnames(decisions.out))]<-""
    decisions.out[decision.omits,grepl("Decision",colnames(decisions.out))[1]]<-"CABIN excluded taxa"
    
  }
  
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
                        taxa.names=NA, #output of benth.taxnames
                        CABIN.taxa=F
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
  
  
  if (CABIN.taxa){
    dat.out<-dat.out[!lpl.omits,]
    dat.fam<-dat.fam[!fam.omits,]
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
                       roll.thresh=0.2, #threshold for rolling up taxa
                       CABIN.taxa=F
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
  
  
  if (CABIN.taxa){
    dat.out<-dat.out[!lpl.omits,]
    dat.fam<-dat.fam[!fam.omits,]
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


########################################################
#Function to calculate benthic endpoints
#
#
########################################################

benth.endpoint<-function(x){
  if (!any(class(x)=="benth.taxroll",class(x)=="benth.rolldown",class(x)=="benth.rollup")) stop("Input dataset must be an output from benth.taxroll(), benth.rolldown() or benth.rollup()")
  
  if (!require(vegan,quietly = T)) install.packages('vegan')
  require(vegan,quietly = T)
  
  if (!require(BenthicAnalysistesting,quietly = T)) {
    install.packages('devtools')
    devtools::install_github('p-schaefer/BenthicAnalysistesting')
  }
  require(BenthicAnalysistesting,quietly = T,warn.conflicts=F)
  
  #browser()
  fam.mat<-x$fam.matrix
  lpl.mat<-x$lpl.matrix
  
  endpoint.out <- matrix(nrow = ncol(lpl.mat), ncol = 42)
  rownames(endpoint.out) <- colnames(lpl.mat)
  colnames(endpoint.out) <-
    c(
      "Richness",
      "Rarefied_Richness",
      "Abundance",
      "Perc_Dominance",
      "Simpsons_Diversity",
      "InvSimpsons_Diversity",
      "Shannons_Diversity",
      "Pielous_Evenness",
      "Simpsons_Dominance",
      "Simpsons_Eveness",
      "EPT_Rich",
      "EPT_Perc",
      "Eph_Rich",
      "Eph_Perc",
      "Tric_Rich",
      "Tric_Perc",
      "Plec_Rich",
      "Plec_Perc",
      "Chiron_Perc",
      "Oligo_Perc",
      "Dipt_Rich",
      "Dipt_Perc",
      "HBI",
      "CEFI",
      "Shredder_Perc",
      "Shredder_Rich",
      "Filterer_Perc",
      "Filterer_Rich",
      "Predator_Perc",
      "Predator_Rich",
      "CollGath_Perc",
      "CollGath_Rich",
      "ScrapGraz_Perc",
      "ScrapGraz_Rich",
      "Swimmer_Perc",
      "Swimmer_Rich",
      "Clinger_Perc",
      "Clinger_Rich",
      "Sprawler_Perc",
      "Sprawler_Rich",
      "Climber_Perc",
      "Climber_Rich"
    )
  
  endpoint.out.fam<-data.frame(endpoint.out)
  endpoint.out.lpl<-data.frame(endpoint.out)
  
  fam.out<-t(fam.mat)
  colnames(fam.out)<-x$fam.taxonomy

  lpl.out<-t(lpl.mat)
  colnames(lpl.out)<-x$lpl.taxonomy
  
  lpl.BA.out<-BenthicAnalysistesting::benth.metUI(x=lpl.out)
  fam.BA.out<-BenthicAnalysistesting::benth.metUI(x=fam.out)
  
  endpoint.out.lpl[,1:10]<-
    data.frame(cbind(vegan::specnumber(lpl.out),
                     vegan::rarefy(ceiling(lpl.out),min(rowSums(lpl.out))),
                     rowSums(lpl.out),
                     apply(lpl.out,1,max)/rowSums(lpl.out),
                     vegan::diversity(lpl.out,index="simpson"),
                     vegan::diversity(lpl.out,index="invsimpson"),
                     vegan::diversity(lpl.out,index="shannon"),
                     vegan::diversity(lpl.out,index="shannon")/log(vegan::specnumber(lpl.out)),
                     1/vegan::diversity(lpl.out,index="simpson"),
                     (vegan::diversity(lpl.out,index="invsimpson"))/vegan::specnumber(lpl.out)
    ))
  
  endpoint.out.fam[,1:10]<-
    data.frame(cbind(vegan::specnumber(fam.out),
                     vegan::rarefy(ceiling(fam.out),min(rowSums(fam.out))),
                     rowSums(fam.out),
                     apply(fam.out,1,max)/rowSums(fam.out),
                     vegan::diversity(fam.out,index="simpson"),
                     vegan::diversity(fam.out,index="invsimpson"),
                     vegan::diversity(fam.out,index="shannon"),
                     vegan::diversity(fam.out,index="shannon")/log(vegan::specnumber(fam.out)),
                     1/vegan::diversity(fam.out,index="simpson"),
                     (vegan::diversity(fam.out,index="invsimpson"))/vegan::specnumber(fam.out)
                     
    ))
  
  endpoint.out.lpl$EPT_Rich <-lpl.BA.out$Summary.Metrics$EPT.Richness
  endpoint.out.lpl$EPT_Perc <-fam.BA.out$Summary.Metrics$Percent.EPT
  endpoint.out.lpl$Eph_Rich <-lpl.BA.out$Summary.Metrics$Ephem.Richness
  endpoint.out.lpl$Eph_Perc <-fam.BA.out$Summary.Metrics$Percent.Ephem
  endpoint.out.lpl$Tric_Rich <-lpl.BA.out$Summary.Metrics$Trich.Richness
  endpoint.out.lpl$Tric_Perc <-fam.BA.out$Summary.Metrics$Percent.Trich
  endpoint.out.lpl$Plec_Rich <-lpl.BA.out$Summary.Metrics$Plec.Richness
  endpoint.out.lpl$Plec_Perc <-fam.BA.out$Summary.Metrics$Percent.Plec
  endpoint.out.lpl$Dipt_Rich <- lpl.BA.out$Summary.Metrics$Dipt.Richness
  endpoint.out.lpl$Dipt_Perc <- fam.BA.out$Summary.Metrics$Percent.Dipt
  endpoint.out.lpl$Chiron_Perc <-lpl.BA.out$Summary.Metrics$Percent.Chironomidae
  endpoint.out.lpl$Oligo_Perc <-lpl.BA.out$Summary.Metrics$Percent.Oligochaeta
  endpoint.out.lpl$HBI <-lpl.BA.out$Summary.Metrics$HBI
  endpoint.out.lpl$CEFI <-lpl.BA.out$Summary.Metrics$CEFI
  endpoint.out.lpl$Shredder_Perc<-lpl.BA.out$Summary.Metrics$Shredder.Percent
  endpoint.out.lpl$Shredder_Rich<-lpl.BA.out$Summary.Metrics$Shredder.Richness
  endpoint.out.lpl$Filterer_Perc<-lpl.BA.out$Summary.Metrics$Filterer.Percent
  endpoint.out.lpl$Filterer_Rich<-lpl.BA.out$Summary.Metrics$Filterer.Richness
  endpoint.out.lpl$Predator_Perc<-lpl.BA.out$Summary.Metrics$Predator.Percent
  endpoint.out.lpl$Predator_Rich<-lpl.BA.out$Summary.Metrics$Predator.Richness
  endpoint.out.lpl$CollGath_Perc<-lpl.BA.out$Summary.Metrics$Gatherer.Percent
  endpoint.out.lpl$CollGath_Rich<-lpl.BA.out$Summary.Metrics$Gatherer.Richness
  endpoint.out.lpl$ScrapGraz_Perc<-lpl.BA.out$Summary.Metrics$ScraperGrazer.Percent
  endpoint.out.lpl$ScrapGraz_Rich<-lpl.BA.out$Summary.Metrics$ScraperGrazer.Richness
  endpoint.out.lpl$Swimmer_Perc<-lpl.BA.out$Summary.Metrics$Swimmer.Percent
  endpoint.out.lpl$Swimmer_Rich<-lpl.BA.out$Summary.Metrics$Swimmer.Richness
  endpoint.out.lpl$Clinger_Perc<-lpl.BA.out$Summary.Metrics$Clinger.Percent
  endpoint.out.lpl$Clinger_Rich<-lpl.BA.out$Summary.Metrics$Clinger.Richness
  endpoint.out.lpl$Sprawler_Perc<-lpl.BA.out$Summary.Metrics$Sprawler.Percent
  endpoint.out.lpl$Sprawler_Rich<-lpl.BA.out$Summary.Metrics$Sprawler.Richness
  endpoint.out.lpl$Climber_Perc<-lpl.BA.out$Summary.Metrics$Clinger.Percent
  endpoint.out.lpl$Climber_Rich<-lpl.BA.out$Summary.Metrics$Clinger.Richness
  
  endpoint.out.fam$EPT_Rich <-fam.BA.out$Summary.Metrics$EPT.Richness
  endpoint.out.fam$EPT_Perc <-fam.BA.out$Summary.Metrics$Percent.EPT
  endpoint.out.fam$Eph_Rich <-fam.BA.out$Summary.Metrics$Ephem.Richness
  endpoint.out.fam$Eph_Perc <-fam.BA.out$Summary.Metrics$Percent.Ephem
  endpoint.out.fam$Tric_Rich <-fam.BA.out$Summary.Metrics$Trich.Richness
  endpoint.out.fam$Tric_Perc <-fam.BA.out$Summary.Metrics$Percent.Trich
  endpoint.out.fam$Plec_Rich <-fam.BA.out$Summary.Metrics$Plec.Richness
  endpoint.out.fam$Plec_Perc <-fam.BA.out$Summary.Metrics$Percent.Plec
  endpoint.out.fam$Dipt_Rich <- fam.BA.out$Summary.Metrics$Dipt.Richness
  endpoint.out.fam$Dipt_Perc <- fam.BA.out$Summary.Metrics$Percent.Dipt
  endpoint.out.fam$Chiron_Perc <-fam.BA.out$Summary.Metrics$Percent.Chironomidae
  endpoint.out.fam$Oligo_Perc <-fam.BA.out$Summary.Metrics$Percent.Oligochaeta
  endpoint.out.fam$HBI <-fam.BA.out$Summary.Metrics$HBI
  endpoint.out.fam$CEFI <-fam.BA.out$Summary.Metrics$CEFI
  endpoint.out.fam$Shredder_Perc<-fam.BA.out$Summary.Metrics$Shredder.Percent
  endpoint.out.fam$Shredder_Rich<-fam.BA.out$Summary.Metrics$Shredder.Richness
  endpoint.out.fam$Filterer_Perc<-fam.BA.out$Summary.Metrics$Filterer.Percent
  endpoint.out.fam$Filterer_Rich<-fam.BA.out$Summary.Metrics$Filterer.Richness
  endpoint.out.fam$Predator_Perc<-fam.BA.out$Summary.Metrics$Predator.Percent
  endpoint.out.fam$Predator_Rich<-fam.BA.out$Summary.Metrics$Predator.Richness
  endpoint.out.fam$CollGath_Perc<-fam.BA.out$Summary.Metrics$Gatherer.Percent
  endpoint.out.fam$CollGath_Rich<-fam.BA.out$Summary.Metrics$Gatherer.Richness
  endpoint.out.fam$ScrapGraz_Perc<-fam.BA.out$Summary.Metrics$ScraperGrazer.Percent
  endpoint.out.fam$ScrapGraz_Rich<-fam.BA.out$Summary.Metrics$ScraperGrazer.Richness
  endpoint.out.fam$Swimmer_Perc<-fam.BA.out$Summary.Metrics$Swimmer.Percent
  endpoint.out.fam$Swimmer_Rich<-fam.BA.out$Summary.Metrics$Swimmer.Percent
  endpoint.out.fam$Clinger_Perc<-fam.BA.out$Summary.Metrics$Clinger.Percent
  endpoint.out.fam$Clinger_Rich<-fam.BA.out$Summary.Metrics$Clinger.Richness
  endpoint.out.fam$Sprawler_Perc<-fam.BA.out$Summary.Metrics$Sprawler.Percent
  endpoint.out.fam$Sprawler_Rich<-fam.BA.out$Summary.Metrics$Sprawler.Richness
  endpoint.out.fam$Climber_Perc<-fam.BA.out$Summary.Metrics$Clinger.Percent
  endpoint.out.fam$Climber_Rich<-fam.BA.out$Summary.Metrics$Clinger.Richness
  
  out<-list()
  #browser()
  out$lpl.endpoints<-endpoint.out.lpl
  out$fam.endpoints<-endpoint.out.fam
  out$lpl.attributes<-lpl.BA.out$Attributes
  out$fam.attributes<-fam.BA.out$Attributes
  out$fam.mat<-fam.mat
  out$lpl.mat<-lpl.mat
  
  return(out)
}

########################################################
#Function to calculate Bray-Curtis distances
#
#
########################################################

benth.bray<-function(ref,test,data){
  #browser()
  #if(any(!test%in%colnames(data)))stop("one or more test sites is not in column names of data")
  #if(any(!ref%in%colnames(data)))stop("one or more reference sites is not in column names of data")
  
  ref<-colnames(data)[grepl(paste0(ref,collapse = "|"),colnames(data))]
  test<-colnames(data)[grepl(paste0(test,collapse = "|"),colnames(data))]
  
  if (!require(vegan,quietly = T)) install.packages('vegan')
  require(vegan,quietly = T)

  site.class<-data.frame(cbind(c(ref,test),c(rep("ref",length(ref)),rep("test",length(test)))),stringsAsFactors = F)
  
  ref.taxa<-apply(data[,colnames(data)%in%ref],1,median)
  
  bc.out<-cbind(ref.taxa,data[,colnames(data)%in%ref],data[,colnames(data)%in%test])
  bc.out<-as.matrix(vegdist(t(bc.out),method="bray"))[-c(1),1]
  site.class<-cbind(site.class,bc.out)
  colnames(site.class)<-c("Site","Class","Bray-Curtis")
  
  return(site.class)
}

########################################################
#Function to calculate data summaries
#
#
########################################################

benth.summaries<-function(benth.endpoint,stations,taxa.summary,abund.thresh=NA) {
  if (!require(reshape2,quietly = T)) install.packages('reshape2')
  require(reshape2,quietly = T)
  
  if (!require(plyr,quietly = T)) install.packages('plyr')
  require(plyr,quietly = T)
  
  for (i in stations){
    benth.endpoint$lpl.endpoints$stations[grep(i,rownames(benth.endpoint$lpl.endpoints))]<-i
    benth.endpoint$fam.endpoints$stations[grep(i,rownames(benth.endpoint$fam.endpoints))]<-i
  }
  
  melted.lpl <- melt(benth.endpoint$lpl.endpoints, id.vars=c("stations"))
  melted.fam <- melt(benth.endpoint$fam.endpoints, id.vars=c("stations"))
  
  lpl.summ<-ddply(melted.lpl, c("variable","stations"), summarise,
                  n=length(value[!is.na(value)]),
                  min=min(value,na.rm=T),
                  mean = mean(value,na.rm=T),
                  median=median(value,na.rm=T),
                  max=max(value,na.rm=T),
                  sd = sd(value,na.rm=T),
                  se = sd(value,na.rm=T)/sqrt(length(value[!is.na(value)])))
  fam.summ<-ddply(melted.fam, c( "variable","stations"), summarise,
                  n=length(value[!is.na(value)]),
                  min=min(value,na.rm=T),
                  mean = mean(value,na.rm=T),
                  median=median(value,na.rm=T),
                  max=max(value,na.rm=T),
                  sd = sd(value,na.rm=T),
                  se = sd(value,na.rm=T)/sqrt(length(value[!is.na(value)])))
  
  lpl.summ[,-c(1:3)]<-round(lpl.summ[,-c(1:3)],digits=3)
  fam.summ[,-c(1:3)]<-round(fam.summ[,-c(1:3)],digits=3)
  
  out<-list()
  out$lpl.summary<-lpl.summ
  out$fam.summary<-fam.summ
  
  if (taxa.summary==T) {
    
    fam.mat<-benth.endpoint$fam.mat
    lpl.mat<-benth.endpoint$lpl.mat
    
    for (i in 1:ncol(fam.mat)){
      s<-sum(fam.mat[,i])
      fam.mat[,i]<-fam.mat[,i]/s
    }
    
    for (i in 1:ncol(lpl.mat)){
      s<-sum(lpl.mat[,i])
      lpl.mat[,i]<-lpl.mat[,i]/s
    }

    fam.mat<-fam.mat[order(rowSums(fam.mat),decreasing = T),]
    lpl.mat<-lpl.mat[order(rowSums(lpl.mat),decreasing = T),]
    
    fam.mat<-round(fam.mat,digits=3)
    lpl.mat<-round(lpl.mat,digits=3)
    
    if (!is.na(abund.thresh)){
      lpl.mat<-lpl.mat[apply(lpl.mat,1,function(x)any(x>abund.thresh)),]
      fam.mat<-fam.mat[apply(fam.mat,1,function(x)any(x>abund.thresh)),]
    }
    
    out$lpl.prop<-lpl.mat
    out$fam.prop<-fam.mat
    
    lpl.mat<-data.frame(t(lpl.mat))
    lpl.mat$stations<-benth.endpoint$lpl.endpoints$stations
    
    fam.mat<-data.frame(t(fam.mat))
    fam.mat$stations<-benth.endpoint$fam.endpoints$stations
    
    melted.lpl <- melt(lpl.mat, id.vars=c("stations"))
    melted.fam <- melt(fam.mat, id.vars=c("stations"))
    
    lpl.summ<-ddply(melted.lpl, c("stations", "variable"), summarise,
                    n=length(value[!is.na(value)&value>0]),
                    min=min(value,na.rm=T),
                    mean = mean(value,na.rm=T),
                    median=median(value,na.rm=T),
                    max=max(value,na.rm=T),
                    sd = sd(value,na.rm=T),
                    se = sd(value,na.rm=T)/sqrt(length(value[!is.na(value)])))
    fam.summ<-ddply(melted.fam, c("stations", "variable"), summarise,
                    n=length(value[!is.na(value)&value>0]),
                    min=min(value,na.rm=T),
                    mean = mean(value,na.rm=T),
                    median=median(value,na.rm=T),
                    max=max(value,na.rm=T),
                    sd = sd(value,na.rm=T),
                    se = sd(value,na.rm=T)/sqrt(length(value[!is.na(value)])))
    
    lpl.summ[,-c(1:3)]<-round(lpl.summ[,-c(1:3)],digits=3)
    fam.summ[,-c(1:3)]<-round(fam.summ[,-c(1:3)],digits=3)
    
    out$lpl.prop.summary<-lpl.summ
    out$fam.prop.summary<-fam.summ
    
    fam.mat.c<-benth.endpoint$fam.mat
    lpl.mat.c<-benth.endpoint$lpl.mat
    
    if (!is.na(abund.thresh)){
      lpl.mat.c<-lpl.mat.c[rownames(lpl.mat.c)%in%colnames(lpl.mat),]
      fam.mat.c<-fam.mat.c[rownames(fam.mat.c)%in%colnames(fam.mat),]
    }

    fam.mat.c<-fam.mat.c[order(rowSums(fam.mat.c),decreasing = T),]
    lpl.mat.c<-lpl.mat.c[order(rowSums(lpl.mat.c),decreasing = T),]
    
    fam.mat.c<-round(fam.mat.c,digits=3)
    lpl.mat.c<-round(lpl.mat.c,digits=3)
    
    out$lpl.dens<-lpl.mat.c
    out$fam.dens<-fam.mat.c
    
    lpl.mat.c<-data.frame(t(lpl.mat.c))
    lpl.mat.c$stations<-benth.endpoint$lpl.endpoints$stations
    
    fam.mat.c<-data.frame(t(fam.mat.c))
    fam.mat.c$stations<-benth.endpoint$fam.endpoints$stations
    
    melted.lpl <- melt(lpl.mat.c, id.vars=c("stations"))
    melted.fam <- melt(fam.mat.c, id.vars=c("stations"))
    
    lpl.summ<-ddply(melted.lpl, c("stations", "variable"), summarise,
                    n=length(value[!is.na(value)&value>0]),
                    min=min(value,na.rm=T),
                    mean = mean(value,na.rm=T),
                    median=median(value,na.rm=T),
                    max=max(value,na.rm=T),
                    sd = sd(value,na.rm=T),
                    se = sd(value,na.rm=T)/sqrt(length(value[!is.na(value)])))
    fam.summ<-ddply(melted.fam, c("stations", "variable"), summarise,
                    n=length(value[!is.na(value)&value>0]),
                    min=min(value,na.rm=T),
                    mean = mean(value,na.rm=T),
                    median=median(value,na.rm=T),
                    max=max(value,na.rm=T),
                    sd = sd(value,na.rm=T),
                    se = sd(value,na.rm=T)/sqrt(length(value[!is.na(value)])))
    
    lpl.summ[,-c(1:3)]<-round(lpl.summ[,-c(1:3)],digits=3)
    fam.summ[,-c(1:3)]<-round(fam.summ[,-c(1:3)],digits=3)
    
    out$lpl.dens.summary<-lpl.summ
    out$fam.dens.summary<-fam.summ
  }
  
  return(out)
}

########################################################
#Function to perform two group comparisions
#
#
########################################################

twogroup_comparision<-function(ref,test,data,alpha=0.1){
  if (!require(car,quietly = T)) install.packages('car')
  require(car,quietly = T)
  
  ref<-rownames(data)[grepl(ref,rownames(data))]
  test<-rownames(data)[grepl(test,rownames(data))]
  
  site.class<-data.frame(cbind(c(ref,test),c(rep("ref",length(ref)),rep("exp",length(test)))),stringsAsFactors = F)
  colnames(site.class)<-c("Station","Area")
  
  site.class<-cbind(site.class,data[match(site.class$Station,rownames(data)),])
  
  metrics<-colnames(data)
  
  out<-data.frame(matrix(nrow=length(metrics),ncol=10))
  colnames(out)<-c("Endpoint","Transformation","Test",
                   "P-value","Reference_Mean", "Exposure_Mean","Observed_ES","Observed_percent","Measurable_ES","Power")
  out$Endpoint<-metrics
  for (i in metrics) {
    dat1<-data.frame(matrix(nrow=7,ncol=3))
    colnames(dat1)<-c("Trans","norm","var")
    dat1$Trans<-c("none","log","squared","square root","forth root",
                  "logit","arcsine")
    
    for (n in dat1$Trans) {
      if(class(try(transform(site.class[,i],n),silent=T))!="try-error"){
        if (all(!is.na(transform(site.class[,i],n))) & all(!is.nan(transform(site.class[,i],n)))){
          dat1$norm[dat1$Trans==n]<-shapiro.test(residuals(lm(transform(site.class[,i],n)~as.factor(site.class$Area))))$p.value
          dat1$var[dat1$Trans==n]<-leveneTest(y=transform(site.class[,i],n),group=as.factor(site.class$Area))$`Pr(>F)`[1]
        }
      }
    }
    
    if (any(dat1$norm>0.05,na.rm = T)){
      out$Transformation[out$Endpoint==i]<-dat1$Trans[which.max(dat1$norm)]
      if (dat1$var[which.max(dat1$norm)]>0.05) {
        out$Test[out$Endpoint==i]<-"T-equal"
        out$'P-value'[out$Endpoint==i]<-t.test(transform(site.class[,i],out$Transformation[out$Endpoint==i])~as.factor(site.class$Area),
                                             var.equal=T)$p.value
      } else {
        out$Test[out$Endpoint==i]<-"T-unequal"
        out$'P-value'[out$Endpoint==i]<-t.test(transform(site.class[,i],out$Transformation[out$Endpoint==i])~as.factor(site.class$Area),
                                             var.equal=F)$p.value
      }
      #browser()
      n<-length(which(site.class$Area=="ref"))
      ta<-qt(1-alpha/2,length(which(site.class$Area=="ref"))-1)
      tb<-qt(1-alpha,length(which(site.class$Area=="ref"))-1)
      sd<-sd(transform(site.class[site.class$Area=="ref",i],out$Transformation[out$Endpoint==i]),na.rm=T)
      mse<-mean(residuals(aov(transform(site.class[,i],out$Transformation[out$Endpoint==i])~as.factor(site.class$Area)))^2)
      
      out$Reference_Mean[out$Endpoint==i]<-mean(site.class[site.class$Area=="ref",i],na.rm=T)
      out$Exposure_Mean[out$Endpoint==i]<-mean(site.class[site.class$Area=="exp",i],na.rm=T)
      out$Observed_ES[out$Endpoint==i]<-((transform(out$Exposure_Mean[out$Endpoint==i],out$Transformation[out$Endpoint==i])-
                                           transform(out$Reference_Mean[out$Endpoint==i],out$Transformation[out$Endpoint==i]))/
        sd)
      out$Observed_percent[out$Endpoint==i]<-((transform(out$Exposure_Mean[out$Endpoint==i],out$Transformation[out$Endpoint==i])-
                                            transform(out$Reference_Mean[out$Endpoint==i],out$Transformation[out$Endpoint==i]))/
                                           transform(out$Reference_Mean[out$Endpoint==i],out$Transformation[out$Endpoint==i]))
      
      #browser()
      
      #out$Power[out$Endpoint==i]<-pt(sqrt(n/(2*(sqrt(mse)/(2*sd))^2))-ta,n-1)
      #out$Measurable_ES[out$Endpoint==i]<-((ta+tb)*(sqrt(mse)*(sqrt(2/n))))/sd
    } else {
      out$Transformation[out$Endpoint==i]<-"rank"
      out$Test[out$Endpoint==i]<-"MW"
      out$'P-value'[out$Endpoint==i]<-wilcox.test(site.class[,i]~as.factor(site.class$Area),exact=T)$p.value
      out$Reference_Mean[out$Endpoint==i]<-median(site.class[site.class$Area=="ref",i],na.rm=T)
      out$Exposure_Mean[out$Endpoint==i]<-median(site.class[site.class$Area=="exp",i],na.rm=T)
      out$Observed_ES[out$Endpoint==i]<-(median(site.class[site.class$Area=="exp",i],na.rm=T)-
                                           median(site.class[site.class$Area=="ref",i],na.rm=T))/
                                           mad(site.class[,i],na.rm=T)
      
      out$Observed_percent[out$Endpoint==i]<-(median(site.class[site.class$Area=="exp",i],na.rm=T)-
                                                 median(site.class[site.class$Area=="ref",i],na.rm=T))/
                                                median(site.class[site.class$Area=="ref",i],na.rm=T)
    }
  }
  #browser()
  
  out[,-c(1:3)]<-round(out[,-c(1:3)],digits=3)
  out$Observed_percent<-paste0(out$Observed_percent*100,"%")
  #out[,-c(1:3)]<-scale::comma(out[,-c(1:3,8)])
  
  out1<-list()
  out1$table<-out
  out1$data<-data
  
  class(out1)<-"twogroup_comparision"
  return(out1)
}

########################################################
#Function to plot two group comparisions
#
#
########################################################

plot.twogroup_comparision<-function(x,...) {
  results<-x$table
  data<-x$data
  browser()
}


