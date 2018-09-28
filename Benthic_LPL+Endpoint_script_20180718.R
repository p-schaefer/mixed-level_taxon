################################################################################
#  Title: Benthic Macroinvertebrate LPL species list correction and endpoint calculation
#  Created by: Patrick Schaefer (pschaefer@minnow.ca)
#  Created on: July 18, 2018
#  R version: 3.4.4
#
#
#  Required input files and data for script
#  1. A matrix of taxonomy data with taxa in rows and stations in columns
#
#  Steps in Script
#  1. Load functions, set working directory, load data set
#  2. Clean taxonomic names and retrieve taxonomic informtion from ITIS
#  3. Perform LPL rollups/rolldowns and calculate family level matrix
#  4. Calculate endpoints at LPL and family level
#
#  Output files
#  1. LPL matrix: lpl_matrix.csv
#  2. Family level matrix: fam_matrix.csv
#  3. Table of decisions and taxonomy: QAQC_output.csv
#  4. Table of endpoins:
#
################################################################################

################################################################################
#  1.  Load functions, set working directory, load data set
rm(list=ls())
library(openxlsx)

source("S:/Templates/Stats Team/SOP/Benthic/benth_functions_20180718.r")

project<-"####"

setwd("S:/Templates/Stats Team/SOP/Benthic/")
dat1 <- read.csv("example.csv")
dat.orig<-dat1
#
################################################################################


################################################################################
#  2. Clean taxonomies and retrieve taxonomic informtion from ITIS

#  Remove empty rows and columns, replace NAs with 0s
dat1 <- dat1[rowSums(dat1[, -c(1)],na.rm=T) > 0, ]
dat1 <- dat1[,c(1,which(colSums(dat1[, -c(1)],na.rm=T) > 0)+1)]
dat1[is.na(dat1)]<-0

#Combine any duplicate taxa and reorder according to input dataset
dat2<-aggregate(dat1[, -c(1)],by=list(dat1[, 1]),sum)
dat2<-dat2[match(unique(dat1[,1]),dat2[,1]),]
rownames(dat2)<-1:nrow(dat2)
dat1<-dat2

tax.names<-benth.taxnames(dat1[, 1], # First column of dataframe should be taxa names
                          oligo.hairs=T, # Should oligochaetes with and without hairs be treated as separate taxa?
                          write.output=T, # should results be saved to an output file for archiving?
                          use.NCBI=T
)

#  To load a benth.taxnames file from an old project use:
#readRDS("filenames.rds") 
#  The filename gets saved as: YYYY-MM-DD_tax_names(IMPORTANT).rds
#  where YYYY-MM-DD is the year month the function was originally run

#  Examine the print out from the function above and the output below to
#  determine if all taxa were correctly identified
tax.names
#  If any taxa were not correctly identified, change them in the original input spreadsheet 
#  and rerun the function

#  This function allows easier visualiziation of the higher level taxonomy
higher.taxa(tax.names)

#  If a taxon classification  MUST be modified, the mod.taxnames() function is an easy way 
#  to do so. It requires the name of the taxa to be replaced, and a new table
#  containing the modified taxonomy. The best way to do so is extract the tables
#  that need to be replaced, and make the changes directly there
taxa.replace<-names(tax.names)[c(1,2)]
taxa.replace.table<-tax.names[taxa.replace]
taxa.replace.table[[1]]$name[taxa.replace.table[[1]]$rank=="superphylum"]<-"modified"
taxa.replace.table[[2]]$name[taxa.replace.table[[2]]$rank=="order"]<-"modified"

tax.names1<-mod.taxnames(taxon=taxa.replace,replacement=taxa.replace.table,benth.taxnames=tax.names)
##################################################################################

##################################################################################
#  3. Perform LPL rollups/rolldowns, calculate endpoints and data summaries

tax.roll1<-benth.taxroll(taxa=dat1, # original input data set
                         taxa.names=tax.names, # output from benth.taxnames()
                         roll.down=T, # apply roll downs in cases of unresolved taxonomy / otherwise will delete higher taxa
                         fam.roll.down=T, # apply roll downs in cases of unresolved taxonomy for family level output
                         assign.undet=T, # if rolling down, sites with no lower level taxa will be assigned to the most abundant taxon in dataset
                         fam.assign.undet=T, # same as assign.undet but for family level output
                         Criteria3.percent=0.2, # critical limit for criteria 3
                         Criteria5a.percent=0.5, # critical limit for criteria 5a
                         Criteria5b.numb=2, # critical limit for criteria 5b
                         CABIN.taxa=F # Exclude taxa based on CABIN standards
)

tax.roll2<-benth.minnowroll(taxa=dat1, # original input data set
                            taxa.names=tax.names, # output from benth.taxnames()
                            CABIN.taxa=F # Exclude taxa based on CABIN standards
                            )

endpoints<-benth.endpoint(tax.roll1)

summaries<-benth.summaries(benth.endpoint=endpoints, # output of benth.endpoint()
                           stations=c("CCDS","CCNF","CCMT"), # station names to group
                           taxa.summary=T, # include taxa summaries
                           abund.thresh=0.05 # for taxa summaries, minimum abundance threshold - per site
                           ) 

bray <-
  benth.bray(
    test = c("CCDS","CCNF"),
    ref = "CCMT", 
    data = tax.roll1$lpl.matrix
  )
###################################################################################

###################################################################################
#  4. Analyze primary endpoints (WIP)

primary.endpoints<-endpoints$lpl.endpoints[,c("Richness","Simpsons_Eveness")]
counts<-colSums(dat1[,-c(1)])
primary.endpoints$Counts<-counts[match(rownames(primary.endpoints),names(counts))]
primary.endpoints$Bray<-bray[match(rownames(primary.endpoints),bray$Site),"Bray-Curtis"]

#two.groups<-twogroup_comparision(ref="CCDS",test="CCMT",primary.endpoints)

plot(two.groups)
###################################################################################
#  5. Write data tables

wb1<-openxlsx::createWorkbook()
openxlsx::addWorksheet(wb1,sheet="raw_data")
writeDataTable(wb1,sheet="raw_data",dat.orig,rowNames=F,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="submitted_data")
writeDataTable(wb1,sheet="submitted_data",dat1,rowNames=F,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="taxa_decisions")
writeDataTable(wb1,sheet="taxa_decisions",tax.roll1$decisions,rowNames=T,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="lpl_matrix")
writeDataTable(wb1,sheet="lpl_matrix",tax.roll1$lpl.matrix,rowNames=T,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="fam_matrix")
writeDataTable(wb1,sheet="fam_matrix",tax.roll1$fam.matrix,rowNames=T,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="lpl_endpoints")
writeDataTable(wb1,sheet="lpl_endpoints",endpoints$lpl.endpoints,rowNames=T,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="fam_endpoints")
writeDataTable(wb1,sheet="fam_endpoints",endpoints$fam.endpoints,rowNames=T,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="lpl_summaries")
writeDataTable(wb1,sheet="lpl_summaries",summaries$lpl.summary,rowNames=F,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="fam_summaries")
writeDataTable(wb1,sheet="fam_summaries",summaries$fam.summary,rowNames=F,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="Bray-Curtis")
writeDataTable(wb1,sheet="Bray-Curtis",bray,rowNames=F,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="lpl_dens_summary")
writeDataTable(wb1,sheet="lpl_dens_summary",summaries$lpl.dens.summary,rowNames=T,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="fam_dens_summary")
writeDataTable(wb1,sheet="fam_dens_summary",summaries$fam.den.summary,rowNames=T,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="lpl_attr")
writeDataTable(wb1,sheet="lpl_attr",endpoints$lpl.attributes,rowNames=F,keepNA=T)

openxlsx::addWorksheet(wb1,sheet="fam_attr")
writeDataTable(wb1,sheet="fam_attr",endpoints$fam.attributes,rowNames=F,keepNA=T)

#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")
openxlsx::saveWorkbook(wb1, paste(project,"_Benthic_report.xlsx"), overwrite = TRUE)

write.csv(tax.roll1$lpl.matrix, paste0(project,"_lpl_matrix.csv"))
write.csv(tax.roll1$fam.matrix, paste0(project,"_fam_matrix.csv"))

wb2<-openxlsx::createWorkbook()
openxlsx::addWorksheet(wb2,sheet="lpl_prop")
writeDataTable(wb2,sheet="lpl_prop",summaries$lpl.prop,rowNames=T,keepNA=T)

openxlsx::addWorksheet(wb2,sheet="fam_prop")
writeDataTable(wb2,sheet="fam_prop",summaries$fam.prop,rowNames=T,keepNA=T)

openxlsx::addWorksheet(wb2,sheet="lpl_prop_summary")
writeDataTable(wb2,sheet="lpl_prop_summary",summaries$lpl.prop.summary,rowNames=T,keepNA=T)

openxlsx::addWorksheet(wb2,sheet="fam_prop_summary")
writeDataTable(wb2,sheet="fam_prop_summary",summaries$fam.prop.summary,rowNames=T,keepNA=T)

#Sys.setenv("R_ZIPCMD" = "C:/Rtools/bin/zip.exe")
openxlsx::saveWorkbook(wb2, paste(project,"_Benthic_proportions.xlsx"), overwrite = TRUE)

##################################################################################

