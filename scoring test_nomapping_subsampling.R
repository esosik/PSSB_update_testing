###TEST Script #3: No mapping, Subsampling

###In this test, James Develle has temporarily changed the subsampling routine in the dev version of PSSB to be non-random.
###The PSSB subsampling routine instead orders taxa first by uniqueness (TRUE to FALSE) and then by *original* TSN (Low to High).
##Though the subsampling occurs after rollup, the subsampling routine uses the TSN of the original taxa, prior to any STE rollup or mapping.
##The result is a list of taxa in order of TSN from low to high. Where two identical taxa are in the sample (as in pupa and larva), 
##the unique taxa appears in the list first, followed immediately by the non-unique taxa, and then continuing the list with the next highest TSN.
##This R code tests that the subsampling routine is functioning as specified, and can replicate scores. Any score mismatches (not already identified in previous tests) would indicate an unexpected behavior in the PSSB subsampling function.
##For this test, I re-released samples from 09SOO1130 (ONLY IN THE DEV VERSION OF PSSB) to trigger the new subsampling routine. I then downloaded the scores from PSSB.
## I replicate in my R code the expected data processing, and compare the resulting scores to the PSSB output.

##Based on this test, I did not identify any new issues.





library(plyr)
library(openxlsx)
library(lubridate)
library(stringr)
library(tidyverse)
library(dplyr)
getwd()

setwd(here::here())




##This function takes the raw taxa data .txt output from PSSB and binds it all into one data object
taxaBind <- function(file.path) {
  
  path.files <- list.files(file.path)
  # read in files
  list.with.each.file <- lapply(paste(file.path, list.files(file.path), sep = ''), function(y) read.delim(y, header=TRUE))
  taxa<-do.call("rbind.data.frame", list.with.each.file)
  return(taxa)
  
  
}

file.path="./taxonomy data/"###create a folder for taxa .txt downloads from PSSB
raw<-taxaBind(file.path)

length(unique(raw$Project))
length(unique(raw$Agency))
length(unique(raw$Site.Code))
length(unique(raw$WRIA.Number))
length(unique(raw$Stream.or.River))
length(unique(raw$Subbasin))

OTU<-raw

###prepare the data
OTU$Visit.Date<-format(as.Date(OTU$Visit.Date, "%m/%d/%Y"))
OTU$Year<-year(OTU$Visit.Date)

OTU<-subset(OTU, is.na(QC.Replicate.Of.Sample.Code)|QC.Replicate.Of.Sample.Code=="")##remove QC replicates

OTU$Unique<-as.logical(OTU$Unique)

##collapse to Visit.ID, because 1998-2015 samples were often three reps of 3 sq ft with different sample names for each rep
# test<-subset(ddply(OTU, .(Visit.ID, WRIA.Number, Agency, Basin, Subbasin, Stream.or.River, Project, Visit.Date, Year, Latitude, Longitude, Lab.Name, Site.Code), summarize, Samples = length(unique(Sample.Code))), Samples>1)

OTU_collapsed<-ddply(OTU, .(Visit.ID, Taxon, WRIA.Number, Agency, Basin, Subbasin, Stream.or.River, Project, Visit.Date, Year, Latitude, Longitude, Lab.Name, Site.Code), summarize, Quantity_OTU = sum(Quantity), Unique_OTU=any(Unique))

OTU_collapsed$Visit.Date<-as.Date(OTU_collapsed$Visit.Date, "%Y-%m-%d")

names(OTU)
unique(OTU[, c(29,48:69)])
OTU[is.na(OTU)]<-""
OTU_collapsed<-merge(OTU_collapsed, unique(OTU[, c(29,48:69)]), by.x="Taxon", by.y="Taxon", all.x=T)

any(is.na(OTU_collapsed$Phylum))
any(OTU_collapsed$Phylum=="")

###################################### Roll-up by broad rules ##########
KC_taxa_coarse<-OTU_collapsed

KC_taxa_coarse$OTU_COARSE<-""
names(KC_taxa_coarse)

coarse_rules<-data.frame(taxa=c("Oligochaeta", "Acari", "Gastropoda","Dytiscidae", "Simuliidae", "Chironomidae", "Trichoptera"), ranktouse=c("Subclass", "Subclass", "Family", "Family", "Family", "Family", "Genus"), rank=c("Subclass", "Subclass", "Class", "Family", "Family", "Family", "Order"))

for (j in 1:nrow(coarse_rules)){
  
  STE_rank<-coarse_rules[j, "ranktouse"]
  rank<-coarse_rules[j, "rank"]
  taxa<-coarse_rules[j, "taxa"]
  index<-which(names(KC_taxa_coarse)== STE_rank)
  rankindex<-which(names(KC_taxa_coarse)== rank)
  halt1<-which(names(KC_taxa_coarse)== "Phylum")
  halt2<-which(names(KC_taxa_coarse)== "Subspecies")
  
  for (i in halt2:halt1){
    if(i > index) {
      KC_taxa_coarse[which(KC_taxa_coarse[,rankindex]==taxa),i]<-""
    } else if (i<= index) {
      
      KC_taxa_coarse[which(KC_taxa_coarse[,rankindex]==taxa&KC_taxa_coarse$OTU_COARSE==""),"OTU_COARSE"]<-KC_taxa_coarse[which(KC_taxa_coarse[,rankindex]==taxa&KC_taxa_coarse$OTU_COARSE==""),i]
    }
    
  }
}
KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==""),"OTU_COARSE"]<-KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==""),"Taxon"]

KC_taxa_coarse[is.na(KC_taxa_coarse)]<-""


library(plyr)
##I changed the following block of code from other tests to include the original taxa names as well as the OTU names. 
##This (temporarily) prevents consolidation of OTUs within samples, but lets us match the correct original TSNs later.
##This step is to help mimic the non-random subsampling routine set up in dev version of PSSB by James Develle for testing purposes.
OTU_collapsed2<-ddply(KC_taxa_coarse, .(Taxon, Visit.ID, OTU_COARSE,Agency, WRIA.Number, Basin, 
                                        Subbasin, Stream.or.River, Project, Visit.Date, 
                                        Year, Latitude, Longitude, Lab.Name, Site.Code
), summarize, Quantity_OTU = sum(Quantity_OTU), Unique_OTU=any(Unique_OTU))

##append on correct hierarchy
###some taxa have multiple unique hierarchies because different entries got rolled up to the same level, and taxonomists have been uneven in entering in 'infraorder', 'suborder' and 'infraclass'. Need to run some loops to clean this up by selecting the most complete hierarchy available
names(KC_taxa_coarse)
new_hierarchy<-unique(KC_taxa_coarse[,c(17:39)])
countlength<-ddply(new_hierarchy, .(OTU_COARSE), summarize, count=length(OTU_COARSE))
fixthese<-countlength[countlength$count>1,]
check<-new_hierarchy[new_hierarchy$OTU_COARSE %in% fixthese$OTU_COARSE,]
for (u in 1:nrow(check)){
  check$sum[u]<-sum(check[u,]!="")
}
for (j in unique(check$OTU_COARSE)) {
  
  test<-check[check$OTU_COARSE==j,]
  summ<-colSums(test == "")
  ilen<-max(summ)
  fixlevels<-summ[which(summ<ilen& summ>0)]
  
  for (z in names(fixlevels)){
    update<-test[which(test[,z]!=""), z]
    if (length(update)>1) print("Error in taxonomy agreement")
    check[check$OTU_COARSE==j,z]<-update
  }
}

check<-check[,-which(names(check) %in% "sum")]
check<-unique(check)

new_hierarchy<-subset(new_hierarchy, ! OTU_COARSE %in% check$OTU_COARSE)
new_hierarchy<-rbind(new_hierarchy, check)
countlength<-ddply(new_hierarchy, .(OTU_COARSE), summarize, count=length(OTU_COARSE))
fixthese<-countlength[countlength$count>1,]

OTU_collapsed3<-merge(OTU_collapsed2, new_hierarchy, by.x="OTU_COARSE", by.y="OTU_COARSE", all.x=T)
any(is.na(OTU_collapsed3$Phylum))
unique(OTU_collapsed3[which(is.na(OTU_collapsed3$Phylum)),]$OTU_COARSE)
any(OTU_collapsed3$Phylum=="")


################read in PSSB attribute table, do rolling lookup between coarse taxa hierarchy and attribute table ############
setwd(here::here())
atts<-read.xlsx("current STEs and attribute.xlsx")
names(OTU_collapsed3)
missing_atts<-unique(subset(OTU_collapsed3, select=c(17:38,1)))
attribs2<-data.frame(Taxon.Name=character(), Predator=character(), Long.Lived=character(), Tolerant=character(), Intolerant=character(), Clinger=character(), OTU_COARSE=character(),  iter=numeric())

for (i in 1:ncol(missing_atts)){
  k<-(ncol(missing_atts)+1)-i
  attribs<-merge(subset(atts, select=c("Taxon.Name", "2012.Clinger", "2012.Intolerant", "2012.Long.Lived", "2012.Predator", "2012.Tolerant")), missing_atts[, c(k, 23)], by.x="Taxon.Name", by.y=names(missing_atts[k]))
  names(attribs)<-str_replace(names(attribs), "2012.", "")
  names(attribs)<-str_replace(names(attribs), ".1", "")
  attribs2<-rbind(attribs, attribs2)
  names(attribs2)<-str_replace(names(attribs2), ".1", "")
  missing_atts<-subset(missing_atts, !OTU_COARSE %in% attribs2$OTU_COARSE)
}

attribs2<-attribs2[, -1]
names(attribs2)[6]<-"Taxon"
attribs<-unique(attribs2)
any(duplicated(attribs$Taxon))
attribs[(which(duplicated(attribs$Taxon))),]
missing<-unique(KC_taxa_coarse$OTU_COARSE)[!unique(KC_taxa_coarse$OTU_COARSE) %in% attribs$Taxon] ##These taxa do not have a match in the PSSB attribute table


OTU_collapsed3<-left_join(OTU_collapsed3, attribs, by=c("OTU_COARSE"="Taxon"))
any(OTU_collapsed3$Clinger=="")
any(is.na(OTU_collapsed3$Clinger))
##Fill in attributes as "FALSE" for the taxa with no attribute matches
OTU_collapsed3[which(is.na(OTU_collapsed3$Clinger)),"Clinger"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Intolerant)),"Intolerant"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Long.Lived)),"Long.Lived"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Predator)),"Predator"]<-FALSE
OTU_collapsed3[which(is.na(OTU_collapsed3$Tolerant)),"Tolerant"]<-FALSE

##Add in the correct TSNs FOR ORIGINAL TAXA 
##For testing purposes, James changed the subsampling routine in dev version of PSSB to order by uniqueness (True to False) and then by *original* TSN (low to high) 
PSSB_taxa_list<-unique(subset(raw, select=c("Taxon.Serial.Number", "Taxon")))
PSSB_taxa_list[which(duplicated(PSSB_taxa_list$Taxon)),]
PSSB_taxa_list<-subset(PSSB_taxa_list, Taxon.Serial.Number !="-40")

OTU_collapsed3<-merge(OTU_collapsed3, PSSB_taxa_list, by.x="Taxon", by.y="Taxon", all.x=T)
any(OTU_collapsed3$Taxon.Serial.Number=="")
any(is.na(OTU_collapsed3$Taxon.Serial.Number))


names(OTU_collapsed3)[which(names(OTU_collapsed3)=="Taxon.Serial.Number")]<-"TSN"

###make sure to exclude taxa before scoring but after mapping

exlude<-read.xlsx("PSSB_exclusions.xlsx")
exlude<-subset(exlude, select=c(Taxon.Name, Excluded))

##Perform a series of rolling lookups to find a match in taxa hierarchies in the Taxon Name column of the exclusion table, starting with direct match, then from lowest hierarchy to highest

exlude2<-data.frame(Taxon.Name=character(), Excluded=character(), OTU_COARSE=character())
names(OTU_collapsed3)
test<-unique(subset(OTU_collapsed3, select=c(18:39,2)))
for (i in 1:ncol(test)){
  k<-(ncol(test)+1)-i ##work backwards through the hierarchy columns
  attribs<-merge(exlude, test[, c(k, 23)], by.x="Taxon.Name", by.y=names(test[k]))
  # names(attribs)<-str_replace(names(attribs), "2012.", "")
  names(attribs)<-str_replace(names(attribs), "OTU_COARSE.1", "OTU_COARSE")
  exlude2<-rbind(attribs, exlude2)
  test<-subset(test, !OTU_COARSE %in% exlude2$OTU_COARSE)
}

exlude2<-subset(exlude2, Excluded==T)
exlude2<-unique(exlude2)
OTU_collapsed3<-merge(OTU_collapsed3, exlude2, by="OTU_COARSE", all.x=T)
OTU_collapsed3<-subset(OTU_collapsed3, is.na(Excluded)|Excluded==F)
OTU_collapsed3<-OTU_collapsed3[,c(1:45)]


setwd(here::here())


write.csv(OTU_collapsed3, "Collapsed_Coarse_Taxa_subsampling.csv")
OTU_collapsed3<-read.csv( "Collapsed_Coarse_Taxa_subsampling.csv")


###This method of subsampling is completely non-random, and is being used for testing purposes only! It sorts taxa by *original* TSN, and by uniqueness.
##In the dev version of PSSB, James Develle temporarily changed the subsample routine to be non-random so we can check data flow through the subsampling step. 
###DO NOT USE THIS FOR GENERATING REAL SCORES!!
rarify<- function (inbug, sample.ID, abund, subsize, TSnumber)
{
  start.time = proc.time()
  outbug <- inbug
  sampid <- unique(inbug[, sample.ID])
  nsamp <- length(sampid)
  outbug[, abund] <- 0
  for (i in 1:nsamp) {
    isamp <- sampid[i]
    utils::flush.console()
    onesamp <- inbug[inbug[, sample.ID] == isamp, ]
    onesamp <- data.frame(onesamp, row.id = seq(1, dim(onesamp)[[1]]))

    samp.expand <- rep(x = onesamp[,c("row.id")], times = onesamp[, abund])

    samp.expand<-merge(samp.expand, onesamp, by.y="row.id", by.x="x")
    samp.expand<-samp.expand[,c("x", "Unique_OTU", TSnumber)]


    nbug <- length(samp.expand$x)

    samp.ex2 <- samp.expand[order(samp.expand$TSN,samp.expand$Unique_OTU,  decreasing=c( FALSE, TRUE)),]
    if (nbug > subsize) {
      subsamp <- samp.ex2[1:(subsize),]
    }
    else {
      subsamp <- samp.ex2
    }
    subcnt <- table(subsamp$x)
    newsamp <- onesamp
    newsamp[, abund] <- 0
    newsamp[match(newsamp$row.id, names(subcnt), nomatch = 0) >
              0, abund] <- as.vector(subcnt)
    outbug[outbug[, sample.ID] == isamp, abund] <- newsamp[,
                                                           abund]
  }
  elaps <- proc.time() - start.time
  cat(c("Rarify of samples complete. \n Number of samples = ",
        nsamp, "\n"))
  cat(c(" Execution time (sec) = ", elaps[1]))
  utils::flush.console()
  return(outbug)
}


KC_rarified<-rarify(inbug=OTU_collapsed3,
                    sample.ID="Visit.ID",
                    abund="Quantity_OTU",
                    subsize<-500,
                    TSnumber="TSN"
)


KC_rarified[duplicated(KC_rarified[,c("Visit.ID", "OTU_COARSE")])&KC_rarified$Site.Code=="09SOO1130",]

KC_rarified_test<-ddply(KC_rarified, .(Visit.ID, OTU_COARSE,Agency, WRIA.Number, Basin, 
                                       Subbasin, Stream.or.River, Project, Visit.Date, 
                                       Year, Latitude, Longitude, Lab.Name, Site.Code), summarize, Quantity_OTU = sum(Quantity_OTU), Unique_OTU=any(Unique_OTU))


detach("package:plyr", unload = TRUE)

str(KC_rarified)
KC_rarified$Clinger<-as.logical(KC_rarified$Clinger)
KC_rarified$Long.Lived<-as.logical(KC_rarified$Long.Lived)
KC_rarified$Intolerant<-as.logical(KC_rarified$Intolerant)
KC_rarified$Tolerant<-as.logical(KC_rarified$Tolerant)
KC_rarified$Predator<-as.logical(KC_rarified$Predator)


KC_rare_names<-names(KC_rarified)[c(2,19:45)]

KC_rarified_test<-merge(KC_rarified_test, unique(KC_rarified[,KC_rare_names]), by="OTU_COARSE")
store_rare<-KC_rarified
KC_rarified<-KC_rarified_test
KC_rarified<-subset(KC_rarified, Quantity_OTU>0)

Tot_Richness<-KC_rarified %>% dplyr::filter(Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(Total_Richness=length(OTU_COARSE))
E_Richness<-KC_rarified %>% dplyr::filter(Order=="Ephemeroptera", Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(E_Richness=length(OTU_COARSE))
P_Richness<-KC_rarified %>% dplyr::filter(Order=="Plecoptera", Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(P_Richness=length(OTU_COARSE))
T_Richness<-KC_rarified %>% dplyr::filter(Order=="Trichoptera", Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(T_Richness=length(OTU_COARSE))
Cling_Richness<-KC_rarified %>% dplyr::filter(`Clinger`==TRUE, Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(Clin_Richness=length(OTU_COARSE))
LL_Richness<-KC_rarified %>% dplyr::filter(`Long.Lived`==TRUE, Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(LL_Richness=length(OTU_COARSE))
Intol_Richness<-KC_rarified %>% dplyr::filter(`Intolerant`==TRUE, Unique_OTU==TRUE)  %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(Intol_Richness=length(OTU_COARSE))
Tot_Abund<-KC_rarified   %>% dplyr::group_by(Visit.ID) %>% dplyr::summarise(Tot_Abund=sum(Quantity_OTU)) 
Tol_Abund<-KC_rarified   %>% dplyr::group_by(Visit.ID)%>% dplyr::filter(`Tolerant`==TRUE) %>% dplyr::summarise(Tol_Abund=sum(Quantity_OTU))
Pred_Abund<-KC_rarified   %>% dplyr::group_by(Visit.ID)%>% dplyr::filter(`Predator`==TRUE) %>% dplyr::summarise(Pred_Abund=sum(Quantity_OTU))
Dom_abund<-KC_rarified   %>% dplyr::group_by(Visit.ID) %>% dplyr::arrange(desc(Quantity_OTU)) %>% dplyr::slice(1:3) %>% dplyr::summarise(Dom_Abund=sum(Quantity_OTU))

KC_results<-Reduce(function(x, y, ...) merge(x,y, all=TRUE, ...), list(Tot_Richness, E_Richness, P_Richness, T_Richness, Cling_Richness, LL_Richness, Intol_Richness, Dom_abund,
                                                                       Tol_Abund, Pred_Abund,Tot_Abund))
KC_results<-KC_results %>% mutate(Tol_Percent=Tol_Abund/Tot_Abund*100, Pred_Percent= Pred_Abund/Tot_Abund*100, Dom_Percent= Dom_Abund/Tot_Abund*100)
KC_results[is.na(KC_results)]<-0

KC_results<-KC_results %>% mutate(Tot_Richness_Score= 10 * (Total_Richness-16)/(37-16))
KC_results<-KC_results %>% mutate(E_Richness_Score= 10 * (E_Richness-1)/(8-1))
KC_results<-KC_results %>% mutate(P_Richness_Score= 10 * (P_Richness-1)/(8-1))
KC_results<-KC_results %>% mutate(T_Richness_Score= 10 * (T_Richness-1)/(9-1))
KC_results<-KC_results %>% mutate(Cling_Richness_Score= 10 * (Clin_Richness-5)/(22-5))
KC_results<-KC_results %>% mutate(LL_Richness_Score= 10 * (LL_Richness-2)/(10-2))
KC_results<-KC_results %>% mutate(Intol_Richness_Score= 10 * (Intol_Richness-0)/(7-0))
KC_results<-KC_results %>% mutate(Dom_Percent_Score= 10-(10 * (Dom_Percent-44)/(82-44)))
KC_results<-KC_results %>% mutate(Pred_Percent_Score= 10 * (Pred_Percent-1)/(21-1))
KC_results<-KC_results %>% mutate(Tol_Percent_Score= 10-(10 * (Tol_Percent-0)/(43-0)))

KC_results[which(KC_results$Tot_Richness_Score>10),"Tot_Richness_Score"]<-10
KC_results[which(KC_results$Tot_Richness_Score<0),"Tot_Richness_Score"]<-0
KC_results[which(KC_results$E_Richness_Score>10),"E_Richness_Score"]<-10
KC_results[which(KC_results$E_Richness_Score<0),"E_Richness_Score"]<-0
KC_results[which(KC_results$P_Richness_Score>10),"P_Richness_Score"]<-10
KC_results[which(KC_results$P_Richness_Score<0),"P_Richness_Score"]<-0
KC_results[which(KC_results$T_Richness_Score>10),"T_Richness_Score"]<-10
KC_results[which(KC_results$T_Richness_Score<0),"T_Richness_Score"]<-0
KC_results[which(KC_results$Cling_Richness_Score>10),"Cling_Richness_Score"]<-10
KC_results[which(KC_results$Cling_Richness_Score<0),"Cling_Richness_Score"]<-0
KC_results[which(KC_results$LL_Richness_Score>10),"LL_Richness_Score"]<-10
KC_results[which(KC_results$LL_Richness_Score<0),"LL_Richness_Score"]<-0
KC_results[which(KC_results$Intol_Richness_Score>10),"Intol_Richness_Score"]<-10
KC_results[which(KC_results$Intol_Richness_Score<0),"Intol_Richness_Score"]<-0
KC_results[which(KC_results$Dom_Percent_Score>10),"Dom_Percent_Score"]<-10
KC_results[which(KC_results$Dom_Percent_Score<0),"Dom_Percent_Score"]<-0
KC_results[which(KC_results$Pred_Percent_Score>10),"Pred_Percent_Score"]<-10
KC_results[which(KC_results$Pred_Percent_Score<0),"Pred_Percent_Score"]<-0
KC_results[which(KC_results$Tol_Percent_Score>10),"Tol_Percent_Score"]<-10
KC_results[which(KC_results$Tol_Percent_Score<0),"Tol_Percent_Score"]<-0

KC_results<-KC_results %>% mutate(Overall.Score=Tot_Richness_Score+ E_Richness_Score+P_Richness_Score+T_Richness_Score+ Cling_Richness_Score+ LL_Richness_Score+ Intol_Richness_Score+ Dom_Percent_Score+ Pred_Percent_Score+ Tol_Percent_Score)

KC_results<-left_join(KC_results, unique(subset(KC_rarified, select=c("Visit.ID","Agency", "WRIA.Number", "Basin", 
                                                                      "Subbasin", "Stream.or.River", "Project", "Visit.Date", 
                                                                      "Year", "Lab.Name", "Site.Code"))), by="Visit.ID")

KC_results<-left_join(KC_results, unique(subset(raw, select=c(Visit.ID, Latitude, Longitude))), by="Visit.ID")

write.csv(KC_results, "B-IBI_results_PSSB_Scores_nomapping_subsampling.csv")

PSSB_scores<-read.xlsx("09SOO1130_testsubsample_scores.xlsx") ##These are scores generated from PSSB with the modified non-random subsampling

test<-merge(KC_results, PSSB_scores, by=c("Visit.ID"))
ggplot(test, aes(x=Overall.Score.x, y=Overall.Score.y))+geom_point()
ggplot(test, aes(x=Total_Richness, y=Taxa.Richness.Quantity))+geom_point()
ggplot(test, aes(x=E_Richness, y=Ephemeroptera.Richness.Quantity))+geom_point()
ggplot(test, aes(x=P_Richness, y=Plecoptera.Richness.Quantity))+geom_point()
ggplot(test, aes(x=T_Richness, y=Trichoptera.Richness.Quantity))+geom_point()
ggplot(test, aes(x=Clin_Richness, y=Clinger.Richness.Quantity))+geom_point()
ggplot(test, aes(x=LL_Richness, y=`Long-Lived.Richness.Quantity`))+geom_point()
ggplot(test, aes(x=Intol_Richness, y=Intolerant.Richness.Quantity))+geom_point()
ggplot(test, aes(x=Tol_Percent, y=Tolerant.Percent.Quantity*100))+geom_point()
ggplot(test, aes(x=Pred_Percent, y=Predator.Percent.Quantity*100))+geom_point()
ggplot(test, aes(x=Dom_Percent, y=Percent.Dominant.Quantity*100))+geom_point()

test$diff<-test$Overall.Score.x-test$Overall.Score.y
write.csv(test, "subsampling_nomapping_scorecomparison.csv")
test$diff<-test$Pred_Percent-(test$Predator.Percent.Quantity*100)
issues<-subset(test, abs(diff)>5)
# issues<-issues[,c(1:26,55:78)]
# zoom<-issues[which(rownames(issues)==1480),]
# zoom<-zoom[,c(1:26,55:78)]
