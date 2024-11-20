###TEST Script #1: No mapping, no subsampling

##In this test, I examine the behavior of the updated STE feature implemented by James Develle in the dev version of PSSB.
##To perform the test, I downloaded all PSSB scores using the criteria of "Use original taxa", and use "Test: Coarse STE", WITHOUT SUBSAMPLING.
##The reason I generated PSSB scores without subsampling or mapping for this test is to better hone in sources of potential error due solely to STE rollup
##Any score mismatches would indicate an unexpected behavior in the PSSB STE rollup function.
##I replicate in my R code the expected data processing. 
##Based on this test, I did not identify any issues in the rollup.
##However, I did observe that intolerant richness scores are including taxa that are non-unique.
##I have also identified some minor rounding issues in %'s-- x.25 is rounding to x.2 in PSSB scoring, rather than to X.3. 

##An additional observation I noted is that PSSB excludes taxa after PSSB rollup occur. So an excluded taxa that rolls up to a non-excluded
##hierarchical level is therefore not excluded. Example: Eiseniella tetrahedra is never excluded from scoring, because even at fine STE,
## it gets rolled to Eiseniella, which is not on the exclusion list. 


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

exam<-subset(raw, Non.B.IBI=="True", select=c("Quantity.Subsampling", "Non.B.IBI", "Visit.ID", "Project", "Agency"))

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

# ##Append correct hierarchy to taxa
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

fine_rules<-data.frame(taxa=c("Oligochaeta", "Acari", "Gastropoda","Dytiscidae", "Simuliidae", "Chironomidae", "Trichoptera"), ranktouse=c("Genus", "Genus", "Genus", "Genus", "Genus", "Species", "Species"), rank=c("Subclass", "Subclass", "Class", "Family", "Family", "Family", "Order"))
coarse_rules<-fine_rules

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
##For fine STE, need to adjust species names in OTU column so that they include genus
KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Species),"OTU_COARSE"]<-paste0(KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Species),"Genus"]," ", KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Species),"Species"])
##For fine STE, need to adjust subgenus names in OTU column so that they include genus
KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Subgenus),"OTU_COARSE"]<-paste0(KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Subgenus),"Genus"]," (", KC_taxa_coarse[which(KC_taxa_coarse$OTU_COARSE==KC_taxa_coarse$Subgenus),"Subgenus"], ")")

library(plyr)
OTU_collapsed2<-ddply(KC_taxa_coarse, .(Visit.ID, OTU_COARSE,Agency, WRIA.Number, Basin, 
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

###make sure to exclude taxa before scoring but after mapping

exlude<-read.xlsx("PSSB_exclusions.xlsx")
exlude<-subset(exlude, select=c(Taxon.Name, Excluded))


##Perform a series of rolling lookups to find a match in taxa hierarchies in the Taxon Name column of the exclusion table, starting with direct match, then from lowest hierarchy to highest

exlude2<-data.frame(Taxon.Name=character(), Excluded=character(), OTU_COARSE=character())
names(OTU_collapsed3)
test<-unique(subset(OTU_collapsed3, select=c(17:38,1)))
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
OTU_collapsed3<-OTU_collapsed3[,c(1:43)]




setwd(here::here())
write.csv(OTU_collapsed3, "Collapsed_Coarse_Taxa.csv")
OTU_collapsed3<-read.csv( "Collapsed_Coarse_Taxa.csv")

# write.csv(OTU_collapsed3, "Collapsed_Fine_Taxa.csv")
# OTU_collapsed3<-read.csv( "Collapsed_Fine_Taxa.csv")


detach("package:plyr", unload = TRUE)
KC_rarified<-OTU_collapsed3
str(KC_rarified)
KC_rarified$Clinger<-as.logical(KC_rarified$Clinger)
KC_rarified$Long.Lived<-as.logical(KC_rarified$Long.Lived)
KC_rarified$Intolerant<-as.logical(KC_rarified$Intolerant)
KC_rarified$Tolerant<-as.logical(KC_rarified$Tolerant)
KC_rarified$Predator<-as.logical(KC_rarified$Predator)


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

##coarse scoring
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


# ####fine scoring
# KC_results<-KC_results %>% mutate(Tot_Richness_Score= 10 * (Total_Richness-27)/(56-27))
# KC_results<-KC_results %>% mutate(E_Richness_Score= 10 * (E_Richness-1)/(8-1))
# KC_results<-KC_results %>% mutate(P_Richness_Score= 10 * (P_Richness-1)/(8-1))
# KC_results<-KC_results %>% mutate(T_Richness_Score= 10 * (T_Richness-1)/(9-1))
# KC_results<-KC_results %>% mutate(Cling_Richness_Score= 10 * (Clin_Richness-7)/(24-7))
# KC_results<-KC_results %>% mutate(LL_Richness_Score= 10 * (LL_Richness-2)/(10-2))
# KC_results<-KC_results %>% mutate(Intol_Richness_Score= 10 * (Intol_Richness-0)/(7-0))
# KC_results<-KC_results %>% mutate(Dom_Percent_Score= 10-(10 * (Dom_Percent-32)/(69-32)))
# KC_results<-KC_results %>% mutate(Pred_Percent_Score= 10 * (Pred_Percent-1)/(21-1))
# KC_results<-KC_results %>% mutate(Tol_Percent_Score= 10-(10 * (Tol_Percent-0)/(43-0)))

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

write.csv(KC_results, "B-IBI_results_PSSB_Scores_nomapping_nosubsampling.csv")
# write.csv(KC_results, "B-IBI_results_PSSB_Scores_nomapping_nosubsampling_fine.csv")

PSSB_scores<-read.xlsx("PSSB_scores_updatedSTEnomapping_nosubsampling.xlsx")
# PSSB_scores<-read.xlsx("PSSB_scores_updatedSTEnomapping_nosubsampling_fine.xlsx")

test<-merge(KC_results, PSSB_scores, by=c("Visit.ID"))
ggplot(test, aes(x=Overall.Score.x, y=Overall.Score.y))+geom_point()

test$diff<-test$Overall.Score.x-test$Overall.Score.y
write.csv(test, "nosubsampling_nomapping_scorecomparison_fine.csv")
issues<-subset(test, abs(diff)>1)
issues<-issues[,c(1:26,55:78)]
zoom<-issues[which(rownames(issues)==1480),]
zoom<-zoom[,c(1:26,55:78)]
