###TEST Script #4: No mapping, variations of subsampling

##The purpose of this test is to evaluate different methods of PSSB subsampling routines. The "business as usual" subsampling method in PSSB vs
##a completely random subsampling method. Business as usual involves assigning individual organisms within a sample to a random number, then sorting first by the random number, and then by the uniqueness. The result is a list of taxa sorted semi-randomly: non-unique taxa always end up at the bottom of the list. The subsampling procedure then skims the first "n" organisms from the top of the list, where "n"= the desired subsample quantity. The "weighting" by uniquness was originally designed to prioritize unique taxa, and preserve sample richness.

##In this test, I use the un-rarified taxa output from the "scoring test.R" script to re-create the expected data flow, and compared the subsequent metrics to those calculated in PSSB. I found a substantial bias in % predator quantities. 

##In examining why, I found that "uniqueness" for the purposes of subsample weighting in PSSB is based only on laboratory defined uniqueness (e.g. TRUE vs "Never"). Uniqueness subsequently defined in PSSB after rollup/mapping (e.g. TRUE vs "Ignored") is not factored into the weighting. In other words, "Never" unique taxa are downweighted in subsampling, but "Ignored" taxa are not downweighted.This can result in skewing of percent based scores (See visit-id:5990 for example). In my R code attempt to replicate PSSB subsampling, but a key difference is that I consolidate identical taxa or identical OTUs prior to subsampling, so that 1 unique Cleptelmis addenda adult and 20 "never" unique Cleptelmis addenda larvae in the same sample aggregate to 21 unique Cleptelmis addenda. Taxa/OTUs that don't consolidate with any unique taxa are left as non-unique for subsampling in the R code. My conclusion is that the weighting process in PSSB creates unintentional artifacts in the data by somewhat arbitrarily skewing the species composition of the sample, preserving diversity metrics at the expense of composition-based metrics.

##As a result, I conducted a second test. I had James remove all weighting from the dev PSSB subsampling routine. I re-released WRIA_08_Survey samples collected in 2013 to use this new subsample procedure. I then compare the PSSB results with results from my own random subsampling routine in R. 
#I then came up with another method of performing random subsampling in R that is more similar to PSSBs random method, and compared the outputs.
#I also compared PSSB "business as usual" subsampling-bsaed scores with PSSB randomized subsampling-based scores.
##From these tests, I find that the PSSB scores based on the new randomized subsampling routine are comparable to scores calculated in R with random subsampling. Randomizing the subsamples in R with different methods had no effect on scores/metric quantities. PSSb scores/metric quantities calculated based on "business as usual" subsampling tended to be higher than those calculated with random subsampling, except for percent tolerant quantities, which tended to be lower in "business as usual" based scores. 


# For a third test, I randomly re-subsampled the same set of taxonomic data (n=10,368 site visits) 50 times and compared the resulting scores to get a sense of how much score "movement" occurs due to the truly random subsampling. I calculated the standard deviation of scores within individual site visits across the 50 different scores, and then removed site visits where SD= 0, because these were visits where fewer than 500 organisms were collected (and thus no subsampling occurred). I then summarized the remaining SD values (n = 7,938 site visits): the minimum SD in the dataset was 0.007, the max SD was 7.79, the median SD was 1.48, and the mean SD was 1.78. I also looked at range of scores by site visit. The minimum range of scores in the dataset was 0.05, the max range was 32.8, the median range was 6.07, and the mean range was 7.59.

##I compared those truly random subsampling-based boot-strapped results to the bootstrapped results of subsampling in the more PSSB-business-as-usual routine-- with unique taxa prioritized for subsampling. As I found in test #1 above, the method I use in R does not exactly replicate the PSSB "business as usual routine", but the results are similar enough in overall score that I felt this was an adequate comparison to use. Running the test exactly as described above but with a unique-prioritized subsampling method, I found: the minimum SD in the dataset was 0.007, the max SD was 7.21, the median SD was 1.36, and the mean SD was 1.68. I also looked at range of scores by site visit. The minimum range of scores in the dataset was 0.05, the max range was 32.7, the median range was 5.63, and the mean range was 7.17. More site visits were screened out for having SD=0 using this subsampling method (n = 7,677), because the subsampling isn't truly random, resulting in identical scores even in samples with >500 organisms. I looked at the results where I screened out only site visits that were excluded from the random bootstrapping results, for a more direct comparison. In this case (n = 7,938 site visits used), I found: the minimum SD in the dataset was 0.0, the max SD was 7.21, the median SD was 1.31, and the mean SD was 1.62. I also looked at range of scores by site visit. The minimum range of scores in the dataset was 0.0, the max range was 32.7, the median range was 5.45, and the mean range was 6.94.


##From all of this I conclude: The "business as usual" subsampling method as currently used by PSSB introduces hard-to-replicate artifacts in the metrics scoring by skewing community composition in a somewhat arbitrary way. The skewing is not very visually obvious except in %predator quantities (probably because of commonly non-unique predator taxa like "Chloroperlidae" or "perlidae" getting flagged as non-unique at the lab bench ("never"), rather than getting flagged as non-unique by the PSSB routine ("ignored")). When the PSSB subsampling routine is set to be fully random, the scores and metrics are without bias from what I would expect with truly random subsampling. Examining the variability in the scores between "business as usual" and fully random subsampling, I find that variability in overall scores is hardly different between the two. This is important because the original rationale for using semi-nonrandom-subsampling was that true richness would be better preserved in the "business as usual" approach, and thus B-IBI scores would be more stable. The results indicate that the semi-random approach does not substantially reduce variability, but it does add unintentional and arbitrary artifacts to the data, resulting in high potential for bias and difficulty in replication. I recommend PSSB switches to the fully random subsampling routine in the production site. Once this push is made, subsamples (and resulting scores) will not change in the production site until/unless site visits get re-released.


##########test 1######

library(dplyr)
library(ggplot2)
library(plotly)

OTU_collapsed3<-read.csv( "Collapsed_Coarse_Taxa.csv")## For testing truly random methods


#####the rarify function below may work more similarly to PSSB current procedure. This method prioritizes unique taxa. This was generated after discussion with James Develle.
##Note-- I discovered PSSB subsampling doesn't work quite like this.
##One issue in PSSB is that "uniqueness" for the purposes of subsample weighting is based on laboratory defined uniqueness (e.g. TRUE vs "Never"). Uniqueness subsequently defined in PSSB after rollup (e.g. TRUE vs "Ignored") is not factored into the weighting. In other words, "Never" unique taxa are downweighted in subsampling, but "Ignored" taxa are not downweighted. This can result in skewing of percent based scores.
##As a result, comparing scores based on the code below to scores based on PSSB shows a very subtle bias toward higher PSSB scores.
##The bias is especially pronounced in % predator quantities.

rarify<- function (inbug, sample.ID, abund, subsize, mySeed = NA)
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
    samp.expand<-samp.expand[,c("x", "Unique_OTU")]


    nbug <- length(samp.expand$x)
    if (!is.na(mySeed))
      set.seed(mySeed, "Mersenne-Twister",normal.kind = "Inversion",  sample.kind="Rounding")
    ranvec <- stats::runif(n = nbug)
    samp.expand<-cbind(samp.expand, ranvec)

    samp.ex2 <- samp.expand[order(samp.expand$Unique_OTU, samp.expand$ranvec, decreasing=c(TRUE, FALSE)),]
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
                    mySeed=17760704)


detach("package:plyr", unload = TRUE)

str(KC_rarified)
KC_rarified$Clinger<-as.logical(KC_rarified$Clinger)
KC_rarified$Long.Lived<-as.logical(KC_rarified$Long.Lived)
KC_rarified$Intolerant<-as.logical(KC_rarified$Intolerant)
KC_rarified$Tolerant<-as.logical(KC_rarified$Tolerant)
KC_rarified$Predator<-as.logical(KC_rarified$Predator)

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

PSSB_scores<-read.xlsx("PSSB_scores_updatedSTEnomapping.xlsx")## These are scores generated with "Business as usual" PSSB subsampling.

##compare R-based subsampled scores against PSSB subsampled scores
test<-merge(KC_results, PSSB_scores, by=c("Visit.ID"))
ggplot(test, aes(x=Overall.Score.x, y=Overall.Score.y))+geom_point()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Total_Richness, y=Taxa.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=E_Richness, y=Ephemeroptera.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=P_Richness, y=Plecoptera.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=T_Richness, y=Trichoptera.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Clin_Richness, y=Clinger.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=LL_Richness, y=`Long-Lived.Richness.Quantity`))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Intol_Richness, y=Intolerant.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Tol_Percent, y=Tolerant.Percent.Quantity*100))+geom_point()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Pred_Percent, y=Predator.Percent.Quantity*100))+geom_point()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Dom_Percent, y=Percent.Dominant.Quantity*100))+geom_point()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")


test$diff<-test$Overall.Score.x-test$Overall.Score.y
# write.csv(test, "subsampling_nomapping_scorecomparison.csv")
test$diff<-test$Pred_Percent-(test$Predator.Percent.Quantity*100)
issues<-subset(test, abs(diff)>5)

############test 2################

# ##The rarify function below is similar to how subsampling was done for the trends report-- and is truly random.
rarify<- function (inbug, sample.ID, abund, subsize, taxa, mySeed = NA)
{
  KC_rarified<-inbug[0,]
  for(visit in unique(inbug[,sample.ID])){
    set.seed(mySeed, "Mersenne-Twister",normal.kind = "Inversion",  sample.kind="Rounding")
    print(visit)
    test<-subset(inbug, get(sample.ID)==visit)
    testsample<-rep(test[,taxa], test[,abund])
    if(sum(test[,abund])>=subsize){
      subsamp<-sample(x = testsample, size = subsize,replace=F)
    }
    else {subsamp<-sample(x = testsample, size = sum(test[,abund]),replace=F)}
    subsamp<-as.data.frame(table(subsamp))
    names(subsamp)<-c(taxa,abund)
    subsamp_meta<-test[c(match(subsamp[,taxa], test[,taxa])),]
    subsamp_meta<-subset(subsamp_meta, select=-c(get(abund)))
    subsamp_comp<-merge(subsamp_meta, subsamp, by=taxa)
    KC_rarified<-rbind(subsamp_comp, KC_rarified)

  }
  return(KC_rarified)
}

KC_rarified<-rarify(inbug=OTU_collapsed3,
                    sample.ID="Visit.ID",
                    abund="Quantity_OTU",
                    subsize<-500,
                    taxa="OTU_COARSE",
                    mySeed=17760704)


detach("package:plyr", unload = TRUE)

str(KC_rarified)
KC_rarified$Clinger<-as.logical(KC_rarified$Clinger)
KC_rarified$Long.Lived<-as.logical(KC_rarified$Long.Lived)
KC_rarified$Intolerant<-as.logical(KC_rarified$Intolerant)
KC_rarified$Tolerant<-as.logical(KC_rarified$Tolerant)
KC_rarified$Predator<-as.logical(KC_rarified$Predator)

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

PSSB_scores<-read.xlsx("WRIA8_2013_PSSB_randomsubsampletest.xlsx")## These are scores generated with PSSB subsampling set to be completely random-- No weighting.
PSSB_scores<-read.xlsx("2015_PSSB_randomsubsampletest.xlsx")## These are scores generated with PSSB subsampling set to be completely random-- No weighting.

##compare R-based subsampled scores against PSSB subsampled scores
test<-merge(KC_results, PSSB_scores, by=c("Visit.ID"))
ggplot(test, aes(x=Overall.Score.x, y=Overall.Score.y))+geom_point()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Total_Richness, y=Taxa.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=E_Richness, y=Ephemeroptera.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=P_Richness, y=Plecoptera.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=T_Richness, y=Trichoptera.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Clin_Richness, y=Clinger.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=LL_Richness, y=`Long-Lived.Richness.Quantity`))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Intol_Richness, y=Intolerant.Richness.Quantity))+geom_jitter()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Tol_Percent, y=Tolerant.Percent.Quantity*100))+geom_point()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Pred_Percent, y=Predator.Percent.Quantity*100))+geom_point()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")
ggplot(test, aes(x=Dom_Percent, y=Percent.Dominant.Quantity*100))+geom_point()+geom_abline(slope=1, intercept = 0)+geom_smooth(method="lm")


###test different ways of doing random subsamples in R#####
OTU_collapsed_my_random<-KC_results

###Another way of making the subsampling totally random that is more similar to PSSB randomized method. For testing purposes.
rarify<- function (inbug, sample.ID, abund, subsize, mySeed = NA)
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
    samp.expand<-samp.expand[,c("x", "Unique_OTU")]


    nbug <- length(samp.expand$x)
    if (!is.na(mySeed))
      set.seed(mySeed, "Mersenne-Twister",normal.kind = "Inversion",  sample.kind="Rounding")
    ranvec <- stats::runif(n = nbug)
    samp.expand<-cbind(samp.expand, ranvec)

    samp.ex2 <- samp.expand[order(samp.expand$ranvec, decreasing=c(FALSE)),]
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
                    mySeed=17760704)


detach("package:plyr", unload = TRUE)

str(KC_rarified)
KC_rarified$Clinger<-as.logical(KC_rarified$Clinger)
KC_rarified$Long.Lived<-as.logical(KC_rarified$Long.Lived)
KC_rarified$Intolerant<-as.logical(KC_rarified$Intolerant)
KC_rarified$Tolerant<-as.logical(KC_rarified$Tolerant)
KC_rarified$Predator<-as.logical(KC_rarified$Predator)

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

OTU_collapsed_my_other_random<-KC_results
names(OTU_collapsed_my_random)
names(OTU_collapsed_my_other_random)
my_random_test<-merge(OTU_collapsed_my_random[c(1:26)], OTU_collapsed_my_other_random[c(1:26)], by="Visit.ID")
ggplot(my_random_test, aes(x=Overall.Score.x, y=Overall.Score.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(my_random_test, aes(x=Total_Richness.x, y=Total_Richness.y))+geom_jitter()+geom_abline(slope=1, intercept = 0)
ggplot(my_random_test, aes(x=E_Richness.x, y=E_Richness.y))+geom_jitter()+geom_abline(slope=1, intercept = 0)
ggplot(my_random_test, aes(x=P_Richness.x, y=P_Richness.y))+geom_jitter()+geom_abline(slope=1, intercept = 0)
ggplot(my_random_test, aes(x=T_Richness.x, y=T_Richness.y))+geom_jitter()+geom_abline(slope=1, intercept = 0)
ggplot(my_random_test, aes(x=Clin_Richness.x, y=Clin_Richness.y))+geom_jitter()+geom_abline(slope=1, intercept = 0)
ggplot(my_random_test, aes(x=LL_Richness.x, y=LL_Richness.y))+geom_jitter()+geom_abline(slope=1, intercept = 0)
ggplot(my_random_test, aes(x=Intol_Richness.x, y=Intol_Richness.y))+geom_jitter()+geom_abline(slope=1, intercept = 0)
ggplot(my_random_test, aes(x=Tol_Percent.x, y=Tol_Percent.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(my_random_test, aes(x=Pred_Percent.x, y=Pred_Percent.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(my_random_test, aes(x=Dom_Percent.x, y=Dom_Percent.y))+geom_point()+geom_abline(slope=1, intercept = 0)


##compare PSSB "business as usual" scores against PSSB "random" subsamples#########
PSSB_scores1<-read.xlsx("PSSB_scores_updatedSTEnomapping.xlsx")## These are scores generated with "Business as usual" PSSB subsampling.
PSSB_scores2<-read.xlsx("WRIA8_2013_PSSB_randomsubsampletest.xlsx")## These are scores generated with PSSB subsampling set to be completely random-- No weighting.
PSSB_scores3<-read.xlsx("2015_PSSB_randomsubsampletest.xlsx")## These are scores generated with PSSB subsampling set to be completely random-- No weighting.
names(PSSB_scores1)
names(PSSB_scores2)
names(PSSB_scores3)
PSSB_scores_test<-merge(PSSB_scores1[c(2,18:40)], PSSB_scores2[c(2,18:40)], by="Visit.ID")
PSSB_scores_test2<-merge(PSSB_scores1[c(2,18:40)], PSSB_scores3[c(2,18:40)], by="Visit.ID")

ggplot(PSSB_scores_test, aes(x=Overall.Score.x, y=Overall.Score.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test, aes(x=Taxa.Richness.Quantity.x, y=Taxa.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test, aes(x=Ephemeroptera.Richness.Quantity.x, y=Ephemeroptera.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test, aes(x=Plecoptera.Richness.Quantity.x, y=Plecoptera.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test, aes(x=Trichoptera.Richness.Quantity.x, y=Trichoptera.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test, aes(x=Clinger.Richness.Quantity.x, y=Clinger.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test, aes(x=`Long-Lived.Richness.Quantity.x`, y=`Long-Lived.Richness.Quantity.y`))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test, aes(x=Intolerant.Richness.Quantity.x, y=Intolerant.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test, aes(x=Tolerant.Percent.Quantity.x, y=Tolerant.Percent.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test, aes(x=Predator.Percent.Quantity.x, y=Predator.Percent.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test, aes(x=Percent.Dominant.Quantity.x, y=Percent.Dominant.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)

ggplot(PSSB_scores_test2, aes(x=Overall.Score.x, y=Overall.Score.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test2, aes(x=Taxa.Richness.Quantity.x, y=Taxa.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test2, aes(x=Ephemeroptera.Richness.Quantity.x, y=Ephemeroptera.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test2, aes(x=Plecoptera.Richness.Quantity.x, y=Plecoptera.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test2, aes(x=Trichoptera.Richness.Quantity.x, y=Trichoptera.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test2, aes(x=Clinger.Richness.Quantity.x, y=Clinger.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test2, aes(x=`Long-Lived.Richness.Quantity.x`, y=`Long-Lived.Richness.Quantity.y`))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test2, aes(x=Intolerant.Richness.Quantity.x, y=Intolerant.Richness.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test2, aes(x=Tolerant.Percent.Quantity.x, y=Tolerant.Percent.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test2, aes(x=Predator.Percent.Quantity.x, y=Predator.Percent.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)
ggplot(PSSB_scores_test2, aes(x=Percent.Dominant.Quantity.x, y=Percent.Dominant.Quantity.y))+geom_point()+geom_abline(slope=1, intercept = 0)


# p<-ggplot(combined2, aes(x=Overall.Score.x, y=Overall.Score.y))+geom_point(aes(x=Overall.Score.x, y=Overall.Score.y, color=Site.Code.x, text=Year))+theme(legend.position = "none")+geom_smooth(method="lm", se=F, color="black")+stat_poly_eq(formula=y~x, aes(label=paste(..eq.label.., ..rr.label.., sep="~~~")), parse=T)+
#   xlab("Homecooked B-IBI")+ylab("PSSB B-IBI")+geom_abline(slope=1, intercept=0)
# p
# library(plotly)
# ggplotly(p)


########test 3 #######
#####random subsample 50 times to see how much scores move from one subsample to the next

rarify<- function (inbug, sample.ID, abund, subsize, taxa, mySeed = NA)
{
  KC_rarified<-inbug[0,]
  for(visit in unique(inbug[,sample.ID])){
    set.seed(mySeed, "Mersenne-Twister",normal.kind = "Inversion",  sample.kind="Rounding")
    print(visit)
    test<-subset(inbug, get(sample.ID)==visit)
    testsample<-rep(test[,taxa], test[,abund])
    if(sum(test[,abund])>=subsize){
      subsamp<-sample(x = testsample, size = subsize,replace=F)
    }
    else {subsamp<-sample(x = testsample, size = sum(test[,abund]),replace=F)}
    subsamp<-as.data.frame(table(subsamp))
    names(subsamp)<-c(taxa,abund)
    subsamp_meta<-test[c(match(subsamp[,taxa], test[,taxa])),]
    subsamp_meta<-subset(subsamp_meta, select=-c(get(abund)))
    subsamp_comp<-merge(subsamp_meta, subsamp, by=taxa)
    KC_rarified<-rbind(subsamp_comp, KC_rarified)
    
  }
  return(KC_rarified)
}

detach("package:plyr", unload = TRUE)

KC_results_test<-as.data.frame(unique(OTU_collapsed3$Visit.ID))
names(KC_results_test)<-"Visit.ID"

for (i in 1:50) {

  KC_rarified<-rarify(inbug=OTU_collapsed3,
                      sample.ID="Visit.ID",
                      abund="Quantity_OTU",
                      subsize<-500,
                      taxa="OTU_COARSE",
                      mySeed=i)
  

  
  # str(KC_rarified)
  KC_rarified$Clinger<-as.logical(KC_rarified$Clinger)
  KC_rarified$Long.Lived<-as.logical(KC_rarified$Long.Lived)
  KC_rarified$Intolerant<-as.logical(KC_rarified$Intolerant)
  KC_rarified$Tolerant<-as.logical(KC_rarified$Tolerant)
  KC_rarified$Predator<-as.logical(KC_rarified$Predator)
  
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
  

  KC_results_test<-merge(KC_results_test, KC_results[,c("Overall.Score", "Visit.ID")], by="Visit.ID")
  names(KC_results_test)[1+i]<-paste0("Overall.Score_",i)

  
}

KC_results_test

long <- KC_results_test %>% 
  tidyr::pivot_longer(
    cols = `Overall.Score_1`:`Overall.Score_50`, 
    names_to = "Iteration",
    values_to = "Overall.Score"
  )
library(plyr)
standev<-ddply(long, .(Visit.ID), summarize, standard.dev=sd(Overall.Score), minr=min(Overall.Score), maxr=max(Overall.Score))
standev<-subset(standev, standard.dev>0)
standev$range<-standev$maxr-standev$minr
summary(standev)

#############semi-random subsample (more PSSB-business-as-usual-like) 50 times to see how much scores move

rarify<- function (inbug, sample.ID, abund, subsize, mySeed = NA)
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
    samp.expand<-samp.expand[,c("x", "Unique_OTU")]
    
    
    nbug <- length(samp.expand$x)
    if (!is.na(mySeed))
      set.seed(mySeed, "Mersenne-Twister",normal.kind = "Inversion",  sample.kind="Rounding")
    ranvec <- stats::runif(n = nbug)
    samp.expand<-cbind(samp.expand, ranvec)
    
    samp.ex2 <- samp.expand[order(samp.expand$Unique_OTU, samp.expand$ranvec, decreasing=c(TRUE, FALSE)),]
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


detach("package:plyr", unload = TRUE)

KC_results_test2<-as.data.frame(unique(OTU_collapsed3$Visit.ID))
names(KC_results_test2)<-"Visit.ID"

for (i in 1:50) {
  
  KC_rarified<-rarify(inbug=OTU_collapsed3,
                      sample.ID="Visit.ID",
                      abund="Quantity_OTU",
                      subsize<-500,
                      mySeed=i)
  
  
  
  # str(KC_rarified)
  KC_rarified$Clinger<-as.logical(KC_rarified$Clinger)
  KC_rarified$Long.Lived<-as.logical(KC_rarified$Long.Lived)
  KC_rarified$Intolerant<-as.logical(KC_rarified$Intolerant)
  KC_rarified$Tolerant<-as.logical(KC_rarified$Tolerant)
  KC_rarified$Predator<-as.logical(KC_rarified$Predator)
  
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
  
  
  KC_results_test2<-merge(KC_results_test2, KC_results[,c("Overall.Score", "Visit.ID")], by="Visit.ID")
  names(KC_results_test2)[1+i]<-paste0("Overall.Score_",i)
  
  
}

KC_results_test2

long2 <- KC_results_test2 %>% 
  tidyr::pivot_longer(
    cols = starts_with("Overall.Score"), 
    names_to = "Iteration",
    values_to = "Overall.Score"
  )
library(plyr)
standev2<-ddply(long2, .(Visit.ID), summarize, standard.dev=sd(Overall.Score), minr=min(Overall.Score), maxr=max(Overall.Score))
# standev2<-subset(standev2, standard.dev>0)
standev2<-subset(standev2, Visit.ID %in% standev$Visit.ID) ##Changed this line to be completely comparable to the random subsampling bootstrapped results. Since unique taxa are prioritized in the subsampling, there are some site visits where SD of overall score is 0 even though there are >500 organisms in the sample.
standev2$range<-standev2$maxr-standev2$minr
summary(standev2)


compare<-merge(KC_results_test[,c(1,2)], KC_results_test2[,c(1,2)], by="Visit.ID")
ggplot(compare, aes(x=Overall.Score_1.x, y=Overall.Score_1.y))+geom_point()+geom_abline(intercept=0, slope=1, color="red")+geom_smooth(method="lm", color="yellow")

PSSB_scores<-read.xlsx("PSSB_scores_updatedSTEnomapping.xlsx")
compare<-merge(KC_results_test[,c(1,2)], PSSB_scores[,c("Visit.ID","Overall.Score", "Agency")], by="Visit.ID")
ggplot(compare, aes(x=Overall.Score_1, y=Overall.Score, color=Agency))+geom_point()+geom_abline(intercept=0, slope=1, color="red")+geom_smooth(method="lm", color="yellow")+facet_wrap(.~Agency)
