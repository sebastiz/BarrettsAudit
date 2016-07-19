library(stringr)
library(ggplot2)
library(lattice)
library(compare)
library(XLConnect)
library(grid)
library (reshape2)
library(plyr)
library(gtools)
library(gplots)
Endo<-read.csv("/Users/sebastianzeki/Dropbox/Work/Medical/Clinical/Gastro/BarrettsAudit/AllBarrettsPatients.txt",sep="\t", header=TRUE,quote = "", row.names = NULL, stringsAsFactors = FALSE)
RFA<-read.csv("/Users/sebastianzeki/Dropbox/Work/Medical/Clinical/Gastro/BarrettsAudit/Barr_All_RFAAndAPC.txt",sep="\t", header=TRUE,quote = "", row.names = NULL, stringsAsFactors = FALSE)
Endo<-rbind.fill(Endo, RFA)
detach("package:plyr", unload=TRUE)
library(dplyr)

#Date of Interest
DOI="2013-08-01"

Endo$VisitDate<-gsub("00:00:00 BST ","",Endo$VisitDate)
Endo$VisitDate<-gsub("00:00:00 GMT ","",Endo$VisitDate)
Endo$VisitDate<-gsub("^[A-Za-z]{3}","",Endo$VisitDate)
Endo$VisitDate <- as.Date(Endo$VisitDate, " %b %d %Y")
Endo$VisitDate<-format(Endo$VisitDate, format = "%d/%m/%y")
Endo$VisitDate<-as.Date(Endo$VisitDate, "%d/%m/%y")

#---------------------------------FILTERING---------------------------------------------------------------------------------------------------
#Make sure the negatives removed
EndoSubset <- Endo[grepl("(?<!egative for |No evidence of |No |either |or )[Bb]arret", Endo$ERDIAGNOSISSTR,perl=TRUE), ]
EndoSubset <- EndoSubset[grepl("(?<!egative for |No evidence of |No |either |or )[Bb]arret", EndoSubset$ERFINDINGSSTR,perl=TRUE), ]
#This removes reports with colon in them as they are giving false positive dysplastics
EndoSubset<-EndoSubset[!grepl("[Cc]olon",EndoSubset$Diagnosis,perl=TRUE),]

#Get the Prague scoring sorted############SORT THIS OUT SO ALL THE SCORES INCLUDED REGARDLESS OF SOURCE OF INFO
#EndoSubset$CStage <-str_extract (EndoSubset$ERFINDINGSSTR,"[1-3]{1}(\\s)*?[cC][mM]")

EndoSubset$CStage <- str_extract(EndoSubset$ERFINDINGSSTR, "(C(\\s|=)*\\d+)")


EndoSubset$CStage<-as.numeric(gsub("C","",EndoSubset$CStage))
EndoSubset$MStage <- str_extract(EndoSubset$ERFINDINGSSTR, "(M(\\s|=)*\\d+)")
EndoSubset$MStage<-as.numeric(gsub("M","",EndoSubset$MStage))

#Get the worst pathology for that sample
EndoSubset$IMorNoIM <- ifelse(grepl("[Ss][Mm]2",EndoSubset$Diagnosis,perl=TRUE),"SM2",
                              ifelse(grepl("[Ss][Mm]1",EndoSubset$Diagnosis,perl=TRUE)|grepl("[Ss][Mm]1",EndoSubset$Histology,perl=TRUE),"SM1",
                                     ifelse(grepl("T1b",EndoSubset$Diagnosis,perl=TRUE)|grepl("T1b",EndoSubset$Histology,perl=TRUE),"T1b_Unspec",
                                            ifelse(grepl("T1a|ntramucosal",EndoSubset$Diagnosis,perl=TRUE)|grepl("T1a|ntramucosal",EndoSubset$Histology,perl=TRUE),"T1a",
                                                   ifelse(grepl("[Hh]igh grade ",EndoSubset$Diagnosis,perl=TRUE)|grepl("[Hh]igh grade ",EndoSubset$Histology,perl=TRUE),"HGD",
                                                          ifelse(grepl("[Ll]ow ",EndoSubset$Diagnosis,perl=TRUE)|grepl("[Ll]ow ",EndoSubset$Histology,perl=TRUE),"LGD",
                                                                 ifelse(grepl("[Ii]ndef",EndoSubset$Diagnosis,perl=TRUE)|grepl("[Ii]ndef",EndoSubset$Histology,perl=TRUE),"IGD",
                                                                        ifelse(grepl("(?<!egative for |No evidence of |[Nn]o |either |or |and )[Ii]ntestinal [Mm]etaplasia",EndoSubset$Histology,perl=TRUE)|
                                                                                 grepl("(?<!egative for |No evidence of |[Nn]o |either |or |and )[Ii]ntestinal [Mm]etaplasia",EndoSubset$Histology,perl=TRUE),"IM",
                                                                               ifelse(grepl("(?<!egative for |No evidence of |[Nn]o |either |or |and )[Ii]ntestinal [Mm]etaplasia",EndoSubset$Histology,perl=TRUE)|
                                                                                        grepl("(?<!egative for |No evidence of |[Nn]o |either |or |and )[Ii]ntestinal [Mm]etaplasia",EndoSubset$Histology,perl=TRUE),"IM",
                                                                                      "No_IM")))))))))
#Get all the EVENTS in:
EndoSubset$EVENT <- ifelse(grepl("HALO|RFA|APC", EndoSubset$ERPROCEDUREPERFORMED), "RFA",
                           ifelse(grepl("HALO|APC", EndoSubset$ERDIAGNOSISSTR), "RFA",
                                  ifelse(grepl("HALO|APC", EndoSubset$ERFINDINGSSTR), "RFA",
                                         ifelse(grepl("[Ee][Mm][Rr]|[Ee]ndoscopic [Mm]ucosal [Rr]esection|ndomucosal", EndoSubset$Diagnosis), "EMR",
                                                ifelse(grepl("[Ee][Mm][Rr]|[Ee]ndoscopic [Mm]ucosal [Rr]esection|ndomucosal", EndoSubset$NatureOfSpec), "EMR",
                                                       "nothing")))))

#Do the follow-up groupings
EndoSubset$MStage<-ifelse(grepl("(C(\\s|=)*\\d+)",EndoSubset$ERFINDINGSSTR,perl=TRUE),str_extract(EndoSubset$ERFINDINGSSTR, "(M(\\s|=)*\\d+)"),
                          ifelse(grepl("[Ss]hort|[Tt]iny|[Tt]ongue|[Ss]mall",EndoSubset$ERFINDINGSSTR,perl=TRUE),"0",
                                 ifelse(grepl("[1-3]{1}(\\s)*?[cC][mM]",EndoSubset$ERFINDINGSSTR,perl=TRUE),
                                        str_extract(EndoSubset$ERFINDINGSSTR, "[1-3]{1}(\\s)*?[cC][mM]"),NA)))

EndoSubset$MStage<-gsub("M|cm|=","",EndoSubset$MStage)

EndoSubset$MStage<-as.integer(EndoSubset$MStage)
EndoSubset$FU_Group<-ifelse(EndoSubset$IMorNoIM=="No_IM" &EndoSubset$MStage<3,"Rule1",
                            ifelse(EndoSubset$IMorNoIM=="IM" &EndoSubset$MStage<3,"Rule2",
                                   ifelse(EndoSubset$MStage>=3,"Rule3","NoRules")))




#Get all the endoscopists with >50 scopes and have a new column called expert in the main and call them expert vs NonExpert

------------------------------------------------------------------------------------------------------------------------------------------------------
  ------------------------------------------------------------------------------------------------------------------------------------------------------
  ------------------------------------------------------------#For the surveillance population------------------------------------------------------------
----------------------------------------#QUESTION- WHAT IS THE AVERAGE BIOPSY SIZE FOR SURVEILLANCE----------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------
  -----------------------------------------------------------------------------------------------------------------------------------------------------
  ------------------------------------------------------------------------------------------------------------------------------------------------------
  ------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #Make sure therapeutics removed
  EndoSubsetBx <- EndoSubset[!grepl("Therapeutic", EndoSubset$ERPROCEDUREPERFORMED), ]
EndoSubsetBx <- EndoSubsetBx[!grepl("EMR|ucosal resection", EndoSubsetBx$Diagnosis), ]
EndoSubsetBx <- EndoSubsetBx[grepl("astroscopy", EndoSubsetBx$ERPROCEDUREPERFORMED), ]

EndoSubsetBx<-EndoSubsetBx[grepl("[Ss]urveill|[Bb]arrett",EndoSubsetBx$INDICATIONS,perl=TRUE),]

EndoSubsetBx$MeasurementLargest<-as.integer(as.character(EndoSubsetBx$MeasurementLargest))
EndoSubsetBx$NumberBx<-as.integer(as.character(EndoSubsetBx$NumberBx))
Audit1 <- data.frame(EndoSubsetBx$ENDOSCOPIST,EndoSubsetBx$MeasurementLargest)
GroupedByEndoscopist <- Audit1 %>% group_by(EndoSubsetBx.ENDOSCOPIST) %>% na.omit() %>% summarise(avg = mean(EndoSubsetBx.MeasurementLargest))
Numbiop<-GroupedByEndoscopist


Experts<-data.frame(table(EndoSubsetOnSurveil$ENDOSCOPIST))
Experts$Status<-ifelse(Experts$Freq>50,"Expert","Non-Expert")
colnames(Experts) <- c("ENDOSCOPIST","Freq","Status") 
colnames(Numbiop) <- c("ENDOSCOPIST","avg") 
Numbiop<-merge(Experts,Numbiop,by=c("ENDOSCOPIST"))



ggplot(Numbiop) + geom_violin(aes(x=Status, y=avg, fill="red"))+
  geom_boxplot(aes(x=Status, y=avg),width=0.1,outlier.shape = NA)+
  labs(title="Average of largest biopsy size for Barrett's Surveillance list: \nExperts vs non-Experts", x="",y="") +
  theme(plot.margin = unit(c(0,0,0,0), "lines"))+
  theme(legend.position="top")+
  xlab("Experts                                                         Non-Expert") + 
  ylab("Shortfall(Obs-expec number Bx)") +
  theme(axis.text.x=element_text(angle=-90,size=3)) +
  theme(axis.text.y=element_text(size=30)) +
  theme(title=element_text(size=30))+
  theme(legend.position="top")




ggplot() + 
  geom_point(aes(Numbiop$ENDOSCOPIST, GroupedByEndoscopist$avg),size= 9) +  
  geom_point(cex=2) +
  labs(title=" Average of largest biopsy size for Barrett's \n surveillance OGD by endoscopist", x="",y="") +
  theme(plot.margin = unit(c(0,0,0,0), "lines"))+
  theme(legend.position="top")+
  xlab("Endoscopist") + 
  ylab("Size in cubic mm") +
  theme(axis.text.x=element_text(angle=-90,size=13)) +
  theme(axis.text.y=element_text(size=30)) +
  theme(title=element_text(size=30))+
  theme(legend.position="top")



-----------------------------------------------------------------------------------------------------------------------
  ----------------------------------#QUESTION- ARE WE TAKING ENOUGH BIOPSIES ON SURVEILLANCE?----------------------------------
-----------------------------------------------------------------------------------------------------------------------
  -----------------------------------------------------------------------------------------------------------------------
  #Filtered for the surveillance population
  
  #On surveillance vs not on surveillance of endoscopies
  EndoSubsetOnSurveil<-EndoSubset[grepl("[Ss]urveill|[Bb]arrett",EndoSubset$INDICATIONS,perl=TRUE),]
EndoSubsetOnSurveil$NumberBx<-as.integer(EndoSubsetOnSurveil$NumberBx)
EndoSubsetOnSurveil$MeasurementLargest<-as.integer(EndoSubsetOnSurveil$MeasurementLargest)


EndoSubsetOnSurveil$TimeZone<-ifelse(EndoSubsetOnSurveil$VisitDate<DOI,"Pre",
                                     ifelse(EndoSubsetOnSurveil$VisitDate >DOI,"Post","NA"))

EndoSubsetOffSurveil<-EndoSubset[!grepl("[Ss]urveill|[Bb]arrett",EndoSubset$INDICATIONS,perl=TRUE),]
EndoSubsetOnSurveilPreDOIShortfall<-subset(EndoSubsetOnSurveil,VisitDate <DOI)
EndoSubsetOnSurveilPostDOIShortfall<-subset(EndoSubsetOnSurveil,VisitDate >=DOI)


#For Pre DOI
GroupedByEndoscopy <- EndoSubsetOnSurveilPreDOIShortfall %>% filter(!is.na(CStage),!is.na(NumberBx))%>%group_by(Endoscopy_ID,ENDOSCOPIST) %>%summarise(Sum = sum(NumberBx),AvgC = mean(CStage))
GroupedByEndoscopy$ExpectedNumber<-(GroupedByEndoscopy$AvgC+1)*2
GroupedByEndoscopy$Difference<-GroupedByEndoscopy$Sum-GroupedByEndoscopy$ExpectedNumber

#Now group the difference by endoscopist
BxShortfallPre <- GroupedByEndoscopy %>% group_by(ENDOSCOPIST) %>% summarise(MeanDiff = mean(Difference))

#e) Then show shortfall of number of biopsies on a graph
ggplot() + 
  geom_point(aes(BxShortfallPre$ENDOSCOPIST, BxShortfallPre$MeanDiff),size= 9) +  
  geom_point(cex=2) +
  labs(title="Shortfall number of biopsies on Barrett's Surveillance list", x="",y="") +
  theme(plot.margin = unit(c(0,0,0,0), "lines"))+
  theme(legend.position="top")+
  xlab("Endoscopist") + 
  ylab("Shortfall(Obs-expec number Bx)") +
  theme(axis.text.x=element_text(angle=-90,size=13)) +
  theme(axis.text.y=element_text(size=30)) +
  theme(title=element_text(size=30))+
  theme(legend.position="top")


#For Post DOI

GroupedByEndoscopy <- EndoSubsetOnSurveilPostDOIShortfall %>% filter(!is.na(CStage),!is.na(NumberBx))%>%group_by(Endoscopy_ID,ENDOSCOPIST,TimeZone) %>%summarise(Sum = sum(NumberBx),AvgC = mean(CStage))
GroupedByEndoscopy$ExpectedNumber<-(GroupedByEndoscopy$AvgC+1)*2
GroupedByEndoscopy$Difference<-GroupedByEndoscopy$Sum-GroupedByEndoscopy$ExpectedNumber

Experts<-data.frame(table(EndoSubsetOnSurveil$ENDOSCOPIST))
Experts$Status<-ifelse(Experts$Freq>50,"Expert","Non-Expert")
colnames(Experts) <- c("ENDOSCOPIST","Freq","Status") 
BxShortfall<-merge(Experts,BxShortfall,by=c("ENDOSCOPIST"))

#Now group the difference by endoscopist
BxShortfallPost <- GroupedByEndoscopy %>% group_by(ENDOSCOPIST) %>% summarise(MeanDiff = mean(Difference))

#e) Then show shortfall of number of biopsies on a graph
ggplot() + 
  geom_point(aes(BxShortfallPost$ENDOSCOPIST, BxShortfallPost$MeanDiff),size= 9,colour="red") + 
  geom_point(aes(BxShortfallPre$ENDOSCOPIST, BxShortfallPre$MeanDiff),size= 9) + 
  geom_point(cex=2) +
  labs(title="Shortfall number of biopsies on Barrett's Surveillance list", x="",y="") +
  theme(plot.margin = unit(c(0,0,0,0), "lines"))+
  theme(legend.position="top")+
  xlab("Endoscopist") + 
  ylab("Shortfall(Obs-expec number Bx)") +
  theme(axis.text.x=element_text(angle=-90,size=3)) +
  theme(axis.text.y=element_text(size=30)) +
  theme(title=element_text(size=30))+
  theme(legend.position="top")




#The comparative boxplot:

GroupedByEndoscopy <- EndoSubsetOnSurveil %>% filter(!is.na(CStage),!is.na(NumberBx))%>%group_by(Endoscopy_ID,ENDOSCOPIST,TimeZone) %>%summarise(Sum = sum(NumberBx),AvgC = mean(CStage))
GroupedByEndoscopy$ExpectedNumber<-(GroupedByEndoscopy$AvgC+1)*2
GroupedByEndoscopy$Difference<-GroupedByEndoscopy$Sum-GroupedByEndoscopy$ExpectedNumber

#Now group the difference by endoscopist
BxShortfall <- GroupedByEndoscopy %>% group_by(ENDOSCOPIST,TimeZone) %>% summarise(MeanDiff = mean(Difference))

ggplot(BxShortfall) + geom_violin(aes(x=TimeZone, y=MeanDiff, fill="red"))+
  geom_boxplot(aes(x=TimeZone, y=MeanDiff),width=0.1)+
  labs(title="Shortfall number of biopsies on Barrett's Surveillance list", x="",y="") +
  theme(plot.margin = unit(c(0,0,0,0), "lines"))+
  theme(legend.position="top")+
  xlab("Post                                                            Pre") + 
  ylab("Shortfall(Obs-expec number Bx)") +
  theme(axis.text.x=element_text(angle=-90,size=3)) +
  theme(axis.text.y=element_text(size=30)) +
  theme(title=element_text(size=30))+
  theme(legend.position="top")

#Graph to show that the more you do, the better the number of biopsies you take
BxShortfall2<-subset(BxShortfall,BxShortfall$Freq>5)
ggplot(BxShortfall2) + geom_point(aes(x=Freq, y=MeanDiff, fill="red",group=1),size=10)+
  labs(title="Shortfall by number of surveillance endoscopies \ndone for each endoscopist", x="",y="") +
  theme(plot.margin = unit(c(0,0,0,0), "lines"))+
  theme(legend.position="top")+
  xlab("Number of surveillance endoscopies done") + 
  ylab("Shortfall(Obs-expec number Bx)") +
  theme(axis.text.x=element_text(angle=-90,size=3)) +
  theme(axis.text.y=element_text(size=30)) +
  theme(title=element_text(size=30))+
  theme(legend.position="top")

EndoSubsetOnSurveilExp<-subset(EndoSubsetOnSurveil,EndoSubsetOnSurveil$TimeZone=="Post")
#Subset to look at experts vs non experts in the POST DOI group
GroupedByEndoscopy <- EndoSubsetOnSurveilExp %>% filter(!is.na(CStage),!is.na(NumberBx))%>%group_by(Endoscopy_ID,ENDOSCOPIST,TimeZone) %>%summarise(Sum = sum(NumberBx),AvgC = mean(CStage))
GroupedByEndoscopy$ExpectedNumber<-(GroupedByEndoscopy$AvgC+1)*2
GroupedByEndoscopy$Difference<-GroupedByEndoscopy$Sum-GroupedByEndoscopy$ExpectedNumber

#Now group the difference by endoscopist
BxShortfall <- GroupedByEndoscopy %>% group_by(ENDOSCOPIST,TimeZone) %>% summarise(MeanDiff = mean(Difference))


Experts<-data.frame(table(EndoSubsetOnSurveil$ENDOSCOPIST))
Experts$Status<-ifelse(Experts$Freq>50,"Expert","Non-Expert")
colnames(Experts) <- c("ENDOSCOPIST","Freq","Status") 
BxShortfall<-merge(Experts,BxShortfall,by=c("ENDOSCOPIST"))


ggplot(BxShortfall) + geom_violin(aes(x=Status, y=MeanDiff, fill="red"))+
  geom_boxplot(aes(x=Status, y=MeanDiff),width=0.1)+
  labs(title="Shortfall number of biopsies on Barrett's Surveillance list\nExpert vs Non-expert", x="",y="") +
  theme(plot.margin = unit(c(0,0,0,0), "lines"))+
  theme(legend.position="top")+
  xlab("Expert                                                            Non-Expert") + 
  ylab("Shortfall(Obs-expec number Bx)") +
  theme(axis.text.x=element_text(angle=-90,size=3)) +
  theme(axis.text.y=element_text(size=30)) +
  theme(title=element_text(size=30))+
  theme(legend.position="top")



#-----------------------------------------------------------------------------------------------------------------------
#----------------------------------#QUESTION- ARE WE ENROLLING PATIENTS FROM INDEX ENDOSCOPY?-------------------------
#-------------------------------------------#THIS HAS TO BE RESTRICTED TO PRAGUE PATIENTS-----------------------------
#-----------------------------------------------------------------------------------------------------------------------
#  ---------------------------------------------------------------------------------------------------------------------



#Need to remove NumberBx column and Measurement column to prevent duplicates here:
EndoSubsetOnSurveil$NumberBx<-NULL
EndoSubsetOnSurveil$MeasurementLargest<-NULL
EndoSubsetOnSurveil<-unique(EndoSubsetOnSurveil)
EndoSubsetOffSurveil$NumberBx<-NULL
EndoSubsetOffSurveil$MeasurementLargest<-NULL
EndoSubsetOffSurveil<-unique(EndoSubsetOffSurveil)

#The pre or post  DOI subset of endoscopies
EndoSubsetOnSurveilPreDOI<-subset(EndoSubsetOnSurveil,VisitDate <=DOI)
EndoSubsetOffSurveilPreDOI<-subset(EndoSubsetOffSurveil,VisitDate <= DOI)
EndoSubsetOffSurveilPostDOI<-subset(EndoSubsetOffSurveil,VisitDate >DOI)
EndoSubsetOnSurveilPostDOI<-subset(EndoSubsetOnSurveil,VisitDate > DOI)


########--------Enrollment number work-up-----------#################
#----------------------------------------------------------------------------------------------------
#  ----------------------------------------------------------------------------------------------------
#  ----------------------------------------------------------------------------------------------------
EndoSubsetOnSurveilPostDOIa<-EndoSubsetOnSurveilPostDOI[1:2]
EndoSubsetOffSurveilPreDOIa<-EndoSubsetOffSurveilPreDOI[1:2]
EndoSubsetOnSurveilPostDOIHospNum<-data.frame(unique(EndoSubsetOnSurveilPostDOIa$HospNum_Id))
EndoSubsetOffSurveilPreDOIHospNum<-data.frame(unique(EndoSubsetOffSurveilPreDOIa$HospNum_Id))

#How many of these unique hospital numbers then pop up in the post DOI as surveillance endoscopy
names(EndoSubsetOnSurveilPostDOIHospNum)<-c("x")
names(EndoSubsetOffSurveilPreDOIHospNum)<-c("x")
PreAndPost<-rbind(EndoSubsetOffSurveilPreDOIHospNum,EndoSubsetOnSurveilPostDOIHospNum)

#Total number of patients enrolled after Index endoscopy
EnrolledAfterIndexEndoscopy<-data.frame(PreAndPost[duplicated(PreAndPost),])

names(EnrolledAfterIndexEndoscopy)<-c("x")
EnrolledAfterIndexEndoscopy$x<-as.character(EnrolledAfterIndexEndoscopy$x)
nrow(EnrolledAfterIndexEndoscopy)

#Proportions
Rule1Prop<-round(nrow(subset(HowManyNonSurveillBarrGetFollowedBySurveill,HowManyNonSurveillBarrGetFollowedBySurveill$FU_Group=="Rule1"))/nrow(HowManyNonSurveillBarrGetFollowedBySurveill),2)
Rule2Prop<-round(nrow(subset(HowManyNonSurveillBarrGetFollowedBySurveill,HowManyNonSurveillBarrGetFollowedBySurveill$FU_Group=="Rule2"))/nrow(HowManyNonSurveillBarrGetFollowedBySurveill),2)
Rule3Prop<-round(nrow(subset(HowManyNonSurveillBarrGetFollowedBySurveill,HowManyNonSurveillBarrGetFollowedBySurveill$FU_Group=="Rule3"))/nrow(HowManyNonSurveillBarrGetFollowedBySurveill),2)

#Who were enrolled- DO BY PATIENT
HowManyNonSurveillBarrGetFollowedBySurveill<-EndoSubsetOffSurveilPreDOI[EndoSubsetOffSurveilPreDOI$HospNum_Id %in% EnrolledAfterIndexEndoscopy$x,]

#Now get a list of the patients who were not enrolled to find out what they were-DO BY PATIENT
HowManyNonSurveillBarrGetDoNotGetEnrolled<-EndoSubsetOffSurveilPreDOI[!(EndoSubsetOffSurveilPreDOI$HospNum_Id %in% EnrolledAfterIndexEndoscopy$x),]

#---------INVESTIGATING THOSE WHO WERE NOT ENROLLED AFTER INDEX ENDOSCOPY#---------#---------#---------#---------#---------
#GROUPING THE INDEX ENDOSCOPY RESULTS

Rule1NotEnrolled<-nrow(subset(HowManyNonSurveillBarrGetDoNotGetEnrolled,HowManyNonSurveillBarrGetDoNotGetEnrolled$FU_Group=="Rule1"))
Rule2NotEnrolled<-nrow(subset(HowManyNonSurveillBarrGetDoNotGetEnrolled,HowManyNonSurveillBarrGetDoNotGetEnrolled$FU_Group=="Rule2"))
Rule3NotEnrolled<-nrow(subset(HowManyNonSurveillBarrGetDoNotGetEnrolled,HowManyNonSurveillBarrGetDoNotGetEnrolled$FU_Group=="Rule3"))
Rule3<-subset(HowManyNonSurveillBarrGetDoNotGetEnrolled,HowManyNonSurveillBarrGetDoNotGetEnrolled$FU_Group=="Rule3")

#Why didn't they get enrolled? DO for Rule 3 only as thats the only one that can be done reliably

library(gplots)
GroupedByRecallReason <- HowManyNonSurveillBarrGetDoNotGetEnrolled %>% group_by(HowManyNonSurveillBarrGetDoNotGetEnrolled$ERRECALLREASON1) %>% summarise(freq = freq(p.MeasurementLargest))
par(oma = c(14, 0, 4, 0))
barplot2(table(Rule3$ERRECALLREASON1),las=3,
         cex.lab = 1.0,cex.axis=2.5,cex.main = 2.5,cex.names=1.5,
         space=1,main = "Recall decision in patients endoscoped pre DOI with Barretts\n but not followed up on surveillance",angle = 45)

barplot(table(EndoSubsetOnSurveilPreDOI$MStage),main="M Stage Patient Frequency\nin the pre-DOI group",xlab="M Stage",ylab="Freq",cex.lab = 2.0,cex.axis=2.5,cex.main = 2.5,cex.names=2.5,)
barplot(table(EndoSubsetOnSurveilPreDOI$FU_Group),xlab = "Group", ylab = "Number of endoscopies",cex.lab = 2.0,cex.axis=2.5,cex.main = 2.5,cex.names=2.5,main = "Number of endoscopies in the pre-DOI group\nfor each rule group")


#Was there a documented decision for taking the patient off repeat endoscopy?

#Numbers for the DiagrammeR plot: This is initially by endoscopy
#PreDOI
#1. All endoscopies- dataset is EndoSubsetOffSurveilPreDOI for all the endoscopies
NumEndoSubsetOffSurveilPreDOI<-nrow(EndoSubsetOffSurveilPreDOI)
#2. Number lost as Prague -ve
EndoSubsetOffSurveilPreDOI_withPrague<-subset(EndoSubsetOffSurveilPreDOI,!(is.na(EndoSubsetOffSurveilPreDOI$MStage)))
NumEndoSubsetOffSurveilPreDOI_withPrague<-nrow(EndoSubsetOffSurveilPreDOI_withPrague)
#3. Classifiable as Rule
EndoSubsetOffSurveilPreDOI_withPragueNoRules<-subset(EndoSubsetOffSurveilPreDOI,(EndoSubsetOffSurveilPreDOI$FU_Group=="NoRules"))
NumEndoSubsetOffSurveilPreDOI_withPragueNoRules<-nrow(EndoSubsetOffSurveilPreDOI_withPragueNoRules)
#4. Number with Rule 1
EndoSubsetOffSurveilPreDOI_withPrague_Rule1<-subset(EndoSubsetOffSurveilPreDOI,(EndoSubsetOffSurveilPreDOI$FU_Group=="Rule1"))
NumEndoSubsetOffSurveilPreDOI_withPrague_Rule1<-nrow(EndoSubsetOffSurveilPreDOI_withPrague_Rule1)
#5. Number with Rule 2
EndoSubsetOffSurveilPreDOI_withPrague_Rule2<-subset(EndoSubsetOffSurveilPreDOI,(EndoSubsetOffSurveilPreDOI$FU_Group=="Rule2"))
NumEndoSubsetOffSurveilPreDOI_withPrague_Rule2<-nrow(EndoSubsetOffSurveilPreDOI_withPrague_Rule2)
#6. Number with Rule 3
EndoSubsetOffSurveilPreDOI_withPrague_Rule3<-subset(EndoSubsetOffSurveilPreDOI,(EndoSubsetOffSurveilPreDOI$FU_Group=="Rule3"))
NumEndoSubsetOffSurveilPreDOI_withPrague_Rule3<-nrow(EndoSubsetOffSurveilPreDOI_withPrague_Rule3)
#7. Percentage followed up
PercentFU_Rule3OffSuveil<-(nrow(Rule3)/nrow(HowManyNonSurveillBarrGetFollowedBySurveill))*100
#8. Reason for no follow-up


library(magrittr)
library(DiagrammeR)

AllIndexPreLabel<-paste("All Index Pre DOI",NumEndoSubsetOffSurveilPreDOI, sep = ": ")
nodes <- create_nodes(nodes = c(AllIndexPreLabel, "Index \nWith Prague", "Rule1","Rule2","Rule3","Rule 1 \nF-up","Rule 2 \nF-up","Rule 3 \nF-up"),                     
                      label = TRUE,
                      fontsize = 55,
                      fontcolour = "White",
                      type = "lower",
                      style = "filled",
                      color = "aqua",
                      shape = c("circle"),
                      x = c(0,0,-80,0,80,-80,0,80),
                      y = c(280,-80,-250,-250,-250,-400,-400,-400))

edges <- create_edges(from = c(AllIndexPreLabel, "Index \nWith Prague", "Index \nWith Prague","Index \nWith Prague","Rule1","Rule2","Rule3"), 
                      to = c("Index \nWith Prague", "Rule1","Rule2","Rule3", "Rule 1 \nF-up","Rule 2 \nF-up","Rule 3 \nF-up"),
                      rel = c(NumEndoSubsetOffSurveilPreDOI_withPrague, NumEndoSubsetOffSurveilPreDOI_withPrague_Rule1, NumEndoSubsetOffSurveilPreDOI_withPrague_Rule2, 
                              NumEndoSubsetOffSurveilPreDOI_withPrague_Rule3, Rule1Prop,Rule2Prop,Rule3Prop),
                      arrowhead = rep("normal", 6),
                      # color = c("red", "red", "red", "red", "red", "red"),
                      length = c(200,200,50,50,200,200),
                      fontsize = 55,
                      width=c(NumEndoSubsetOnSurveilPreDOI_withPrague/10,NumEndoSubsetOnSurveilPreDOI_withPrague_Rule1/10,NumEndoSubsetOnSurveilPreDOI_withPrague_Rule2/10,NumEndoSubsetOnSurveilPreDOI_withPrague_Rule3/10,Rule1Prop*10,Rule2Prop*10,Rule3Prop*10))


graph <-
  create_graph(
    nodes_df = nodes,
    edges_df = edges,
    graph_attrs <-c("layout = visNetwork","overlap = FALSE","outputorder = edgesfirst"),
    edge_attrs = "color = white")

# View the graph
render_graph(graph,output = "visNetwork")

#PostDOI
#1. All endoscopies- dataset is EndoSubsetOffSurveilPreDOI for all the endoscopies
NumEndoSubsetOffSurveilPostDOI<-nrow(EndoSubsetOffSurveilPostDOI)
#2. Number lost as Prague -ve
EndoSubsetOffSurveilPostDOI_withPrague<-subset(EndoSubsetOffSurveilPostDOI,!(is.na(EndoSubsetOffSurveilPostDOI$MStage)))
NumEndoSubsetOffSurveilPostDOI_withPrague<-nrow(EndoSubsetOffSurveilPostDOI_withPrague)
#3. Classifiable as Rule
EndoSubsetOffSurveilPostDOI_withPragueNoRules<-subset(EndoSubsetOffSurveilPostDOI,(EndoSubsetOffSurveilPostDOI$FU_Group=="NoRules"))
NumEndoSubsetOffSurveilPostDOI_withPragueNoRules<-nrow(EndoSubsetOffSurveilPostDOI_withPragueNoRules)
#4. Number with Rule 1
EndoSubsetOffSurveilPostDOI_withPrague_Rule1<-subset(EndoSubsetOffSurveilPostDOI,(EndoSubsetOffSurveilPostDOI$FU_Group=="Rule1"))
NumEndoSubsetOffSurveilPostDOI_withPrague_Rule1<-nrow(EndoSubsetOffSurveilPostDOI_withPrague_Rule1)
#5. Number with Rule 2
EndoSubsetOffSurveilPostDOI_withPrague_Rule2<-subset(EndoSubsetOffSurveilPostDOI,(EndoSubsetOffSurveilPostDOI$FU_Group=="Rule2"))
NumEndoSubsetOffSurveilPostDOI_withPrague_Rule2<-nrow(EndoSubsetOffSurveilPostDOI_withPrague_Rule2)
#6. Number with Rule 3
EndoSubsetOffSurveilPostDOI_withPrague_Rule3<-subset(EndoSubsetOffSurveilPostDOI,(EndoSubsetOffSurveilPostDOI$FU_Group=="Rule3"))
NumEndoSubsetOffSurveilPostDOI_withPrague_Rule3<-nrow(EndoSubsetOffSurveilPostDOI_withPrague_Rule3)

#7. Percentage followed up- CAN'T HAVE THIS AS DONT HAVE THE FOLLOW-UP PERIOD
#8. Reason for no follow-up


library(magrittr)
library(DiagrammeR)

AllIndexLabel<-paste("All Index Post DOI",NumEndoSubsetOffSurveilPostDOI, sep = ": ")
nodes <- create_nodes(nodes = c(AllIndexLabel, "Index \nWith Prague", "Rule1","Rule2","Rule3","Followed-up"),
                      label = TRUE,
                      fontsize = 55,
                      fontcolour = "White",
                      type = "lower",
                      style = "filled",
                      color = "aqua",
                      shape = c("circle"),
                      x = c(0,0,-80,0,80,80,0),
                      y = c(280,-80,-250,-250,-250,-400))

edges <- create_edges(from = c(AllIndexLabel, "Index \nWith Prague", "Index \nWith Prague","Index \nWith Prague","Rule3"), 
                      to = c("Index \nWith Prague", "Rule1","Rule2","Rule3", "Followed-up"),
                      rel = c(NumEndoSubsetOffSurveilPostDOI_withPrague, NumEndoSubsetOffSurveilPostDOI_withPrague_Rule1, NumEndoSubsetOffSurveilPostDOI_withPrague_Rule2, 
                              NumEndoSubsetOffSurveilPostDOI_withPrague_Rule3, "To be determined next audit"),
                      arrowhead = rep("normal", 6),
                      # color = c("red", "red", "red", "red", "red", "red"),
                      length = c(200,200,50,50,200,200),
                      fontsize = 55,
                      width=c(NumEndoSubsetOnSurveilPostDOI_withPrague/10,NumEndoSubsetOnSurveilPostDOI_withPrague_Rule1/10,NumEndoSubsetOnSurveilPostDOI_withPrague_Rule2/10,NumEndoSubsetOnSurveilPostDOI_withPrague_Rule3/10,10,100))


graph <-
  create_graph(
    nodes_df = nodes,
    edges_df = edges,
    graph_attrs <-c("layout = visNetwork","overlap = FALSE","outputorder = edgesfirst"),
    edge_attrs = "color = white")

# View the graph
render_graph(graph,output = "visNetwork")

#-----------------------------------------------------------------------------------------------------------------------
#  ----------------------------------#QUESTION- ARE WE SURVEYING PEOPLE PROPERLY?-----------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#  ---------------------------------------------------------------------------------------------------------------------


#Is Surveillance documentation done properly?
#PreDOI
PragueSubsetOnSurveilPreDOI<-subset(EndoSubsetOnSurveilPreDOI,!is.na(EndoSubsetOnSurveilPreDOI$MStage))
PragueSubsetOffSurveilPreDOI<-subset(EndoSubsetOffSurveilPreDOI,!is.na(EndoSubsetOffSurveilPreDOI$MStage))

#Presence of islands        
IslandSubsetOnSurveilPreDOI <- EndoSubsetOnSurveilPreDOI[grepl("[Ii]sland", EndoSubsetOnSurveilPreDOI$ERFINDINGSSTR), ]
IslandSubsetOffSurveilPreDOI<-EndoSubsetOffSurveilPreDOI[grepl("[Ii]sland", EndoSubsetOffSurveilPreDOI$ERFINDINGSSTR), ]

#Hiatus hernia (top of gastric folds)
HerniaSubsetOnSurveilPreDOI <- EndoSubsetOnSurveilPreDOI[grep("[Hh]iat|astric fold|[Pp]inch", EndoSubsetOnSurveilPreDOI$ERFINDINGSSTR), ]
HerniaSubsetOffSurveilPreDOI<-EndoSubsetOffSurveilPreDOI[grep("[Hh]iat|astric fold|[Pp]inch", EndoSubsetOffSurveilPreDOI$ERFINDINGSSTR), ]


#Visible lesions- should also describe the absence of visible lesions explicitly
LesionSubsetOnSurveilPreDOI <- EndoSubsetOnSurveilPreDOI[grep("esion|odule|lcer", EndoSubsetOnSurveilPreDOI$ERFINDINGSSTR), ]
LesionSubsetOffSurveilPreDOI<- EndoSubsetOffSurveilPreDOI[grep("esion|odule|lcer", EndoSubsetOffSurveilPreDOI$ERFINDINGSSTR), ]

#Classification of lesions
#On surveillance vs not on surveillance
#Thie one is done as part of the Therapeutic survey so a different dataset.

#Biopsies (location and samples taken)
#Decided not to do this as no point as all biopsies are labelled at time of pathology so don't see why they should be on the form.
#On surveillance vs not on surveillance


#Pre DOI
n = c(as.numeric(nrow(PragueSubsetOnSurveilPreDOI)/nrow(EndoSubsetOnSurveilPreDOI)), as.numeric(nrow(PragueSubsetOffSurveilPreDOI)/nrow(EndoSubsetOffSurveilPreDOI)), 
      as.numeric(nrow(IslandSubsetOnSurveilPreDOI)/nrow(EndoSubsetOnSurveilPreDOI)),as.numeric(nrow(IslandSubsetOffSurveilPreDOI)/nrow(EndoSubsetOffSurveilPreDOI)),as.numeric(nrow(HerniaSubsetOnSurveilPreDOI)/nrow(EndoSubsetOnSurveilPreDOI)),
      as.numeric(nrow(HerniaSubsetOffSurveilPreDOI)/nrow(EndoSubsetOffSurveilPreDOI)),as.numeric(nrow(LesionSubsetOnSurveilPreDOI)/nrow(EndoSubsetOnSurveilPreDOI)),
      as.numeric(nrow(LesionSubsetOffSurveilPreDOI))/nrow(EndoSubsetOffSurveilPreDOI) )
s = c("On", "Off", "On", "Off","On", "Off","On", "Off") 
b=c("Prague","Prague","Island","Island","Hernia","Hernia","Lesion","Lesion")
EndoMinDataSet<-data.frame(s,b,n)

barchart(b~n,data=EndoMinDataSet,groups=s, 
         scales=list(x=list(rot=0,cex=1.8),y=list(cex=1.8),main=list("Endoscopy documentation for Barrett's PRE DOI",cex=2.8)),   
         key=list(space="right",
                  lines=list(col=c("purple","lightgreen"), lty=c(2,2), lwd=16),
                  text=list(c("On surveillance","Off surveillance"),cex=1.8)
         ))

#PostDOI
#Is Surevillance documentation done properly?

PragueSubsetOnSurveilPostDOI<-subset(EndoSubsetOnSurveilPostDOI,!is.na(EndoSubsetOnSurveilPostDOI$MStage))
PragueSubsetOffSurveilPostDOI<-subset(EndoSubsetOffSurveilPostDOI,!is.na(EndoSubsetOffSurveilPostDOI$MStage))


#Presence of islands        
IslandSubsetOnSurveilPostDOI <- EndoSubsetOnSurveilPostDOI[grepl("[Ii]sland", EndoSubsetOnSurveilPostDOI$ERFINDINGSSTR), ]
IslandSubsetOffSurveilPostDOI<-EndoSubsetOffSurveilPostDOI[grepl("[Ii]sland", EndoSubsetOffSurveilPostDOI$ERFINDINGSSTR), ]

#Hiatus hernia (top of gastric folds)
HerniaSubsetOnSurveilPostDOI <- EndoSubsetOnSurveilPostDOI[grep("[Hh]iat|astric fold|[Pp]inch", EndoSubsetOnSurveilPostDOI$ERFINDINGSSTR), ]
HerniaSubsetOffSurveilPostDOI<-EndoSubsetOffSurveilPostDOI[grep("[Hh]iat|astric fold|[Pp]inch", EndoSubsetOffSurveilPostDOI$ERFINDINGSSTR), ]

#Visible lesions- should also describe the absence of visible lesions explicitly
LesionSubsetOnSurveilPostDOI <- EndoSubsetOnSurveilPostDOI[grep("esion|odule|lcer", EndoSubsetOnSurveilPostDOI$ERFINDINGSSTR), ]
LesionSubsetOffSurveilPostDOI<- EndoSubsetOffSurveilPostDOI[grep("esion|odule|lcer", EndoSubsetOffSurveilPostDOI$ERFINDINGSSTR), ]

#Classification of lesions
#On surveillance vs not on surveillance
#Thie one is done as part of the Therapeutic survey so a different dataset.

#Biopsies (location and samples taken)
#Decided not to do this as no point as all biopsies are labelled at time of pathology so don't see why they should be on the form.
#On surveillance vs not on surveillance


n = c(as.numeric(nrow(PragueSubsetOnSurveilPostDOI)/nrow(EndoSubsetOnSurveilPostDOI)), as.numeric(nrow(PragueSubsetOffSurveilPostDOI)/nrow(EndoSubsetOffSurveilPostDOI)), 
      as.numeric(nrow(IslandSubsetOnSurveilPostDOI)/nrow(EndoSubsetOnSurveilPostDOI)),as.numeric(nrow(IslandSubsetOffSurveilPostDOI)/nrow(EndoSubsetOffSurveilPostDOI)),as.numeric(nrow(HerniaSubsetOnSurveilPostDOI)/nrow(EndoSubsetOnSurveilPostDOI)),
      as.numeric(nrow(HerniaSubsetOffSurveilPostDOI)/nrow(EndoSubsetOffSurveilPostDOI)),as.numeric(nrow(LesionSubsetOnSurveilPostDOI)/nrow(EndoSubsetOnSurveilPostDOI)),
      as.numeric(nrow(LesionSubsetOffSurveilPostDOI))/nrow(EndoSubsetOffSurveilPostDOI) )
s = c("On", "Off", "On", "Off","On", "Off","On", "Off") 
b=c("Prague","Prague","Island","Island","Hernia","Hernia","Lesion","Lesion")
EndoMinDataSet<-data.frame(s,b,n)

barchart(b~n,data=EndoMinDataSet,groups=s, 
         scales=list(x=list(rot=0,cex=1.8),y=list(cex=1.8)),main = ("Endoscopy documentation for Barrett's POST DOI"),
         cex.main=2.8,
         ylab="Documentation",cex.lab=1.8,
         key=list(space="right",
                  lines=list(col=c("purple","lightgreen"), lty=c(2,2), lwd=16),
                  text=list(c("On surveillance","Off surveillance"),cex=1.8)
         ))






#PreDOI
#1. All endoscopies- dataset is EndoSubsetOffSurveilPreDOI for all the endoscopies
NumEndoSubsetOnSurveilPreDOI<-nrow(EndoSubsetOnSurveilPreDOI)
#2. Number lost as Prague -ve
EndoSubsetOnSurveilPreDOI_withPrague<-subset(EndoSubsetOnSurveilPreDOI,!(is.na(EndoSubsetOnSurveilPreDOI$MStage)))
NumEndoSubsetOnSurveilPreDOI_withPrague<-nrow(EndoSubsetOnSurveilPreDOI_withPrague)
#3. Classifiable as Rule
EndoSubsetOnSurveilPreDOI_withPragueNoRules<-subset(EndoSubsetOnSurveilPreDOI,(EndoSubsetOnSurveilPreDOI$FU_Group=="NoRules"))
NumEndoSubsetOnSurveilPreDOI_withPragueNoRules<-nrow(EndoSubsetOnSurveilPreDOI_withPragueNoRules)
#4. Number with Rule 1
EndoSubsetOnSurveilPreDOI_withPrague_Rule1<-subset(EndoSubsetOnSurveilPreDOI,(EndoSubsetOnSurveilPreDOI$FU_Group=="Rule1"))
NumEndoSubsetOnSurveilPreDOI_withPrague_Rule1<-nrow(EndoSubsetOnSurveilPreDOI_withPrague_Rule1)
#5. Number with Rule 2
EndoSubsetOnSurveilPreDOI_withPrague_Rule2<-subset(EndoSubsetOnSurveilPreDOI,(EndoSubsetOnSurveilPreDOI$FU_Group=="Rule2"))
NumEndoSubsetOnSurveilPreDOI_withPrague_Rule2<-nrow(EndoSubsetOnSurveilPreDOI_withPrague_Rule2)
#6. Number with Rule 3
EndoSubsetOnSurveilPreDOI_withPrague_Rule3<-subset(EndoSubsetOnSurveilPreDOI,(EndoSubsetOnSurveilPreDOI$FU_Group=="Rule3"))
NumEndoSubsetOnSurveilPreDOI_withPrague_Rule3<-nrow(EndoSubsetOnSurveilPreDOI_withPrague_Rule3)
#7. Percentage followed up- Get this from below from the non-Rule3NoRepeats field
#8. Reason for no follow-up

#SEE DIAGRAMMER GRAPH AT THE END OF THE SURVEILLANCE SECTION SO CAN SEE ALL THE FOLLOW_UP FOR PRE-DOI AS WELL (can't have this for post DOI as not enough follow-up time)


#PostDOI
#1. All endoscopies- dataset is EndoSubsetOffSurveilPreDOI for all the endoscopies
NumEndoSubsetOnSurveilPostDOI<-nrow(EndoSubsetOnSurveilPostDOI)
#2. Number lost as Prague -ve
EndoSubsetOnSurveilPostDOI_withPrague<-subset(EndoSubsetOnSurveilPostDOI,!(is.na(EndoSubsetOnSurveilPostDOI$MStage)))
NumEndoSubsetOnSurveilPostDOI_withPrague<-nrow(EndoSubsetOnSurveilPostDOI_withPrague)
#3. Classifiable as Rule
EndoSubsetOnSurveilPostDOI_withPragueNoRules<-subset(EndoSubsetOnSurveilPostDOI,(EndoSubsetOnSurveilPostDOI$FU_Group=="NoRules"))
NumEndoSubsetOnSurveilPostDOI_withPragueNoRules<-nrow(EndoSubsetOnSurveilPostDOI_withPragueNoRules)
#4. Number with Rule 1
EndoSubsetOnSurveilPostDOI_withPrague_Rule1<-subset(EndoSubsetOnSurveilPostDOI,(EndoSubsetOnSurveilPostDOI$FU_Group=="Rule1"))
NumEndoSubsetOnSurveilPostDOI_withPrague_Rule1<-nrow(EndoSubsetOnSurveilPostDOI_withPrague_Rule1)
#5. Number with Rule 2
EndoSubsetOnSurveilPostDOI_withPrague_Rule2<-subset(EndoSubsetOnSurveilPostDOI,(EndoSubsetOnSurveilPostDOI$FU_Group=="Rule2"))
NumEndoSubsetOnSurveilPostDOI_withPrague_Rule2<-nrow(EndoSubsetOnSurveilPostDOI_withPrague_Rule2)
#6. Number with Rule 3
EndoSubsetOnSurveilPostDOI_withPrague_Rule3<-subset(EndoSubsetOnSurveilPostDOI,(EndoSubsetOnSurveilPostDOI$FU_Group=="Rule3"))
NumEndoSubsetOnSurveilPostDOI_withPrague_Rule3<-nrow(EndoSubsetOnSurveilPostDOI_withPrague_Rule3)
#7. Percentage followed up- CANT HAVE THIS AS DONT HAVE THE FOLLOW-UP PERIOD
#8. Reason for no follow-up

library(magrittr)
library(DiagrammeR)

AllSurveillanceLabel<-paste("All Surveillance Post DOI",NumEndoSubsetOnSurveilPostDOI, sep = ": ")
nodes <- create_nodes(nodes = c(AllSurveillanceLabel, "Surveillance \nWith Prague", "Rule1","Rule2","Rule3","Followed-up"),
                      label = TRUE,
                      fontsize = 55,
                      fontcolour = "White",
                      type = "lower",
                      style = "filled",
                      color = "aqua",
                      shape = c("circle"),
                      x = c(0,0,-80,0,80,80,0),
                      y = c(280,-80,-250,-250,-250,-400))

edges <- create_edges(from = c(AllSurveillanceLabel, "Surveillance \nWith Prague", "Surveillance \nWith Prague","Surveillance \nWith Prague","Rule3"), 
                      to = c("Surveillance \nWith Prague", "Rule1","Rule2","Rule3", "Followed-up"),
                      rel = c(NumEndoSubsetOnSurveilPostDOI_withPrague, NumEndoSubsetOnSurveilPostDOI_withPrague_Rule1, NumEndoSubsetOnSurveilPostDOI_withPrague_Rule2, 
                              NumEndoSubsetOnSurveilPostDOI_withPrague_Rule3, "To be determined next audit"),
                      arrowhead = rep("normal", 6),
                      # color = c("red", "red", "red", "red", "red", "red"),
                      length = c(200,200,50,50,200,200),
                      fontsize = 55,
                      width=c(NumEndoSubsetOnSurveilPostDOI_withPrague/10,NumEndoSubsetOnSurveilPostDOI_withPrague_Rule1/10,NumEndoSubsetOnSurveilPostDOI_withPrague_Rule2/10,NumEndoSubsetOnSurveilPostDOI_withPrague_Rule3/10,10,100))


graph <-
  create_graph(
    nodes_df = nodes,
    edges_df = edges,
    graph_attrs <-c("layout = visNetwork","overlap = FALSE","outputorder = edgesfirst"),
    edge_attrs = "color = white")

# View the graph
render_graph(graph,output = "visNetwork")


#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#------------------------#How many people are getting repeat endoscopy according to guidelines in the surveillance population:--------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------



#Do for rule 3 only
EndoSubsetOnSurveilPostDOIHospNum<-subset(EndoSubsetOnSurveilPostDOI,EndoSubsetOnSurveilPostDOI$FU_Group=="Rule3")
EndoSubsetOnSurveilPreDOIHospNum<-subset(EndoSubsetOnSurveilPreDOI,EndoSubsetOnSurveilPreDOI$FU_Group=="Rule3")
EndoSubsetOnSurveilPostDOIHospNum<-data.frame(unique(EndoSubsetOnSurveilPostDOIHospNum$HospNum_Id))
EndoSubsetOnSurveilPreDOIHospNum<-data.frame(unique(EndoSubsetOnSurveilPreDOIHospNum$HospNum_Id))

names(EndoSubsetOnSurveilPostDOIHospNum)<-c("x")
names(EndoSubsetOnSurveilPreDOIHospNum)<-c("x")
PreAndPostSurveillance<-rbind(EndoSubsetOnSurveilPreDOIHospNum,EndoSubsetOnSurveilPostDOIHospNum)
Rule3SurveillanceRepeatEndoscopy<-data.frame(PreAndPostSurveillance[duplicated(PreAndPostSurveillance),])
names(Rule3SurveillanceRepeatEndoscopy)<-c("x")
Rule3PropScoped<-round(nrow(Rule3SurveillanceRepeatEndoscopy)/nrow(EndoSubsetOnSurveilPreDOIHospNum),2)



#Do for rule 1 only
EndoSubsetOnSurveilPostDOIHospNum_R1<-subset(EndoSubsetOnSurveilPostDOI,EndoSubsetOnSurveilPostDOI$FU_Group=="Rule1")
EndoSubsetOnSurveilPreDOIHospNum_R1<-subset(EndoSubsetOnSurveilPreDOI,EndoSubsetOnSurveilPreDOI$FU_Group=="Rule1")
EndoSubsetOnSurveilPostDOIHospNum_R1<-data.frame(unique(EndoSubsetOnSurveilPostDOIHospNum_R1$HospNum_Id))
EndoSubsetOnSurveilPreDOIHospNum_R1<-data.frame(unique(EndoSubsetOnSurveilPreDOIHospNum_R1$HospNum_Id))

names(EndoSubsetOnSurveilPreDOIHospNum_R1)<-c("x")
names(EndoSubsetOnSurveilPostDOIHospNum_R1)<-c("x")
PreAndPostSurveillance_R1<-rbind(EndoSubsetOnSurveilPreDOIHospNum_R1,EndoSubsetOnSurveilPostDOIHospNum_R1)
Rule1SurveillanceRepeatEndoscopy<-data.frame(PreAndPostSurveillance_R1[duplicated(PreAndPostSurveillance_R1),])
names(Rule1SurveillanceRepeatEndoscopy)<-c("x")
Rule1PropScoped<-round(nrow(Rule1SurveillanceRepeatEndoscopy)/nrow(EndoSubsetOnSurveilPreDOIHospNum_R1),2)


#Do for rule 2 only
EndoSubsetOnSurveilPostDOIHospNum_R2<-subset(EndoSubsetOnSurveilPostDOI,EndoSubsetOnSurveilPostDOI$FU_Group=="Rule2")
EndoSubsetOnSurveilPreDOIHospNum_R2<-subset(EndoSubsetOnSurveilPreDOI,EndoSubsetOnSurveilPreDOI$FU_Group=="Rule2")
EndoSubsetOnSurveilPostDOIHospNum_R2<-data.frame(unique(EndoSubsetOnSurveilPostDOIHospNum_R2$HospNum_Id))
EndoSubsetOnSurveilPreDOIHospNum_R2<-data.frame(unique(EndoSubsetOnSurveilPreDOIHospNum_R2$HospNum_Id))

names(EndoSubsetOnSurveilPreDOIHospNum_R2)<-c("x")
names(EndoSubsetOnSurveilPostDOIHospNum_R2)<-c("x")
PreAndPostSurveillance_R2<-rbind(EndoSubsetOnSurveilPreDOIHospNum_R2,EndoSubsetOnSurveilPostDOIHospNum_R2)
Rule2SurveillanceRepeatEndoscopy<-data.frame(PreAndPostSurveillance_R2[duplicated(PreAndPostSurveillance_R2),])
names(Rule2SurveillanceRepeatEndoscopy)<-c("x")
Rule2PropScoped<-round(nrow(Rule2SurveillanceRepeatEndoscopy)/nrow(EndoSubsetOnSurveilPreDOIHospNum_R2),2)

n = c(Rule1PropScoped,Rule2PropScoped,Rule3PropScoped )
s = c("Rule1","Rule2","Rule3") 
EndoRuleScopeDataSet<-data.frame(s,n)

barchart(s~n,data=EndoRuleScopeDataSet, 
         main = ("Follow-up endoscopy for each group in surveillance populations"),
         scales=list(x=list(rot=0,cex=2.8),y=list(cex=2.8),main=list(cex=5.8)),
)

library(magrittr)
library(DiagrammeR)

AllSurveillancePreLabel<-paste("All Surveillance Pre DOI",NumEndoSubsetOnSurveilPreDOI, sep = ": ")
AllSurveillanceWithPraguePreLabel<-paste("Surveillance \nWith Prague",NumEndoSubsetOnSurveilPreDOI_withPrague, sep = ": ")
nodes <- create_nodes(nodes = c(AllSurveillancePreLabel, AllSurveillanceWithPraguePreLabel, "Rule1","Rule2","Rule3",
                                "Rule 1 \nF-up","Rule 2 \nF-up","Rule 3 \nF-up"),
                      label = TRUE,
                      fontsize = 55,
                      type = "lower",
                      style = "filled",
                      color = "aqua",
                      shape = c("circle", "circle",
                                "rectangle", "rectangle"),
                      data = c(30.5, 2.6, 9.4, 2.7),
                      x = c(0,0,-80,0,80,-80,0,80),
                      y = c(280,-80,-250,-250,-250,-400,-400,-400))

edges <- create_edges(from = c(AllSurveillancePreLabel, AllSurveillanceWithPraguePreLabel,AllSurveillanceWithPraguePreLabel,AllSurveillanceWithPraguePreLabel,"Rule1","Rule2","Rule3"), 
                      to = c(AllSurveillanceWithPraguePreLabel, "Rule1","Rule2","Rule3","Rule 1 \nF-up","Rule 2 \nF-up","Rule 3 \nF-up"),
                      rel = c(NumEndoSubsetOnSurveilPreDOI_withPrague, NumEndoSubsetOnSurveilPreDOI_withPrague_Rule1, NumEndoSubsetOnSurveilPreDOI_withPrague_Rule2, NumEndoSubsetOnSurveilPreDOI_withPrague_Rule3, Rule1PropScoped,Rule2PropScoped,
                              Rule3PropScoped),
                      arrowhead = rep("normal", 6),
                      #color = c("red", "red", "red", "red", "red", "red"),
                      length = c(200,200,50,50,200,200),
                      width=c(NumEndoSubsetOnSurveilPreDOI_withPrague/10, NumEndoSubsetOnSurveilPreDOI_withPrague_Rule1/10, NumEndoSubsetOnSurveilPreDOI_withPrague_Rule2/10, NumEndoSubsetOnSurveilPreDOI_withPrague_Rule3/10, Rule1PropScoped*10,Rule2PropScoped*10,
                              Rule3PropScoped*10))


graph <-
  create_graph(
    nodes_df = nodes,
    edges_df = edges,
    graph_attrs <-
      c("layout = dot","overlap = FALSE","outputorder = edgesfirst"),
    node_attrs <-
      c("shape = circle",
        "fixedsize = TRUE",
        "width = 100",
        "penwidth = 1",
        "color = DodgerBlue",
        "style = filled",
        "fillcolor = Aqua",
        "alpha_fillcolor = 0.5",
        "fontname = Helvetica",
        "fontcolor = Black"),
    edge_attrs = "color = white")

# View the graph
render_graph(graph,layout=constant,output="visNetwork")


#Now get a list of the patients who were not surveyed to find out what they were-DO BY PATIENT

(PreAndPostSurveillanced<-data.frame(unique(PreAndPostSurveillance$x)))
xu <- PreAndPostSurveillance[duplicated(PreAndPostSurveillance$x),]
PreAndPostSurveillance<-data.frame(unique(PreAndPostSurveillance$x))
names(PreAndPostSurveillance)<-c("x")

names(Rule3SurveillanceRepeatEndoscopyMissed)<-c("x")

Rule3NoRepeats<-data.frame(EndoSubsetOnSurveilPreDOIHospNum[which(!EndoSubsetOnSurveilPreDOIHospNum$x %in% EndoSubsetOnSurveilPostDOIHospNum$x),])












#?---------------------------?---------------------------?---------------------------?---------------------------?---------------------------
#?---------------------------?---------------------------?---------------------------?---------------------------?---------------------------
#?---------------------------Whats the dysplasia detection rate on surveillance pre and post DOI?---------------------------?---------------------------?---------------------------
#?---------------------------?---------------------------?---------------------------?---------------------------?---------------------------
#  ?---------------------------?---------------------------?---------------------------?---------------------------?---------------------------

PreLGD<-EndoSubsetOnSurveilPreDOI[grep("LGD",EndoSubsetOnSurveilPreDOI$IMorNoIM),]
PreHGD<-EndoSubsetOnSurveilPreDOI[grep("HGD",EndoSubsetOnSurveilPreDOI$IMorNoIM),]
PreOAC<-EndoSubsetOnSurveilPreDOI[grep("SM1|SM2|T1b|T1a",EndoSubsetOnSurveilPreDOI$IMorNoIM),]
PostLGD<-EndoSubsetOnSurveilPostDOI[grep("LGD",EndoSubsetOnSurveilPostDOI$IMorNoIM),]
PostHGD<-EndoSubsetOnSurveilPostDOI[grep("HGD",EndoSubsetOnSurveilPostDOI$IMorNoIM),]
PostOAC<-EndoSubsetOnSurveilPostDOI[grep("SM1|SM2|T1b|T1a",EndoSubsetOnSurveilPostDOI$IMorNoIM),]


n = c(nrow(PreLGD), nrow(PreHGD), nrow(PreOAC),nrow(PostLGD), nrow(PostHGD), nrow(PostOAC))
s = c("Pre", "Pre", "Pre", "Post","Post", "Post") 
b=c("LGD","HGD","OAC","LGD","HGD","OAC")
EndoMinDataSet<-data.frame(s,b,n)

barchart(b~n,data=EndoMinDataSet,groups=s, 
         scales=list(x=list(rot=0,cex=1.8),y=list(cex=1.8)),main = ("Dysplasia detection for Surveillance lists"),
         cex.main=2.8,
         ylab="Documentation",cex.lab=1.8,
         key=list(space="right",
                  lines=list(col=c("purple","lightgreen"), lty=c(2,2), lwd=16),
                  text=list(c("Pre 2013","Post 2013"),cex=1.8)
         ))

#--------------------How does the dysplasia detection rate differ between surveillance list vs index list?------------------------------------------------
OnLGD<-EndoSubsetOnSurveilPreDOI[grep("LGD",EndoSubsetOnSurveilPreDOI$IMorNoIM),]
OnHGD<-EndoSubsetOnSurveilPreDOI[grep("HGD",EndoSubsetOnSurveilPreDOI$IMorNoIM),]
OnOAC<-EndoSubsetOnSurveilPreDOI[grep("SM1|SM2|T1b|T1a",EndoSubsetOnSurveilPreDOI$IMorNoIM),]
OffLGD<-EndoSubsetOffSurveilPostDOI[grep("LGD",EndoSubsetOffSurveilPostDOI$IMorNoIM),]
OffHGD<-EndoSubsetOffSurveilPostDOI[grep("HGD",EndoSubsetOffSurveilPostDOI$IMorNoIM),]
OffOAC<-EndoSubsetOffSurveilPostDOI[grep("SM1|SM2|T1b|T1a",EndoSubsetOffSurveilPostDOI$IMorNoIM),]


n = c(nrow(OnLGD), nrow(OnHGD), nrow(OnOAC),nrow(OffLGD), nrow(OffHGD), nrow(OffOAC))
s = c("On", "On", "On", "Off","Off", "Off") 
b=c("LGD","HGD","OAC","LGD","HGD","OAC")
EndoMinDataSet<-data.frame(s,b,n)

barchart(b~n,data=EndoMinDataSet,groups=s, 
         scales=list(x=list(rot=0,cex=1.8),y=list(cex=1.8)),main = ("Dysplasia detection for Surveillance lists"),
         cex.main=2.8,
         ylab="Documentation",cex.lab=1.8,
         key=list(space="right",
                  lines=list(col=c("purple","lightgreen"), lty=c(2,2), lwd=16),
                  text=list(c("Pre 2013","Post 2013"),cex=1.8)
         ))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#  ----------------------------------------------------------------------------------------------------------------------------------------------------------------
#  ------------------------------------------------Therapeutics section---------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
#For EMRs-

#Tissue results breakdown
EndoSubset$MeasurementLargest<-NULL
EndoSubset$NumberBx<-NULL
EndoSubsetEMR<-unique(subset(EndoSubset,EndoSubset$EVENT=="EMR"))
#EndoSubsetEMR<-EndoSubset

#Number of EMRs over time:
EMRWork <- data.frame(EndoSubsetEMR$VisitDate,EndoSubsetEMR$EVENT)
Tots<-EMRWork %>%
  mutate(year = format(EndoSubsetEMR.VisitDate, "%Y")) %>%
  group_by(year)

EMRWorka<-data.frame(table(Tots$year))
EMRWorka$Var1<-as.character(EMRWorka$Var1)


AllEMRs<-nrow(EndoSubsetEMR)
SM2 <- nrow(EndoSubsetEMR[EndoSubsetEMR$IMorNoIM=="SM2",])
SM1 <-nrow(EndoSubsetEMR[EndoSubsetEMR$IMorNoIM=="SM1",])
T1b <-nrow(EndoSubsetEMR[EndoSubsetEMR$IMorNoIM=="T1b",])
T1a <-nrow(EndoSubsetEMR[EndoSubsetEMR$IMorNoIM=="T1a",])
HGD <-nrow(EndoSubsetEMR[EndoSubsetEMR$IMorNoIM=="HGD",])
LGD <-nrow(EndoSubsetEMR[EndoSubsetEMR$IMorNoIM=="LGD",])
IM <-nrow(EndoSubsetEMR[EndoSubsetEMR$IMorNoIM=="IM",])
NoIM <-nrow(EndoSubsetEMR[EndoSubsetEMR$IMorNoIM=="No_IM",])

n = c(as.numeric(AllEMRs), as.numeric(SM2), as.numeric(SM1), as.numeric(T1b),as.numeric(T1a),as.numeric(HGD),as.numeric(LGD),as.numeric(IM),as.numeric(NoIM))
s = c("AllEMRs", "SM2", "SM1","T1b_Unspecified","T1a","HGD","LGD","IM","No IM") 
EMRResult<-data.frame(s,n)
axis(1, at=mids, labels=EMRResult%s)
barplot(EMRResult$n,names.arg=c("AllEMRs", "SM2", "SM1","T1b_Unspec","T1a","HGD","LGD","IM","No IM"),
        xlab = "Tissue grade", ylab = "Number of EMRs",
        cex.lab = 2.0,cex.axis=2.5,cex.main = 2.5,cex.names=2.5,main = "EMR Tissue pathology results")

#Number of Rescue EMRs over time:
RescueEMRs<-EndoSubset %>% 
  arrange(HospNum_Id, as.Date(EndoSubset$VisitDate, '%d/%m/%y')) %>% 
  group_by(HospNum_Id) %>% 
  filter(first(EVENT == "RFA") & EVENT == "EMR")

#Number of EMRs per person:
EMRFreq<-data.frame(table(EndoSubsetEMR$HospNum_Id))
barplot(table(EMRFreq$Freq),xlab = "Number of EMRs", ylab = "Frequency of patients",
        cex.lab = 2.0,cex.axis=2.5,cex.main = 2.5,cex.names=2.5,main = "Frequency distribution of number pf EMRs \nper Barrett's patient reuiqring ablation"  )


#Upstaging of EMRs over time:
library(dplyr)
Upstage<-EndoSubset %>%
  group_by(HospNum_Id) %>% 
  mutate(ind = EVENT=="nothing" & lead(EVENT)=="EMR") %>% 
  slice(sort(c(which(ind),which(ind)+1)))




#RFA
EndoSubsetRFA<-EndoSubset
EndoSubsetRFA$MeasurementLargest<-NULL
EndoSubsetRFA$NumberBx<-NULL
EndoSubsetRFA<-unique(subset(EndoSubsetRFA,EndoSubsetRFA$EVENT=="RFA"))
AllRFAs<-data.frame(unique(EndoSubsetRFA$HospNum_Id))
nrow(AllRFAs)

#Number of RFA's needed per patient
NumRFAPerPt<-EndoSubsetRFA %>% 
  arrange(HospNum_Id) %>% 
  filter(EVENT == "RFA") 
RFAFreq<-data.frame(table(NumRFAPerPt$HospNum_Id))
barplot(table(RFAFreq$Freq),xlab = "Number of RFAs", ylab = "Frequency of patients",
        cex.lab = 2.0,cex.axis=2.5,cex.main = 2.5,cex.names=2.5,main = "Frequency distribution of number pf RFAs \nper Barrett's patient requiring ablation"  )

#RFA catheter breakdown

HALO90a <- EndoSubsetRFA[grepl("90", EndoSubsetRFA$ERDIAGNOSISSTR,perl=TRUE), ]
HALO90b <- EndoSubsetRFA[grepl("90", EndoSubsetRFA$ERFINDINGSSTR,perl=TRUE), ]
HALO90c<-rbind(HALO90a,HALO90b)

HALO360a <- EndoSubsetRFA[grepl("360", EndoSubsetRFA$ERDIAGNOSISSTR,perl=TRUE), ]
HALO360b <- EndoSubsetRFA[grepl("360", EndoSubsetRFA$ERFINDINGSSTR,perl=TRUE), ]
HALO360c<-rbind(HALO360a,HALO360b)

HALO60a <- EndoSubsetRFA[grepl("HALO60| 60", EndoSubsetRFA$ERDIAGNOSISSTR,perl=TRUE), ]
HALO60b <- EndoSubsetRFA[grepl("HALO60| 60", EndoSubsetRFA$ERFINDINGSSTR,perl=TRUE), ]
HALO60c<-rbind(HALO60a,HALO60b)

HALOTTSa <- EndoSubsetRFA[grepl("TTS|[Cc]hannel", EndoSubsetRFA$ERDIAGNOSISSTR,perl=TRUE), ]
HALOTTSb <- EndoSubsetRFA[grepl("TTS|[Cc]hannel", EndoSubsetRFA$ERFINDINGSSTR,perl=TRUE), ]
HALOTTSc<-rbind(HALOTTSa,HALOTTSb)

HALOAPCa <- EndoSubsetRFA[grepl("APC", EndoSubsetRFA$ERDIAGNOSISSTR,perl=TRUE), ]
HALOAPCb <- EndoSubsetRFA[grepl("APC", EndoSubsetRFA$ERFINDINGSSTR,perl=TRUE), ]
HALOAPCc<-rbind(HALOAPCa,HALOAPCb)

TotalNumRFASessions<-EndoSubsetRFA[EndoSubsetRFA$]
n = c(nrow(HALO90c), nrow(HALO360c), nrow(HALO60c), nrow(HALOTTSc),nrow(HALOAPCc))
s = c("HALO 90", "HALO 360", "HALO 60","HALO TTS","APC") 
EMRResult<-data.frame(s,n)
axis(1, at=mids, labels=EMRResult%s)
barplot(EMRResult$n,names.arg=c("HALO 90", "HALO 360", "HALO 60","HALO TTS","APC"),xlab = "Catheter type", 
        ylab = "Number of RFA's",cex.lab = 2.0,cex.axis=2.5,cex.main = 2.5,cex.names=2.5,main = "RFA Catheter type usage"  )





#Number of patients where RFA is finished so assumed CRIM

library(dplyr)
CRIM<-EndoSubset %>%
  group_by(HospNum_Id) %>% 
  mutate(ind = EVENT=="RFA" & lead(EVENT)=="nothing") %>% 
  slice(sort(c(which(ind),which(ind)+1)))



#Work done:
#Number of EMRs over time:
EndoSubsetWorkEMR<-EndoSubsetEMR<-unique(subset(EndoSubset,EndoSubset$EVENT=="EMR"))
TotsWorkEMR<-EndoSubsetWorkEMR %>%
  mutate(year = format(VisitDate, "%Y")) %>%
  select(EVENT,year)%>%
  group_by(year,EVENT)

EndoSubsetWorkEMR<-data.frame(table(TotsWorkEMR$year))
EndoSubsetWorkEMR$Var1<-as.character(EndoSubsetWorkEMR$Var1)

EndoSubsetWorkRFA<-EndoSubsetRFA<-unique(subset(EndoSubset,EndoSubset$EVENT=="RFA"))
TotsWorkRFA<-EndoSubsetWorkRFA %>%
  mutate(year = format(VisitDate, "%Y")) %>%
  select(EVENT,year)%>%
  group_by(year,EVENT)

EndoSubsetWorkRFA<-data.frame(table(TotsWorkRFA$year))
EndoSubsetWorkRFA$Var1<-as.character(EndoSubsetWorkRFA$Var1)

ggplot() + 
  geom_line(data=EndoSubsetWorkEMR,aes(EndoSubsetWorkEMR$Var1, EndoSubsetWorkEMR$Freq,group = 1),colour="blue",size=8)+
  geom_line(data=EndoSubsetWorkRFA,aes(EndoSubsetWorkRFA$Var1, EndoSubsetWorkRFA$Freq,group = 1),colour="red",size=8)+ 
  labs(title="Number of Barrett's Ablations") +
  #geom_point(aes(mergedGroup4$Rptname, mergedGroup4$PercentageChangeForWholeSV), mergedGroup4) + theme(axis.text.x=element_text(angle=-90)) + 
  xlab("Year") + 
  ylab("Number of Ablations") +
  theme(axis.text.x=element_text(angle=-90)) +
  theme(legend.position="top")

#Paris lesion recognition Matrix
EndoSubsetWorkEMR2<-EndoSubsetEMR<-unique(subset(EndoSubset,EndoSubset$EVENT=="EMR"))
EndoSubsetWorkEMR2$ParisClass<- ifelse(grepl("11a_c|2a_c|[Ii][Ii]a_c",EndoSubsetWorkEMR$ERDIAGNOSISSTR,perl=TRUE)|grepl("11a_c|2a_c|[Ii][Ii]a_c",EndoSubsetWorkEMR$ERFINDINGSSTR,perl=TRUE),"2a_c",
                                       ifelse(grepl("[Ii][Ii]a|2a|11a",EndoSubsetWorkEMR$ERDIAGNOSISSTR,perl=TRUE)|grepl("[Ii][Ii]a|2a|11a",EndoSubsetWorkEMR$ERFINDINGSSTR,perl=TRUE),"2a",
                                              ifelse(grepl("[Ii][Ii]b|2b|11b",EndoSubsetWorkEMR$ERDIAGNOSISSTR,perl=TRUE)|grepl("[Ii][Ii]b|2b|11b",EndoSubsetWorkEMR$ERFINDINGSSTR,perl=TRUE),"2b",
                                                     ifelse(grepl("[Ii][Ii]c|2c|11c",EndoSubsetWorkEMR$ERDIAGNOSISSTR,perl=TRUE)|grepl("[Ii][Ii]c|2c|11c",EndoSubsetWorkEMR$ERFINDINGSSTR,perl=TRUE),"2c",
                                                            ifelse(grepl("[Ii][Ii][Ii]|III",EndoSubsetWorkEMR$ERDIAGNOSISSTR,perl=TRUE)|grepl("[Ii][Ii][Ii]|III",EndoSubsetWorkEMR$ERFINDINGSSTR,perl=TRUE),"3",
                                                                   ifelse(grepl("Paris [Tt]ype [Ii]s|1s ",EndoSubsetWorkEMR$ERDIAGNOSISSTR,perl=TRUE)|grepl("Paris [Tt]ype [Ii]s|1s",EndoSubsetWorkEMR$ERFINDINGSSTR,perl=TRUE),"1s",
                                                                          ifelse(grepl("[Ii]p|1p",EndoSubsetWorkEMR$ERDIAGNOSISSTR,perl=TRUE)|grepl("[Ii]p|1p",EndoSubsetWorkEMR$ERFINDINGSSTR,perl=TRUE),"1p",
                                                                                 "No_Paris")))))))
#Create the matrix
df3 <- data.frame(EndoSubsetWorkEMR2$ParisClass,EndoSubsetWorkEMR2$IMorNoIM)
#Reorganise the column names and rows
#Get rid of no Paris EMR's
dfy<-df3[!df3$EndoSubsetWorkEMR2.ParisClass=="No_Paris",]
#Get the histology proportions by the Paris grade
tr4<-as.data.frame.matrix(prop.table(table(dfy),1) )
tr4<-tr4[c("T1b_Unspec","SM1","SM2","T1a","HGD","LGD","IGD","IM","No_IM")]

tr5<-as.matrix(tr4)
tr5<-head(tr5,-1)
#Create  the heatmap
heatmap.2(tr5,trace="none",dendrogram='none', Rowv=FALSE, Colv=FALSE)




#Heatmap by column 
trb<-as.data.frame.matrix(prop.table(table(df3),2) )
trb<-trb[c("T1b_Unspec","SM1","SM2","T1a","HGD","LGD","IGD","IM","No_IM")]
trb<-as.matrix(trb)

heatmap.2(trb,trace="none",dendrogram='none', Rowv=FALSE, Colv=FALSE)
#Follow-up timings of patients undergoing therapy


#Missed cancer rates


library(XLConnect)
wk = loadWorkbook("/Users/sebastianzeki/Dropbox/Work/Medical/Clinical/Gastro/BarrettsAudit/Upper GI cancer data sebastian.xls")
dfw = readWorksheet(wk, sheet="MergedCancersOfInterest",header=TRUE)
#Rename column with HospitalNumber to allow merge later on with Endoscopy dataset
#I need the cancers that develop for people who have had a surveillance endoscopy before the cancer developed but also any endoscopy before the 
#cancer developed
names(dfw)[1] <- c("HospNum_Id")


#Clean my dataset so no therapeutics just index and surveillance endoscopies
EndoSubsetForMissedCancerDetection <- EndoSubset[!grepl("Therapeutic", EndoSubset$ERPROCEDUREPERFORMED), ]
EndoSubsetForMissedCancerDetection <- EndoSubsetForMissedCancerDetection[!grepl("EMR|ucosal resection", EndoSubsetBx$Diagnosis), ]
EndoSubsetForMissedCancerDetection <- EndoSubsetForMissedCancerDetection[grepl("astroscopy", EndoSubsetBx$ERPROCEDUREPERFORMED), ]

#Get rid of MeasurementLargest and NumberBx
EndoSubsetForMissedCancerDetection$NumberBx<-NULL
EndoSubsetForMissedCancerDetection$MeasurementLargest<-NULL

#Remove duplicates
EndoSubsetForMissedCancerDetection<-unique(EndoSubsetForMissedCancerDetection)

#Now merge the hospital numbers from the CancerRegistry
#Need to do this as a full outer merge so all records merged regardless of whether meatched or not
PotentialMissedCancers <- merge(EndoSubsetForMissedCancerDetection,dfw,by=c("HospNum_Id"))


#The dysplasia detection rate:
#Need to include indefinite for dysplasia
EndoSubsetOffSurveilPreDOIDDR <- EndoSubsetOffSurveilPreDOI[!grepl("Therapeutic", EndoSubsetOffSurveilPreDOI$ERPROCEDUREPERFORMED), ]
EndoSubsetOffSurveilPreDOIDDR <- EndoSubsetOffSurveilPreDOIDDR[!grepl("EMR|ucosal resection", EndoSubsetOffSurveilPreDOIDDR$Diagnosis), ]
EndoSubsetOffSurveilPreDOIDDR <- EndoSubsetOffSurveilPreDOIDDR[grepl("astroscopy", EndoSubsetOffSurveilPreDOIDDR$ERPROCEDUREPERFORMED), ]
EndoSubsetOffSurveilPreDOIDDR_HGD<-nrow(subset(EndoSubsetOffSurveilPreDOIDDR,EndoSubsetOffSurveilPreDOIDDR$IMorNoIM=="HGD"))
EndoSubsetOffSurveilPreDOIDDR_LGD<-nrow(subset(EndoSubsetOffSurveilPreDOIDDR,EndoSubsetOffSurveilPreDOIDDR$IMorNoIM=="LGD"))
EndoSubsetOffSurveilPreDOIDDR_IGD<-nrow(subset(EndoSubsetOffSurveilPreDOIDDR,EndoSubsetOffSurveilPreDOIDDR$IMorNoIM=="IGD"))
EndoSubsetOffSurveilPreDOIDDR_Malig<-nrow(EndoSubsetOffSurveilPreDOIDDR[grepl("SM2|SM1|T1a|T1b_Unspec",EndoSubsetOffSurveilPreDOIDDR$IMorNoIM),])
EndoSubsetOffSurveilPreDOIDDRProp=(EndoSubsetOffSurveilPreDOIDDR_HGD+EndoSubsetOffSurveilPreDOIDDR_IGD+EndoSubsetOffSurveilPreDOIDDR_LGD+EndoSubsetOffSurveilPreDOIDDR_Malig)/nrow(EndoSubsetOffSurveilPreDOIDDR)*100

#PreDOI Off Surveillance
EndoSubsetOnSurveilPreDOIDDR <- EndoSubsetOnSurveilPreDOI[!grepl("Therapeutic", EndoSubsetOnSurveilPreDOI$ERPROCEDUREPERFORMED), ]
EndoSubsetOnSurveilPreDOIDDR <- EndoSubsetOnSurveilPreDOIDDR[!grepl("EMR|ucosal resection", EndoSubsetOnSurveilPreDOIDDR$Diagnosis), ]
EndoSubsetOnSurveilPreDOIDDR <- EndoSubsetOnSurveilPreDOIDDR[grepl("astroscopy", EndoSubsetOnSurveilPreDOIDDR$ERPROCEDUREPERFORMED), ]
EndoSubsetOnSurveilPreDOIDDR_HGD<-nrow(subset(EndoSubsetOnSurveilPreDOIDDR,EndoSubsetOnSurveilPreDOIDDR$IMorNoIM=="HGD"))
EndoSubsetOnSurveilPreDOIDDR_LGD<-nrow(subset(EndoSubsetOffSurveilPreDOIDDR,EndoSubsetOffSurveilPreDOIDDR$IMorNoIM=="LGD"))
EndoSubsetOnSurveilPreDOIDDR_IGD<-nrow(subset(EndoSubsetOnSurveilPreDOIDDR,EndoSubsetOnSurveilPreDOIDDR$IMorNoIM=="IGD"))

EndoSubsetOnSurveilPreDOIDDR_Malig<-nrow(EndoSubsetOnSurveilPreDOIDDR[grepl("SM2|SM1|T1a|T1b_Unspec",EndoSubsetOnSurveilPreDOIDDR$IMorNoIM),])
EndoSubsetOnSurveilPreDOIDDRProp=(EndoSubsetOnSurveilPreDOIDDR_HGD+EndoSubsetOnSurveilPreDOIDDR_IGD+EndoSubsetOnSurveilPreDOIDDR_LGD+EndoSubsetOnSurveilPreDOIDDR_Malig)/nrow(EndoSubsetOnSurveilPreDOIDDR)*100

#PostDOI On Surveillance
EndoSubsetOffSurveilPostDOIDDR <- EndoSubsetOffSurveilPostDOI[!grepl("Therapeutic", EndoSubsetOffSurveilPostDOI$ERPROCEDUREPERFORMED), ])
EndoSubsetOffSurveilPostDOIDDR <- EndoSubsetOffSurveilPostDOIDDR[!grepl("EMR|ucosal resection", EndoSubsetOffSurveilPostDOIDDR$Diagnosis), ])
EndoSubsetOffSurveilPostDOIDDR <- EndoSubsetOffSurveilPostDOIDDR[grepl("astroscopy", EndoSubsetOffSurveilPostDOIDDR$ERPROCEDUREPERFORMED), ]
EndoSubsetOffSurveilPostDOIDDR_HGD<-nrow(subset(EndoSubsetOffSurveilPostDOIDDR,EndoSubsetOffSurveilPostDOIDDR$IMorNoIM=="HGD"))
EndoSubsetOffSurveilPostDOIDDR_LGD<-nrow(subset(EndoSubsetOffSurveilPostDOIDDR,EndoSubsetOffSurveilPostDOIDDR$IMorNoIM=="LGD"))
EndoSubsetOffSurveilPostDOIDDR_IGD<-nrow(subset(EndoSubsetOffSurveilPostDOIDDR,EndoSubsetOffSurveilPostDOIDDR$IMorNoIM=="IGD"))
EndoSubsetOffSurveilPostDOIDDR_Malig<-nrow(EndoSubsetOffSurveilPreDOIDDR[grepl("SM2|SM1|T1a|T1b_Unspec",EndoSubsetOffSurveilPreDOIDDR$IMorNoIM),])
EndoSubsetOffSurveilPostDOIDDRProp=(EndoSubsetOffSurveilPostDOIDDR_HGD+EndoSubsetOffSurveilPostDOIDDR_IGD+EndoSubsetOffSurveilPostDOIDDR_LGD+EndoSubsetOffSurveilPostDOIDDR_Malig)/nrow(EndoSubsetOffSurveilPostDOIDDR)*100


#PostDOI Off Surveillance
EndoSubsetOnSurveilPostDOIDDR <- EndoSubsetOnSurveilPostDOI[!grepl("Therapeutic", EndoSubsetOnSurveilPostDOI$ERPROCEDUREPERFORMED), ]
EndoSubsetOnSurveilPostDOIDDR <- EndoSubsetOnSurveilPostDOIDDR[!grepl("EMR|ucosal resection", EndoSubsetOnSurveilPostDOIDDR$Diagnosis), ]
EndoSubsetOnSurveilPostDOIDDR <- EndoSubsetOnSurveilPostDOIDDR[grepl("astroscopy", EndoSubsetOnSurveilPostDOIDDR$ERPROCEDUREPERFORMED), ]
EndoSubsetOnSurveilPostDOIDDR_HGD<-nrow(subset(EndoSubsetOnSurveilPostDOIDDR,EndoSubsetOnSurveilPostDOIDDR$IMorNoIM=="HGD"))
EndoSubsetOnSurveilPostDOIDDR_LGD<-nrow(subset(EndoSubsetOnSurveilPostDOIDDR,EndoSubsetOnSurveilPostDOIDDR$IMorNoIM=="LGD"))
EndoSubsetOnSurveilPostDOIDDR_IGD<-nrow(subset(EndoSubsetOnSurveilPostDOIDDR,EndoSubsetOnSurveilPostDOIDDR$IMorNoIM=="IGD"))
EndoSubsetOnSurveilPostDOIDDR_Malig<-nrow(EndoSubsetOnSurveilPreDOIDDR[grepl("SM2|SM1|T1a|T1b_Unspec",EndoSubsetOnSurveilPreDOIDDR$IMorNoIM),])
EndoSubsetOnSurveilPostDOIDDRProp=(EndoSubsetOnSurveilPostDOIDDR_HGD+EndoSubsetOnSurveilPostDOIDDR_IGD+EndoSubsetOnSurveilPostDOIDDR_LGD+EndoSubsetOnSurveilPostDOIDDR_Malig)/nrow(EndoSubsetOnSurveilPostDOIDDR)*100


#Present as a heatmap
#Create a matrix with the data

B = matrix( 
  +   c(EndoSubsetOffSurveilPreDOIDDR_Malig, EndoSubsetOnSurveilPreDOIDDR_Malig, EndoSubsetOffSurveilPostDOIDDR_Malig, EndoSubsetOnSurveilPostDOIDDR_Malig, EndoSubsetOffSurveilPreDOIDDR_HGD, EndoSubsetOnSurveilPreDOIDDR_HGD,EndoSubsetOffSurveilPostDOIDDR_HGD,EndoSubsetOnSurveilPostDOIDDR_HGD,EndoSubsetOffSurveilPreDOIDDR_LGD,EndoSubsetOnSurveilPreDOIDDR_LGD,
        EndoSubsetOffSurveilPostDOIDDR_LGD,EndoSubsetOnSurveilPostDOIDDR_LGD,EndoSubsetOffSurveilPreDOIDDR_IGD,EndoSubsetOnSurveilPreDOIDDR_IGD,EndoSubsetOffSurveilPostDOIDDR_IGD,EndoSubsetOnSurveilPostDOIDDR_IGD),nrow=4,ncol=4)


rownames(B)  <- c("Pre_OffSurv","Pre_OnSurv","Post_OffSurv","Post_OnSurv")
colnames(B)  <- c("Malig","HGD","LGD","IGD")
heatmap.2(B,trace="none",dendrogram='none', Rowv=TRUE, Colv=TRUE)


C = matrix( 
  +   c(EndoSubsetOffSurveilPreDOIDDR_Malig/nrow(EndoSubsetOnSurveilPostDOIDDR)*100, 
        EndoSubsetOnSurveilPreDOIDDR_Malig/nrow(EndoSubsetOnSurveilPostDOIDDR)*100, 
        EndoSubsetOffSurveilPostDOIDDR_Malig/nrow(EndoSubsetOnSurveilPostDOIDDR)*100, 
        EndoSubsetOnSurveilPostDOIDDR_Malig/nrow(EndoSubsetOnSurveilPostDOIDDR)*100, 
        EndoSubsetOffSurveilPreDOIDDR_HGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100, 
        EndoSubsetOnSurveilPreDOIDDR_HGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100,
        EndoSubsetOffSurveilPostDOIDDR_HGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100,
        EndoSubsetOnSurveilPostDOIDDR_HGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100,
        EndoSubsetOffSurveilPreDOIDDR_IGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100, 
        EndoSubsetOnSurveilPreDOIDDR_IGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100,
        EndoSubsetOffSurveilPostDOIDDR_IGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100,
        EndoSubsetOnSurveilPostDOIDDR_IGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100,
        EndoSubsetOffSurveilPreDOIDDR_LGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100,
        EndoSubsetOnSurveilPreDOIDDR_LGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100,
        EndoSubsetOffSurveilPostDOIDDR_LGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100,
        EndoSubsetOnSurveilPostDOIDDR_LGD/nrow(EndoSubsetOnSurveilPostDOIDDR)*100),nrow=4,ncol=4)
rownames(C)  <- c("Pre_OffSurv","Pre_OnSurv","Post_OffSurv","Post_OnSurv")
colnames(C)  <- c("Malig","HGD","IGD","LGD")

#Overall dysplasia detection rate:
n = c(EndoSubsetOffSurveilPreDOIDDRProp, EndoSubsetOnSurveilPreDOIDDRProp, EndoSubsetOffSurveilPostDOIDDRProp, EndoSubsetOnSurveilPostDOIDDRProp)
s = c("Pre_Off", "Pre_On","Post_Off","Post_On") 
DDR<-data.frame(s,n)

barplot(t(C),names.arg=c("Pre_Off", "Pre_On","Post_Off","Post_On"),
        xlab = "OGD type and date", ylab = "Dysplasia detection proportion (%)",
        cex.lab = 2.0,cex.axis=2.5,cex.main = 2.5,cex.names=2.5,main = "Dysplasia detection proportion \n2010-2013 vs 2013-2016 by OGD type")
legend("topleft", 
       legend = c("IGD","LGD","HGD","Adenocarcinoma"), 
       fill = c("lightgray","lightgray","darkgray","black"))


#Sankey diagram creation

#Get the procedure done- probably just do this for people who have had procedures
#Group by patient




library(dplyr)
EndoSubset<-unique(EndoSubset)
EndoSubset <- EndoSubset[order(EndoSubset$HospNum_Id) , ]

EndoSubsetOffSurveilPreDOI[!(EndoSubsetOffSurveilPreDOI$HospNum_Id %in% EnrolledAfterIndexEndoscopy$x),]

TherapySet<-EndoSubset
TherapySet$EVENT<-gsub("nothing","",TherapySet$EVENT)
#TherapySet<-TherapySet[TherapySet$EVENT!="nothing",]
Sankey<-dcast(setDT(TherapySet)[, if(any(EVENT=="EMR"|EVENT=="RFA")) .SD, HospNum_Id], HospNum_Id~rowid(HospNum_Id), value.var ="EVENT")
library(data.table)
PtFlow<-Sankey
PtFlow<-data.frame(PtFlow)
names(PtFlow)<-c("ord1","ord2","ord3","ord4","ord5","ord6","ord7","ord8","ord9","ord10","ord11","ord12")

orders <- PtFlow %>%
  select(ord1, ord2, ord3, ord4, ord5,ord6,ord7,ord8,ord9,ord10,ord11,ord12)

orders.plot <- data.frame()
library(googleVis)
for (i in 3:ncol(orders)) {
  
  ord.cache <- orders %>%
    group_by(orders[ , i-1], orders[ , i]) %>%
    summarise(n=n())
  
  colnames(ord.cache)[1:2] <- c('from', 'to')
  
  # adding tags to carts
  ord.cache$from <- paste(ord.cache$from, '(', i-1, ')', sep='')
  ord.cache$to <- paste(ord.cache$to, '(', i, ')', sep='')
  
  orders.plot <- rbind(orders.plot, ord.cache) 
  
}

orders.plot<-data.frame(orders.plot)
orders.plot<-orders.plot[grepl("[A-Z]",orders.plot$from)&grepl("[A-Z]",orders.plot$to), ]
orders.plot<-orders.plot[!grepl("NA",orders.plot$from)&!grepl("NA",orders.plot$to), ]
plot(gvisSankey(orders.plot, from='from', to='to', weight='n',
                options=list(height=900, width=1800, 
                             sankey="{link:{color:{fill:'black',stroke: 'black', strokeWidth: 1 }},
node: { color: { fill: '#a61d4c' },
                                     label: { color: '#871b47',fontName: 'Open Sans',fontSize: 35 } }}")))


#Then transpose converting empty cells to NA
#Then label each source
#Then provide the targets and the value which is the summarised value for each source type


