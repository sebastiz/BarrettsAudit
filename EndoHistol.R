#Personalised lesion recognition script
library(XLConnect)
library(stringr)
#Endoscopy query spreadsheet- should be filtered by all endoscopies done by a certain person
#Histology spreadsheet should be filtered by organ eg either by oesophagus of by stomach or both etc

#Import the spreadsheets

Histol = loadWorkbook("/Users/sebastianzeki/Documents/PhysJava/Dumper_EndoscOrHistol/BarrettsSurveillanceReport.xls")
Histol = readWorksheet(Histol, sheet="Data",header=TRUE)

#Firstly split the histology spreadsheet up so that the histology report chopped up into Clinical Details etc in new columns

Histol$ClinDetails<- str_extract_all(Histol$ResultText,regex("(CLINICAL DETAILS.*?)MACROSCOPICAL DESCRIPTION", dotall = TRUE,multiline = TRUE))
Histol$MacroDesc<-str_extract_all(Histol$ResultText,regex("(MACROSCOPICAL DESCRIPTION.*?)HISTOLOGY", dotall = TRUE,multiline = TRUE))
Histol$HistolPath<-str_extract_all(Histol$ResultText,regex("(HISTOLOGY.*?)DIAGNOSIS", dotall = TRUE,multiline = TRUE))
Histol$Dx<-str_extract_all(Histol$ResultText,regex("(DIAGNOSIS.*)", dotall = TRUE,multiline = TRUE))

#Rename the column to merge on with endoscopy
names(Histol)[names(Histol) == 'ResultPerformed'] <- 'VisitDate'
names(Histol)[names(Histol) == 'PatientID'] <- 'HospNum_Id'

Histol$VisitDate<-as.Date(Histol$VisitDate , format="%y/%m/%d")


Endoscop = loadWorkbook("/Users/sebastianzeki/Documents/PhysJava/Dumper_EndoscOrHistol/BarrettEndoscopyatients.xls")
Endoscop =  readWorksheet(Endoscop, sheet="Sheet1",header=TRUE)
names(Endoscop)[names(Endoscop) == 'PATIENTRECORDNUMBER'] <- 'HospNum_Id'
names(Endoscop)[names(Endoscop) == 'DATEOFPROCEDURE'] <- 'VisitDate'Endoscop$VisitDate<-as.Date(Endoscop$VisitDate , format="%d/%m/%Y")


#merge the Endoscopy and Histology spreadsheets by Hospital Number and VisitDate
EndoscHistol<-merge(Endoscop,Histol,by=c("HospNum_Id","VisitDate"))

#---------------------#-----------------------------------------------#-------------------------------------------------------------------------
#---------------------#-----------------------------------------------#-------------------------------------------------------------------------
#---------------------#-----------------------------------------------#-------------------------------------------------------------------------
#---------------------#-----------------------------------------------#-------------------------------------------------------------------------
#---------------------#-----------------------------------------------#-------------------------------------------------------------------------
#---------------------#-----------------------------------------------#-------------------------------------------------------------------------
#---------------------#-----------------------------------------------#-------------------------------------------------------------------------
#---------------------#-----------------------------------------------#-------------------------------------------------------------------------
#---------------------#-----------------------------------------------#-------------------------------------------------------------------------
#---------------------#-----------------------------------------------#-------------------------------------------------------------------------
#---------------------#-----------------------------------------------#-------------------------------------------------------------------------
#---------------------#-----------------------------------------------#-------------------------------------------------------------------------

#Then do the analysis
EndoscHistolQuickView<-data.frame(EndoscHistol["ER_FINDINGS_STR"],EndoscHistol["ER_DIAGNOSIS_STR"],EndoscHistol["HospNum_Id"],EndoscHistol["VisitDate"],EndoscHistol["Dx"])

df<-data.frame(table(Endoscop$INDICATIONS))

df<-subset(df,df$Freq>10)


barplot(df$Freq,las=3,
cex.lab = 1.0,cex.axis=2.5,cex.main = 2.5,cex.names=1.5,names.arg=df$Var1,
space=1,main = "Reasons for Endoscopy",angle = 45)

#Endoscopic parameters audit of self

