#' ---
#' title: "Analysis of norms of science Version 2 and political views and IQ "
#' ---
rm(list=ls())
library(lattice)
library(stargazer)
library(tidyverse)
library(dplyr)
library(lme4)
library(lavaan)
library(Hmisc)
library(psych)
library(haven) #to read SPSS
library(semPlot)
library(corrplot)


#center and average indicators for a construct ------------------------------------------------------------
formcomp <- function(indicators) {
  centind <- as.data.frame(scale(indicators))
  return(apply(centind,1,mean))
}


##--- Begin by reading data and variable names -----------------------------------------------------------------------
#This does not work from Sweave or Knitr: setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#Define working directory without changing it for weaving from other folder with paper in it.
wd <- "C:/Users/Lewan/Documents/Research Projects/Climate Change/religion human special status and science Klaus/Study 2 -- for revision/qualtrics details and data study 2"

source(paste(wd,"religionQ2funcs.R",sep="/"))
#reading SPSS file for a change
newNorms <- read_sav(paste(wd,"Full+sample+IQ+and+norms_March+25,+2019_13.44.sav",sep="/")) 

#read variable names for table of raw responses
vn <- read.csv(paste(wd,"varnamesnorms2Q.csv",sep="/"),header=TRUE,row.names=NULL,stringsAsFactors = FALSE) 
names(vn) <- c("shortvn","longvn")

#AFQs are satisfied:
table(newNorms$AFQ1)
table(newNorms$AFQ2)
newNorms1.5 <- newNorms  #potential subsetting step with elimination of variables and so on
summary(newNorms1.5 %>% select(NOR_COM1:SDO8_REV))
x11()
hist(newNorms1.5 %>% select(starts_with("IQ")))


##--- identify keyhitters before reverse-scoring ----------------------------------------------------------------------
neutral <- 0 #if set to zero, any sequence of identical keys is eliminated. If set to 4, only non-neutral responses are dropped
keyhitters <- NULL
keyhitters <- newNorms1.5 %>% select(starts_with("NOR_")) %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0))
keyhitters <- cbind(keyhitters, newNorms1.5 %>% select(starts_with("POL_")) %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0)))
keyhitters <- cbind(keyhitters, newNorms1.5 %>% select(starts_with("IQ_"))  %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0)))
keyhitters <- cbind(keyhitters, newNorms1.5 %>% select(starts_with("CLIM"))      %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0)))
keyhitters <- cbind(keyhitters, newNorms1.5 %>% select(starts_with("VAX"))     %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0)))
keyhitters <- cbind(keyhitters, newNorms1.5 %>% select(starts_with("SDO"))   %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0)))

table(rowSums(keyhitters))
newNorms1.5 <-  newNorms1.5  %>% filter(!(rowSums(keyhitters)>1))

males2 <- table(newNorms1.5$Q67)["1"]
females2 <- table(newNorms1.5$Q67)["2"]
nonbin2 <- table(newNorms1.5$Q67)["3"]
mage2 <- round(mean(newNorms1.5$Age_1),1)
mdage2 <- round(median(newNorms1.5$Age_1),1)
minage2 <- min(newNorms1.5$Age_1)
maxage2 <- max(newNorms1.5$Age_1)

#**----------- get raw responses before reverse scoring
#function to grab a row and interleave with parentheses for printing
interleave <- function(x){
  hlx <- length(x)/2
  retstr <- paste(sapply(1:hlx, FUN=function(i) paste(" & ",as.character(x[i])," & (",as.character(x[i+hlx]),")",sep="")),collapse="",sep="")
}
itemResppercent <- newNorms1.5 %>% select(NOR_COM1:SDO8_REV) %>% lapply(table) %>% lapply(as.numeric) %>%
  lapply(FUN=function(x) c(x,round(x/sum(x)*100))) 


#this generates latex code for insertion into document by first massaging variable names for paper
vnlookup <- read.table(paste(wd,"vartranslation.txt",sep="/"),stringsAsFactors = FALSE, sep=",")
t4l <- NULL
for (i in 1:length(itemResppercent)) {
  massagedvarnames <- vnlookup[vnlookup[,1]==names(itemResppercent[i]),2]
  massagedvarnames <- str_replace_all(massagedvarnames,"_","\\\\_")
  mychar <-paste(massagedvarnames, interleave(itemResppercent[[i]]), "\\","\\", sep="")
  t4l<-rbind(t4l,mychar,deparse.level = 0)
}
writeLines(sort(t4l[  substr(names(itemResppercent),1,2)=="NO"  ,]),con="_t.norms2.tex") #sort by subscale
writeLines(t4l[  substr(names(itemResppercent),1,2)=="PO"  ,][-5],con="_t.conservatism2.tex")  #drop slider item with different scale
writeLines(t4l[  substr(names(itemResppercent),1,2)=="IQ"  ,],con="_t.malleability2.tex")
writeLines(t4l[  substr(names(itemResppercent),1,2)=="CL"  ,], con="_t.climate2.tex")
writeLines(t4l[  substr(names(itemResppercent),1,2)=="VA"  ,],con="_t.vax2.tex")
writeLines(t4l[  substr(names(itemResppercent),1,2)=="SD"  ,],con="_t.sdo2.tex")


#* --------------------- reverse score and summarize ---------------------------
newNorms2 <- newNorms1.5 %>% mutate_at(vars(contains("_REV")),revscore,mm=7) 
# compute pairwise correlations within each cluster
for (cluster in c("POL","NOR_COM","NOR_DIS","NOR_UNI","IQ","CLIM","VAX","SDO")) {
  newNorms2 %>% select(starts_with(cluster)) %>% cor(use="complete.obs") %>% print
}


# plot distributions of all composites before potentially dropping indicator variables
pdf(file="histoSummary2.pdf",height=9,width=6.5) #this will go into paper directory that is being weaved (since no setwd)
par(mfrow=c(3,2))
newNorms2 %>% select(starts_with("NOR")) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="Norms of science",col="light gray")
newNorms2 %>% select(starts_with("POL")) %>% select(-contains("SLI")) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="Conservatism",col="light gray")
newNorms2 %>% select(starts_with("IQ")) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="IQ largely heritable",col="light gray")
newNorms2 %>% select(starts_with("CLI")) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="Climate science",col="light gray")
newNorms2 %>% select(starts_with("VAX")) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="Vaccinations",col="light gray")
newNorms2 %>% select(starts_with("SDO")) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="Social Dominance",col="light gray")
dev.off()



#*--------- compute measurement models for all constructs so we can use single-indicator models later --------------

# Polarity: Norms all pointing towards greater endorsement of norms of science
NormvarsAll <- c(names(newNorms2 %>% select(starts_with("NOR_UNI"))))
fmmod <- c("NORuniversal =~ ",paste(NormvarsAll,collapse=" + "))
invisible(fmgof<-fitMM (fmmod,newNorms2)) 

NormvarsAll <- c(names(newNorms2 %>% select(starts_with("NOR_DIS"),-NOR_DIS2)))
fmmod <- c("NORdisinterest =~ ",paste(NormvarsAll,collapse=" + "))
invisible(fmgof<-fitMM (fmmod,newNorms2)) 

NormvarsAll <- c(names(newNorms2 %>% select(starts_with("NOR_COM"))))
fmmod <- c("NORcommun =~ ",paste(NormvarsAll,collapse=" + "),"NOR_COM1 ~~     NOR_COM2")
invisible(fmgof<-fitMM (fmmod,newNorms2)) 

#* ---- now compute composites for sub-norms that can yield measurement model for common Norms factor
composites <- NULL
composites$POLconscomp <- newNorms2 %>% select(starts_with("POL")) %>% apply(1,mean,na.rm=TRUE)
composites$NORunicomp <- newNorms2 %>% select(starts_with("NOR_UNI")) %>% apply(1,mean)
composites$NORdiscomp <- newNorms2 %>% select(starts_with("NOR_DIS")) %>% apply(1,mean)
composites$NORcomcomp <- newNorms2 %>% select(starts_with("NOR_COM")) %>% apply(1,mean)
composites$NORscepcomp <- newNorms2 %>% select(starts_with("NOR_SCEP")) %>% apply(1,mean)
composites$AllNorms <- newNorms2 %>% select(starts_with("NOR_")) %>% apply(1,mean)
composites <- as.data.frame(composites)

normcorrs <- cor(composites)
cor.mtest(composites)$p < .1

#* ---- print table with component correlations involving norms
t4l <- NULL
normcorvarnames<- c("Conservatism","Universalism","Disinterestedness","Communism","Scepticism","All norms")
for (i in 2:dim(normcorrs)[1]) {
  rowofcors <- paste(sapply(1:i-1, FUN=function(j) paste(sprintf("%5.3f",normcorrs[i,j])," & ",sep="")),collapse="",sep="")
  mychar <-paste(normcorvarnames[i], rowofcors, "\\","\\", sep="")
  t4l<-rbind(t4l,mychar,deparse.level = 0)
}
writeLines(t4l,con="_t.normcors.tex")



NORvars <- c(names(composites %>% select(starts_with("NOR"),-AllNorms)))
fmmod <- c("NOR =~ ",paste(NORvars,collapse=" + "))
invisible(norgof<-fitMM (fmmod, composites)) 

#* ---- now look at other constructs -----------------------------------------------

# Polarity: pointing towards increased conservatism
POLvars <- c(names(newNorms2 %>% select(starts_with("POL"))))
fmmod <- c("POLcons =~ ",paste(POLvars,collapse=" + "),"POL_CONS2 ~~      POL_CONS4_REV")
invisible(polgof<-fitMM (fmmod,newNorms2)) 

# Polarity: pointing towards greater heritability
IQvars <- c(names(newNorms2 %>% select(starts_with("IQ"),-IQ_GEN5_REV))) 
fmmod <- c("IQ =~ ",paste(IQvars,collapse=" + "),"IQ_GEN2_REV ~~ IQ_GEN4_REV")
invisible(IQgof2<-fitMM (fmmod,newNorms2)) 

climvars <- c(names(newNorms2 %>% select(starts_with("CLIM"))))
fmmod <- c("clim =~ ",paste(climvars,collapse=" + "),"CLIM1_REV ~~ CLIM5_REV")
invisible(climgof2<-fitMM (fmmod,newNorms2)) 

# CLIM1_REV ==  CNatFluct 
# CLIM2     ==  CdueGHG 
# CLIM3     ==  CseriousDamage 
# CLIM4     ==  CO2causesCC  
# CLIM5_REV ==  HumansInsign  


vaxvars <- c(names(newNorms2 %>% select(starts_with("VAX"))))
fmmod <- c("vax =~ ",paste(vaxvars,collapse=" + "),"VAX2_REV ~~ VAX4_REV")
invisible(vaxgof2<-fitMM (fmmod,newNorms2)) 

# VAX1	    ==    VaxSafe  
# VAX2_REV	==  	VaxNegSide  
# VAX3	    ==  	VaxTested  
# VAX4_REV	==  	VaxRisky 
# VAX5	    ==  	VaxContribHealth  


#* ---- SDO stuff from C:/Users/Lewan/Documents/MikTek/bibtex/papers/Womick/data for Womick18 and code
#items explained in Ho15
#The items used here are items 1-8 (using same numbering) as in the SDO_7(s) scale, Appendix B, Ho15
# Polarity: pointing towards greater social dominance
sdo7s <- NULL
sdo7s$SDOcomp1 <- apply(cbind(newNorms2$SDO1,newNorms2$SDO4_REV),1,mean)
sdo7s$SDOcomp2 <- apply(cbind(newNorms2$SDO2,newNorms2$SDO3_REV),1,mean)
sdo7s$SDOcomp3 <- apply(cbind(newNorms2$SDO5,newNorms2$SDO7_REV),1,mean)
sdo7s$SDOcomp4 <- apply(cbind(newNorms2$SDO6,newNorms2$SDO8_REV),1,mean)
sdo7s <- as.data.frame(sdo7s)

cor(sdo7s)

SDOvars <- c(names(sdo7s %>% select(starts_with("SDOcomp"))))
fmmod <- c("SDO =~ ",paste(SDOvars,collapse=" + "),"SDOcomp1 ~~ SDOcomp2")
invisible(sdogof<-fitMM (fmmod,sdo7s)) 

#* ---- combine all variables into new data set for main analysis
newNorms3 <- cbind(newNorms2,composites,sdo7s) %>% 
    select(c(contains("comp"),contains("POL"),contains("IQ"),contains("CLIM"),contains("VAX"))) %>%
    select(-IQ_GEN5_REV) %>% select(-POLconscomp)

#* ---- compute single-indicators models
(IQ.SI   <- singleindmodel(IQvars,list(c("IQ_GEN2_REV", "IQ_GEN4_REV")),newNorms3))
(clim.SI <- singleindmodel(climvars, list(c("CLIM1_REV","CLIM5_REV")),newNorms3))
(vax.SI  <- singleindmodel(vaxvars,list(c("VAX2_REV","VAX4_REV")),newNorms3))
(SDO.SI  <- singleindmodel(SDOvars,list(c("SDOcomp1","SDOcomp2")),newNorms3))
(NOR.SI  <- singleindmodel(NORvars,NULL,newNorms3))
(POL.SI  <- singleindmodel(POLvars,list(c("POL_CONS2","POL_CONS4_REV")),newNorms3))

#put error estimates for single indicators into named array for use in function
eSImods <- c(IQ.SI$eSImod, clim.SI$eSImod, vax.SI$eSImod, SDO.SI$eSImod,  NOR.SI$eSImod, POL.SI$eSImod)
names(eSImods) <- c ("IQ","clim","vax","SDO","NOR","POL")

#compute composite scores for the SI models
compositeNewNorms <- data.frame ( 
  IQ   = apply(newNorms3[,IQvars], 1,mean),
  clim = apply(newNorms3[,climvars], 1,mean),
  vax  = apply(newNorms3[,vaxvars], 1,mean),
  SDO  = apply(newNorms3[,SDOvars], 1,mean),
  NOR  = apply(newNorms3[,NORvars], 1,mean),
  POL  = apply(newNorms3[,POLvars], 1,mean), 
               newNorms3)

#nice plots for correlations
RM  <- cor(select(compositeNewNorms,IQ:POL))
RM2 <- cor.mtest(select(compositeNewNorms,IQ:POL), conf.level = .95)
diag(RM)<-NA

x11(width=11,height=10)
colnames(RM) <- rownames(RM) <- c("IQ herit","Climate","Vax","SDO","Norms","Conservatism")
corrplot.mixed(RM, lower.col = "black", number.cex = 1.2,insig = "blank",
               upper="ellipse",p.mat = RM2$p,sig.level=.05,
               tl.pos="lt",tl.col="black",na.label = ".",bg="lightgray",
               tl.cex=1.5,tl.srt=60)
#this goes into "current directory" which will be the tex directory when script is called from .Rnw file
dev.print(pdf,"cormat.pdf")




#* ----  correlation structure among 6 unidimensional latent constructs
modelCorrel <- c("
                 IQfac   =~ IQ
                 climFac =~ clim
                 vaxFac  =~ vax
                 SDOFac  =~ SDO
                 NORFac  =~ NOR
                 POLFac  =~ POL
                 
                  IQ ~~ ",   IQ.SI$eSImod,   "*IQ",
                 "clim ~~",  clim.SI$eSImod, "*clim", 
                 "vax ~~ ",  vax.SI$eSImod,  "*vax",
                 "SDO ~~ ",  SDO.SI$eSImod,  "*SDO", 
                 "NOR ~~ ",  NOR.SI$eSImod,  "*NOR",
                 "POL ~~ ",  POL.SI$eSImod,  "*POL" )
fitCorrel <- sem(modelCorrel, compositeNewNorms, std.lv=TRUE, estimator="ML")
summary(fitCorrel,standardized=TRUE, fit.measures=TRUE)     
x11()
semPaths(fitCorrel, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE,
         structural=TRUE, layout="circle",rotation=1,edge.label.cex=1.5)

#get correlation matrix for latent variables
lvcormat <- lavInspect(fitCorrel, what = "cor.lv")
colnames(lvcormat) <- rownames(lvcormat) <- (c("IQ heritable", "Climate", "Vaccinations", "SDO", "Norms of science", "Conservatism"))
cormat <- stargazer(lvcormat, title="Correlations among 6 unidimensional latent variables") 
write.table(cormat[10:18],file="_t.lvcor2.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)


corrplot.mixed(lvcormat, lower.col = "black", number.cex = 1.2,insig = "blank",
               upper="ellipse",
               tl.pos="lt",tl.col="black",na.label = ".",bg="lightgray",
               tl.cex=1.5,tl.srt=60)


#**--------------- single factor models ----------------------------------------------
modClim <- c("climFac  ~ SDOFac + POLFac + NORFac   

                  climFac =~ clim
                  SDOFac  =~ SDO
                  NORFac  =~ NOR
                  POLFac  =~ POL

                  clim ~~",  clim.SI$eSImod, "*clim", 
                 "SDO ~~ ",  SDO.SI$eSImod,  "*SDO", 
                 "NOR ~~ ",  NOR.SI$eSImod,  "*NOR",
                 "POL ~~ ",  POL.SI$eSImod,  "*POL")

fitClim  <- sem(modClim, compositeNewNorms, std.lv=TRUE, estimator="ML")
summary(fitClim,standardized=TRUE, fit.measures=TRUE)
smclimgof2 <- fitMeasures(fitClim)
#parameterEstimates(fitall2ffull, standardized=TRUE)
modindices(fitClim, sort.=TRUE, maximum.number=4)
semPaths(fitClim, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE,
         structural=TRUE, layout="tree",rotation=1,edge.label.cex=1.5)

modVax <- c("vaxFac  ~ SDOFac +  NORFac   

             vaxFac =~ vax
             SDOFac  =~ SDO
             NORFac  =~ NOR
             POLFac  =~ POL
             
             vax ~~ ",  vax.SI$eSImod,  "*vax",
             "SDO ~~ ",  SDO.SI$eSImod,  "*SDO", 
             "NOR ~~ ",  NOR.SI$eSImod,  "*NOR",
             "POL ~~ ",  POL.SI$eSImod,  "*POL")

fitVax  <- sem(modVax, compositeNewNorms, std.lv=TRUE, estimator="ML")
summary(fitVax,standardized=TRUE, fit.measures=TRUE)
smvaxgof2 <- fitMeasures(fitVax)
#parameterEstimates(fitall2ffull, standardized=TRUE)
modindices(fitVax, sort.=TRUE, maximum.number=4)
semPaths(fitVax, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE,
         structural=TRUE, layout="tree",rotation=1,edge.label.cex=1.5)

#mediate climate
mediateVax <- c("NORFac  ~ alpha2 * SDOFac 

                climFac   ~ beta1 * NORFac +  direct3 * POLFac

                  climFac  =~ clim
                  SDOFac  =~ SDO
                  NORFac  =~ NOR
                POLFac  =~ POL

                 clim ~~ ",  vax.SI$eSImod,  "*clim",
                "SDO ~~ ",  SDO.SI$eSImod,  "*SDO",
                "NOR ~~ ",  NOR.SI$eSImod,  "*NOR",
                "POL ~~ ",  POL.SI$eSImod,  "*POL")
mediateVaxfit  <- sem(mediateVax, compositeNewNorms, std.lv=TRUE, estimator="ML")
summary(mediateVaxfit,standardized=TRUE, fit.measures=TRUE)
modindices(mediateVaxfit, sort.=TRUE, maximum.number=4)

modIQ <- c("IQFac  ~ SDOFac + NORFac   

            IQFac =~ IQ
            SDOFac  =~ SDO
            NORFac  =~ NOR
            POLFac  =~ POL
            
            IQ ~~ ",  IQ.SI$eSImod,  "*IQ",
            "SDO ~~ ",  SDO.SI$eSImod,  "*SDO", 
            "NOR ~~ ",  NOR.SI$eSImod,  "*NOR",
            "POL ~~ ",  POL.SI$eSImod,  "*POL")

fitIQ  <- sem(modIQ, compositeNewNorms, std.lv=TRUE, estimator="ML")
summary(fitIQ,standardized=TRUE, fit.measures=TRUE)
smIQgof2 <- fitMeasures(fitIQ)
#parameterEstimates(fitall2ffull, standardized=TRUE)
modindices(fitIQ, sort.=TRUE, maximum.number=4)
semPaths(fitIQ, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE,
         structural=TRUE, layout="tree",rotation=1,edge.label.cex=1.5) 


#**--------------- Now fit all scientific constructs together
allfull2fmod <- c("
                  climFac  ~ SDOFac + POLFac + NORFac
                  vaxFac ~ NORFac
                  IQfac ~  SDOFac  
                  IQfac ~~ 0*vaxFac
                  
                  IQfac   =~ IQ
                  climFac =~ clim
                  vaxFac  =~ vax
                  SDOFac  =~ SDO
                  NORFac  =~ NOR
                  POLFac  =~ POL
                  
                  IQ ~~ ",   IQ.SI$eSImod,   "*IQ",
                  "clim ~~",  clim.SI$eSImod, "*clim", 
                  "vax ~~ ",  vax.SI$eSImod,  "*vax",
                  "SDO ~~ ",  SDO.SI$eSImod,  "*SDO", 
                  "NOR ~~ ",  NOR.SI$eSImod,  "*NOR",
                  "POL ~~ ",  POL.SI$eSImod,  "*POL")

fitall2ffull  <- sem(allfull2fmod, compositeNewNorms, std.lv=TRUE, estimator="ML")
summary(fitall2ffull,standardized=TRUE, fit.measures=TRUE)
fitallgof2 <- fitMeasures(fitall2ffull)
#parameterEstimates(fitall2ffull, standardized=TRUE)
modindices(fitall2ffull, sort.=TRUE, maximum.number=4)
semPaths(fitall2ffull, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE,
                   structural=TRUE, layout="tree",rotation=1,edge.label.cex=1.5) 


#experiment with mediation models
# indirectclim := alpha2 * beta1
# indirectvax := alpha2 * beta2
# totalclim:= indirectclim + direct1
# proportionclim := (indirectclim)/totalclim
mediateAll <- c("NORFac  ~ alpha2 * SDOFac 

                climFac  ~  beta1 * NORFac  + direct22 * POLFac 
                  
                vaxFac   ~ beta2 * NORFac 
                IQfac    ~ direct2 * SDOFac 
                

                IQfac ~~ 0*vaxFac
                IQfac ~~ 0*climFac
                
                  vaxFac  =~ vax
                  IQfac   =~ IQ
                  climFac =~ clim
                  SDOFac  =~ SDO
                  NORFac  =~ NOR
                  POLFac  =~ POL
                  
                  clim ~~",  clim.SI$eSImod, "*clim", 
                  "IQ ~~ ",   IQ.SI$eSImod,   "*IQ", 
                  "vax ~~ ",  vax.SI$eSImod,  "*vax",
                  "SDO ~~ ",  SDO.SI$eSImod,  "*SDO", 
                  "NOR ~~ ",  NOR.SI$eSImod,  "*NOR",
                  "POL ~~ ",  POL.SI$eSImod,  "*POL")
mediateAllfit  <- sem(mediateAll, compositeNewNorms, std.lv=TRUE, estimator="ML")
summary(mediateAllfit,standardized=TRUE, fit.measures=TRUE)
modindices(mediateAllfit, sort.=TRUE, maximum.number=4)
x11()
 semPaths(mediateAllfit, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE,
          structural=TRUE, layout="tree",rotation=1,edge.label.cex=1.5) 

 
mediateAll2 <- c("NORFac  ~ alpha2 * SDOFac 

                 climFac  ~ direct1 * POLFac + beta1 * NORFac 
                 vaxFac   ~  beta2 * NORFac  +  direct11 * SDOFac 
                 
                 indirectclim := alpha2 * beta1
                 indirectvax := alpha2 * beta2
                 totalclim:= indirectclim + direct1
                 proportionclim := (indirectclim)/totalclim
                 
                 vaxFac  =~ vax
  
                 climFac =~ clim
                 SDOFac  =~ SDO
                 NORFac  =~ NOR
                 POLFac  =~ POL
                 
                 clim ~~",  clim.SI$eSImod, "*clim", 
              
                 "vax ~~ ",  vax.SI$eSImod,  "*vax",
                 "SDO ~~ ",  SDO.SI$eSImod,  "*SDO", 
                 "NOR ~~ ",  NOR.SI$eSImod,  "*NOR",
                 "POL ~~ ",  POL.SI$eSImod,  "*POL")
 mediateAllfit2  <- sem(mediateAll2, compositeNewNorms, std.lv=TRUE, estimator="ML")
 summary(mediateAllfit2,standardized=TRUE, fit.measures=TRUE)
 modindices(mediateAllfit2, sort.=TRUE, maximum.number=4)
 x11()
 semPaths(mediateAllfit2, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE,
          structural=TRUE, layout="tree",rotation=1,edge.label.cex=1.5) 
 
