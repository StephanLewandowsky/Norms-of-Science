## ---- qchunk
#analyze qualtrics religion representative sample 
rm(list=ls())

library(lattice)
library(stargazer)
library(tidyverse)
library(lme4)
library(lavaan)
library(semPlot)


#**-------------- define working directory but don't change to it so all output files end up in folder where paper is being compiled
wd <- "C:/Users/lewan/Documents/Research Projects/Climate Change/religion human special status and science Klaus/Q rep sample/Q religion data"
source(paste(wd,"religionQfuncs.R",sep="/"))

#This does not work from Sweave or Knitr: setwd(dirname(rstudioapi::getSourceEditorContext()$path))
relig <- read.csv(paste(wd,"Attitudes_towards_science_for_Qualtrics_Representative_Sample (2).csv",sep="/"),header=TRUE,row.names=NULL) 
#row 2 (with verbose names) manually deleted from Q file. 
# 10.3.17: Also deleted last record with single stray number in one column (after soft launch)

#read variable names for table of raw responses
vn <- read.csv(paste(wd,"varNames.csv",sep="/"),header=TRUE,row.names=NULL,stringsAsFactors = FALSE) 
vnfin <- vn[grep("Q3",vn$qvarname)[-length(grep("Q3",vn$qvarname))],]

#**----------- clean up data
relig15 <- relig %>% filter(Q2.1 == Q4.2) %>%
  filter(Q3.40 == 4) %>%              #table is not an animal
  select(contains("Q")) %>% 
  select(-contains("Q_Tot"))          #drop the Q that ain't a q.

# first fix the Qualtrics-induced scale problems
relig15 <- relig15 %>% mutate_at(c(paste("Q3.",c(2:19),sep=""),
                                   paste("Q3.",c(25:39),sep=""),"Q4.1"),fixscore,mm=14) %>%
  mutate_at("Q2.4",fixscore,mm=15) %>% mutate_at("Q3.1",fixscore,mm=22) %>%
  mutate_at("Q3.20",fixscore,mm=28)

# identify people who hit the same key (excepting neutral) for all items in a cluster
neutral <- 0 #if set to zero, any sequence of identical keys is eliminated. If set to 4, only non-neutral responses are dropped
keyhitters <- cbind(
  (relig15 %>% select(num_range("Q3.",1:5)) %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0))),  #humans special
  (relig15 %>% select(num_range("Q3.",6:12)) %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0))),  #nationalism is good 
  (relig15 %>% select(num_range("Q3.",13:19)) %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0))), #IQ environmental
  (relig15 %>% select(num_range("Q3.",20:24)) %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0))), #religiosity
  (relig15 %>% select(num_range("Q3.",25:29)) %>%apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0))), #free markets
  (relig15 %>% select(num_range("Q3.",30:34)) %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0))), #climate
  (relig15 %>% select(num_range("Q3.",35:39)) %>% apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0)))) #vax
#eliminate the key hitters
table(rowSums(keyhitters))
relig15 <- filter(relig15,!(rowSums(keyhitters)>1))

#demographics
males <- table(relig15$Q2.2)["1"]
females <- table(relig15$Q2.2)["2"]
mage <- round(mean(relig15$Q2.1),1)
mdage <- round(median(relig15$Q2.1),1)
minage <- min(relig15$Q2.1)
maxage <- max(relig15$Q2.1)
lateEnglish <- table(relig15$Q2.3)["5"]

#**----------- get raw responses before reverse scoring
#function to grab a row and interleave with parentheses for printing
interleave <- function(x){
  hlx <- length(x)/2
  retstr <- paste(sapply(1:hlx, FUN=function(i) paste(" & ",as.character(x[i])," & (",as.character(x[i+hlx]),")",sep="")),collapse="",sep="")
}
itemResppercent <- relig15 %>% select(num_range("Q3.",1:39)) %>% lapply(table) %>% lapply(as.numeric) %>%
  lapply(FUN=function(x) c(x,round(x/sum(x)*100))) 

#this generates latex code for insertion into document
t4l <- NULL
for (i in 1:length(itemResppercent)) {
  mychar <-paste(vnfin$shortname[i], interleave(itemResppercent[[i]]), "\\","\\", sep="")
  t4l<-rbind(t4l,mychar,deparse.level = 0)
}
write.table(t4l[1:5,],file="_t.exceptionalism.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[6:12,],file="_t.nationalism.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[13:19,],file="_t.malleability.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[20:24,],file="_t.religiosity.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[25:29,],file="_t.freemarket.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[30:34,],file="_t.climate.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[35:39,],file="_t.vax.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)

#now reverse score such that polarity is:
# humans are special 
# nationalism is great
# IQ is environmentally determined and malleable
# being religious
# free market endorsement
# accept climate change
# accept vaccinations
relig2 <- relig15 %>% mutate_at(c("Q3.4",
                                  "Q3.6", "Q3.7", "Q3.10", 
                                  "Q3.14", "Q3.15", "Q3.16", "Q3.19",  #polarity of IQ towards heritability
                                  "Q3.26", "Q3.28", "Q3.29",
                                  "Q3.30", "Q3.34",
                                  "Q3.36", "Q3.38"),revscore,mm=7) %>% 
  mutate_at(c("Q3.21","Q3.23","Q3.24"),revscore,mm=5)

# compute pairwise correlations within each cluster
relig2 %>% select(num_range("Q3.",1:5)) %>% cor(use="complete.obs")   #humans special
relig2 %>% select(num_range("Q3.",6:12)) %>% cor(use="complete.obs")  #nationalism is good 
relig2 %>% select(num_range("Q3.",13:19)) %>% cor(use="complete.obs") #IQ environmental
relig2 %>% select(num_range("Q3.",20:24)) %>% cor(use="complete.obs") #religiosity
relig2 %>% select(num_range("Q3.",25:29)) %>% cor(use="complete.obs") #free markets
relig2 %>% select(num_range("Q3.",30:34)) %>% cor(use="complete.obs") #climate
relig2 %>% select(num_range("Q3.",35:39)) %>% cor(use="complete.obs") #vax


# construct histogram for summary statistics
pdf(file="histoSummary.pdf",height=9,width=6.5) #this will go into paper directory that is being weaved (since no setwd)
par(mfrow=c(4,2))
relig2 %>% select(num_range("Q3.",1:5)) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="Exceptionalism",col="light gray")   #humans special
relig2 %>% select(num_range("Q3.",6:12)) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="Nationalism",col="light gray")
relig2 %>% select(num_range("Q3.",13:19)) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="IQ largely heritable",col="light gray")  #IQ environmental
relig2 %>% select(num_range("Q3.",20:24)) %>% rowMeans %>% hist(las=1,xlim=c(1,6),xlab="Average score",main="Religiosity",col="light gray")  #religiosity
relig2 %>% select(num_range("Q3.",25:29)) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="Free market",col="light gray") #free markets
relig2 %>% select(num_range("Q3.",30:34)) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="Climate science",col="light gray")  #climate
relig2 %>% select(num_range("Q3.",35:39)) %>% rowMeans %>% hist(las=1,xlim=c(1,7),xlab="Average score",main="Vaccinations",col="light gray")  #vax
dev.off()

#write.csv(relig2,paste(wd,"QreligProcessed.csv",sep="/"),row.names=FALSE)

#**-------------- now compute measurement models for all constructs so we can use single-indicator models later
#exceptionalism
humvars <- paste("Q3.",c(1:5),sep="")
humspecmod <- c("humspec =~ ",paste(humvars,collapse=" + "),"Q3.3 ~~ Q3.5")
humspecgof <- fitMM(humspecmod,relig2)

#nationalism
intvars <- paste("Q3.",c(7,8,9,10,11,12),sep="") 
internatmod <- c("nationalism   =~ ", paste(intvars,collapse=" + "),"Q3.7 ~~ Q3.10")
internatgof <- fitMM (internatmod,relig2)

# separate IQ models for positive and reverse-scored items
IQ.env <- paste("Q3.",c(19,14:16),sep="")
IQmodel.env <-c("IQ.env =~ ",paste(IQ.env,collapse=" + "),NULL)
fitMM (IQmodel.env,relig2)

IQ.her <- paste("Q3.",c(13,17,18),sep="")
IQmodel.her <-c("IQ.her =~ ",paste(IQ.her,collapse=" + "), NULL)
fitMM (IQmodel.her,relig2)

# two-correlated factors model
IQmodel.twocf <- c("IQ.her =~ ", paste(IQ.her, collapse=" + "), "\n", 
                   "IQ.env =~ ", paste(IQ.env, collapse=" + "))
twocfacIQfit <- sem(IQmodel.twocf,relig2)
twocfacIQcor <- lavInspect(twocfacIQfit, what = "cor.lv")
twocfacIQp   <- pnorm(abs(inspect(twocfacIQfit,what="est")$psi/inspect(twocfacIQfit,what="se")$psi),lower.tail = FALSE)*2

# hierarchical two-factor model
IQmodel.two <- c("IQ.her =~ ", paste(IQ.her, collapse=" + "), "\n", 
                 "IQ.env =~ ", paste(IQ.env, collapse=" + "), "\n",
                 "IQ.two =~ a*IQ.her + a*IQ.env")
twofacIQfit <- sem(IQmodel.two,relig2)
summary(twofacIQfit, standardized=TRUE, fit.measures=TRUE)
twofacIQgof <- fitmeasures(twofacIQfit)
semPaths(twofacIQfit, "std", title =FALSE, curvePivot = TRUE)

# bi-factor IQ model for positive and reverse-scored items
IQ.bi <- paste("Q3.", c(13:19), sep="") 
IQ.her <- paste("Q3.",c(13,17,18),sep="")
IQmodel.bi <- c("IQ.bi =~ ", paste(IQ.bi, collapse=" + "), "\n", 
                "IQ.her =~ ",  paste(IQ.her, collapse= " + "), "\n", 
                "IQ.bi ~~ 0*IQ.her")
IQbigof <- fitMM (IQmodel.bi,relig2)

# only consider items that are tentative in their propositions
IQ.tentative <- paste("Q3.", c(14,15,17,19), sep="") 
IQmodel.tent <- c("IQ.tent =~ ", paste(IQ.tentative, collapse=" + "))
IQtentgof <- fitMM (IQmodel.tent,relig2)

#religiosity
Relvars <- paste("Q3.",c(20:24),sep="")
religmodel <- c("religiosity =~ ",paste(Relvars,collapse=" + "))
religgof <- fitMM (religmodel,relig2)

#FM 
FMvars <- paste("Q3.",c(25:29),sep="")
FMmodel <- c("FM =~ ", paste(FMvars,collapse=" + "),"Q3.25 ~~ Q3.27")   #as before for PLOS ONE
fmgof <- fitMM (FMmodel,relig2)

# hierarchical two-factor model for conservatism
Consmodel.two <- c("C.FM =~ ", paste(FMvars,  collapse=" + "),"Q3.25 ~~ Q3.27", "\n", 
                   "C.Rel =~ ", paste(Relvars, collapse=" + "), "\n",
                   "C.Int =~ ", paste(intvars,collapse=" + "),"Q3.7 ~~ Q3.10","\n",
                   "C.hum =~ ", paste(humvars,collapse=" + "),"Q3.3 ~~ Q3.5","\n",
                   "C.two =~ NA*C.FM + C.Rel + C.Int + C.hum","\n",
                   "C.two ~~ 1*C.two") 
C2fit <- sem(Consmodel.two,relig2)
summary(C2fit, standardized=TRUE, fit.measures=TRUE)


#climate 
Climvars <- paste("Q3.",c(30:34),sep="")
climmodel <- c("climate =~ ", paste(Climvars,collapse=" + "), "Q3.30 ~~ Q3.34")   #as before for PLOS ONE
climategof <- fitMM (climmodel,relig2)

#vaccination 
Vaxvars <- paste("Q3.",c(35:39),sep="")
vaxmodel <- c("vax =~ ",paste(Vaxvars,collapse=" + "), "Q3.36 ~~ Q3.38")   #as before for PLOS ONE
vaxgof <- fitMM (vaxmodel,relig2)


#**----------- compute single-indicators models
(humSI  <- singleindmodel(humvars,list(c("Q3.3","Q3.5")),relig2))
(intSI  <- singleindmodel(intvars, list(c("Q3.7","Q3.10")),relig2))
(RelSI  <- singleindmodel(Relvars,NULL,relig2))
(FMSI   <- singleindmodel(FMvars,list(c("Q3.25","Q3.27")),relig2))
(ClimSI <- singleindmodel(Climvars,list(c("Q3.30","Q3.34")),relig2))
(VaxSI  <- singleindmodel(Vaxvars,list(c("Q3.36","Q3.38")),relig2))
(IQSI   <- singleindmodel(IQ.tentative,NULL,relig2))
#put error estimates for single indicators into named array for use in function
eSImods <- c(humSI$eSImod, intSI$eSImod, RelSI$eSImod, FMSI$eSImod,  ClimSI$eSImod,  VaxSI$eSImod)
names(eSImods) <- c ("hum","int","Rel","FM","Clim","Vax")


#compute composite scores for the SI models
compositeRelig <- data.frame ( 
  hum =  apply(relig2[,humvars], 1,mean),
  int =  apply(relig2[,intvars], 1,mean),
  Rel =  apply(relig2[,Relvars], 1,mean),
  FM =   apply(relig2[,FMvars], 1,mean),
  Clim = apply(relig2[,Climvars], 1,mean),
  Vax =  apply(relig2[,Vaxvars], 1,mean),
  IQ  =  apply(relig2[,IQ.tentative], 1,mean),
  select(relig2,num_range("Q3.",13:29)))

#**------ correlation structure among 6 unidimensional latent constructs
modelCorrel <- c("
                 humFac =~ hum
                 intFac =~ int
                 RelFac =~ Rel
                 FMFac =~ FM
                 ClimFac =~ Clim
                 VaxFac =~ Vax
                 
                 hum ~~ ", humSI$eSImod, "*hum",
                 "int ~~",  intSI$eSImod, "*int", 
                 "Rel ~~ ", RelSI$eSImod, "*Rel",
                 "FM ~~ ",  FMSI$eSImod,  "*FM", 
                 "Clim ~~ ",ClimSI$eSImod,"*Clim",
                 "Vax ~~ ", VaxSI$eSImod, "*Vax" 
)
fitCorrel <- sem(modelCorrel, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fitCorrel,standardized=TRUE, fit.measures=TRUE)

#**------ correlation structure among all 7 unidimensional latent constructs (including tentative IQ items)
modelxCorrel <- c("
               humFac =~ hum
               intFac =~ int
               RelFac =~ Rel
               FMFac =~ FM
               ClimFac =~ Clim
               VaxFac =~ Vax
                IQFac =~ IQ

                  hum ~~ ", humSI$eSImod, "*hum",
                  "IQ ~~", IQSI$eSImod, "*IQ",
                  "int ~~",  intSI$eSImod, "*int", 
                  "Rel ~~ ", RelSI$eSImod, "*Rel",
                  "FM ~~ ",  FMSI$eSImod,  "*FM", 
                  "Clim ~~ ",ClimSI$eSImod,"*Clim",
                  "Vax ~~ ", VaxSI$eSImod, "*Vax" 
)
fitCorrelx <- sem(modelxCorrel, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fitCorrelx,standardized=TRUE, fit.measures=TRUE)
#parameterEstimates(fitCorrelx, standardized=TRUE)

#get correlation matrix for latent variables
lvcormat <- lavInspect(fitCorrel, what = "cor.lv")
colnames(lvcormat) <- rownames(lvcormat) <- (c("Exceptionalism", "Nationalism", "Religiosity", "Free market", "Climate", "Vaccinations"))
cormat <- stargazer(lvcormat, title="Correlations among 6 unidimensional latent variables") 
write.table(cormat[10:18],file="_t.lvcor.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)

#now compute p-values
pvals2tailed <- pnorm(abs(inspect(fitCorrel,what="est")$psi/inspect(fitCorrel,what="se")$psi),lower.tail = FALSE)*2
colnames(pvals2tailed) <- rownames(pvals2tailed) <- (c("Exceptionalism", "Nationalism", "Religiosity", "Free market", "Climate", "Vaccinations"))
pvals2tailed[upper.tri(pvals2tailed)] <- 0
maxloc <- which(pvals2tailed == max(pvals2tailed), arr.ind = TRUE) 
maxpval <- max(pvals2tailed)
rownames(pvals2tailed)[maxloc[1]]
colnames(pvals2tailed)[maxloc[2]]


#**----------- correlations including the two-factor hierarchical model for IQ
modelCorrel3 <- c("
               humFac =~ hum
               intFac =~ int
               RelFac =~ Rel
               FMFac =~ FM
               ClimFac =~ Clim
               VaxFac =~ Vax
               IQ.her =~ Q3.13 + Q3.17 + Q3.18 
               IQ.env =~ Q3.19 + Q3.14 + Q3.15 + Q3.16
               IQ.two =~ a*IQ.her + a*IQ.env

               hum ~~ ", humSI$eSImod, "*hum",
                  "int ~~",  intSI$eSImod, "*int", 
                  "Rel ~~ ", RelSI$eSImod, "*Rel",
                  "FM ~~ ",  FMSI$eSImod,  "*FM", 
                  "Clim ~~ ",ClimSI$eSImod,"*Clim",
                  "Vax ~~ ", VaxSI$eSImod, "*Vax" 
)
fitCorrel3 <- sem(modelCorrel3, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fitCorrel3,standardized=TRUE, fit.measures=TRUE)

#get correlation matrix for latent variables
lvcormat3 <- lavInspect(fitCorrel3, what = "cor.lv")
colnames(lvcormat3) <- rownames(lvcormat3) <- (c("Exceptionalism", "Nationalism", 
                                                 "Religiosity", "Free market", "Climate", 
                                                 "Vaccinations", "junk1", "junk2", "IQ Heritable"))
cormat3 <- stargazer(lvcormat3[,1:6], title="") #row 18 contains corrs with main factor in IQ model. Chop junk columns of subordinate factors
write.table(cormat3[20],file="_t.IQcor.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)

#now compute p-values
pvals2tailed3 <- pnorm(abs(inspect(fitCorrel3,what="est")$psi/inspect(fitCorrel3,what="se")$psi),lower.tail = FALSE)*2
twofacsigcorrs <- (colnames(lvcormat3)[1:6]) [pvals2tailed3[9,1:6]<.05]



#**-------------- now turn to specific prediction models


#**---------- IQ hierarchical two-factor: start with full model ...
fullmodel2fIQ <- c("
               IQ.two  ~ FMFac + intFac + RelFac + humFac
               humFac =~ hum
               intFac =~ int 
               RelFac =~ Rel
               FMFac =~ FM
               IQ.her =~ Q3.13 + Q3.17 + Q3.18 
               IQ.env =~ Q3.19 + Q3.14 + Q3.15 + Q3.16
               IQ.two =~ a*IQ.her + a*IQ.env
               
               hum ~~ ",  humSI$eSImod, "*hum",
                   "int ~~",  intSI$eSImod, "*int", 
                   "Rel ~~ ", RelSI$eSImod, "*Rel",
                   "FM ~~ ",  FMSI$eSImod,  "*FM" 
)
fitfull2fIQ <- sem(fullmodel2fIQ, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fitfull2fIQ,standardized=TRUE, fit.measures=TRUE)
modificationindices(fitfull2fIQ,sort.=TRUE,maximum.number=20)
full2fIQgof <- fitMeasures(fitfull2fIQ)


#.... smaller model for comparison ...
smmodel2fIQ <- c("
               IQ.two  ~ intFac 
               humFac =~ hum
               intFac =~ int
               RelFac =~ Rel
               FMFac =~ FM
               IQ.her =~ Q3.13 + Q3.17 + Q3.18 
               IQ.env =~ Q3.19 + Q3.14 + Q3.15 + Q3.16
               IQ.two =~ a*IQ.her + a*IQ.env
                 
                 hum ~~ ",  humSI$eSImod, "*hum",
                 "int ~~",  intSI$eSImod, "*int", 
                 "Rel ~~ ", RelSI$eSImod, "*Rel",
                 "FM ~~ ",  FMSI$eSImod,  "*FM" 
)
fitsm2fIQ <- sem(smmodel2fIQ, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fitsm2fIQ,standardized=TRUE, fit.measures=TRUE)
fitsm2fIQgof <- fitMeasures(fitsm2fIQ)
# .... which turns out to fit equally well ...
anova(fitfull2fIQ,fitsm2fIQ)
# .... so we can extract a small model without the other constructs
tinymodel2fIQ <- c("
               IQ.two  ~  intFac
               intFac =~ int
               IQ.her =~ Q3.13 + Q3.17 + Q3.18 
               IQ.env =~ Q3.19 + Q3.14 + Q3.15 + Q3.16
               IQ.two =~ a*IQ.her + a*IQ.env
                 
               int ~~",  intSI$eSImod, "*int"
                   
)
fittiny2fIQ <- sem(tinymodel2fIQ, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fittiny2fIQ,standardized=TRUE, fit.measures=TRUE)



#**---------- climate: start with full model...
bigmodelclim <- c("
               ClimFac  ~ FMFac + intFac + RelFac + humFac
               ClimFac =~ Clim               
               intFac =~ int
               RelFac =~ Rel
               humFac =~ hum
               FMFac =~ FM
               Rel ~~ ", RelSI$eSImod, "*Rel",
                  "hum ~~ ",  humSI$eSImod, "*hum",
                  "int ~~",  intSI$eSImod, "*int", 
                  "FM ~~ ",  FMSI$eSImod,  "*FM", 
                  "Clim ~~ ",ClimSI$eSImod,"*Clim"
)
fitbigmodclim <- sem(bigmodelclim, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fitbigmodclim,standardized=TRUE, fit.measures=TRUE)
fullclimgof <- fitMeasures(fitbigmodclim)
#.... eliminate predictors ....
modelclim <- c("
               ClimFac  ~ FMFac +  RelFac + humFac
               ClimFac =~ Clim               
               intFac =~ int
               FMFac =~ FM
               RelFac =~ Rel
               humFac =~ hum
               Rel ~~ ", RelSI$eSImod, "*Rel",
               "hum ~~ ",  humSI$eSImod, "*hum",
               "int ~~",  intSI$eSImod, "*int", 
               "FM ~~ ",  FMSI$eSImod,  "*FM", 
               "Clim ~~ ",ClimSI$eSImod,"*Clim"
)

fitmodclim <- sem(modelclim, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fitmodclim,standardized=TRUE, fit.measures=TRUE)
smclimgof <- fitMeasures(fitmodclim)
#parameterEstimates(fitmodclim, standardized=TRUE)
modindices(fitmodclim,sort. = TRUE, maximum.number = 4)
# ... which turns out not to make a difference ...
anova(fitbigmodclim,fitmodclim)

# and we have a final small model
smmodelclim <- c("
               ClimFac  ~ FMFac + intFac 
               ClimFac =~ Clim               
               intFac =~ int
               FMFac =~ FM
               int ~~",  intSI$eSImod, "*int", 
                 "FM ~~ ",  FMSI$eSImod,  "*FM", 
                 "Clim ~~ ",ClimSI$eSImod,"*Clim"
)
fitsmmodclim <- sem(smmodelclim, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fitsmmodclim,standardized=TRUE, fit.measures=TRUE)

# mediation?
mediateclimmodel <- c("
               humFac  ~ alpha1 * FMFac 
               intFac ~ alpha2 * FMFac

               ClimFac ~ direct * FMFac + beta1 *humFac + beta2 *intFac

              indirect1 := alpha1 * beta1
              indirect2 := alpha2 * beta2
              total:= indirect1 + indirect2 + direct
              proportion := (indirect1 + indirect2)/total

                  ClimFac =~ Clim               
                  humFac =~ hum
                  intFac =~ int
                  FMFac =~ FM
                  
                  hum ~~ ",  humSI$eSImod, "*hum",
                      "int ~~",  intSI$eSImod, "*int",     
                      "FM ~~ ",  FMSI$eSImod,  "*FM", 
                      "Clim ~~ ",ClimSI$eSImod,"*Clim"
)
fitmediateclimmodel <- sem(mediateclimmodel, compositeRelig, std.lv=TRUE, estimator="ML")
semPaths(fitmediateclimmodel, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE,
         structural=TRUE, layout="tree2",rotation=2)
summary(fitmediateclimmodel,standardized=TRUE, fit.measures=TRUE)
fullclimgof <- fitMeasures(fitmediateclimmodel)

# doube mediation for both climate and vax?
mediate2climmodel <- c("
                      ClimFac ~ d1 * FMFac + d2 * RelFac + b2 *humFac + b1 *intFac
                      VaxFac ~ d4 * FMFac + d3 * RelFac + b4 *humFac + b3 *intFac
                      humFac  ~ a4 * FMFac  + a2 * RelFac
                      intFac ~  a3 * FMFac + a1 * RelFac


                      indirectInt := a1 * b1 + a3 * b1 + a1 * b3 + a3 * b3 
                      indirecthum := a2 * b2 + a4 * b2 + a2 * b4 + a2 * b4 
                      total:= indirectInt + indirecthum + d1 + d2
                      proportion := (indirectInt + indirecthum)/total
                      
                      ClimFac =~ Clim               
                      humFac =~ hum
                      intFac =~ int
                      FMFac =~ FM
                      RelFac =~ Rel
                      VaxFac =~ Vax
                      
                      hum ~~ ",  humSI$eSImod, "*hum",
                       "int ~~",  intSI$eSImod, "*int",     
                       "FM ~~ ",  FMSI$eSImod,  "*FM", 
                       "Rel ~~ ", RelSI$eSImod, "*Rel",
                       "Clim ~~ ",ClimSI$eSImod,"*Clim",
                       "Vax ~~ ", VaxSI$eSImod, "*Vax"
)
fitmediate2climmodel <- sem(mediate2climmodel, compositeRelig, std.lv=TRUE, estimator="ML")
semPaths(fitmediate2climmodel, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE,
         structural=TRUE, layout="tree2",rotation=2)
summary(fitmediate2climmodel,standardized=TRUE, fit.measures=TRUE)
fullclimgof <- fitMeasures(fitmediate2climmodel)



#**------------ vax:start with full model
fullmodelvax <- c("
               VaxFac  ~ FMFac + intFac + RelFac + humFac
               VaxFac =~ Vax
               intFac =~ int
               RelFac =~ Rel
               FMFac =~ FM
               humFac =~ hum

               int ~~",  intSI$eSImod, "*int", 
                  "Rel ~~ ", RelSI$eSImod, "*Rel",
                  "hum ~~ ",  humSI$eSImod, "*hum",
                  "FM ~~ ",  FMSI$eSImod,  "*FM", 
                  "Vax ~~ ", VaxSI$eSImod, "*Vax"
)
fitfullmodvax <- sem(fullmodelvax, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fitfullmodvax,standardized=TRUE, fit.measures=TRUE)
fullvaxgof <- fitMeasures(fitfullmodvax)
#parameterEstimates(fitfullmodvax, standardized=TRUE)
modindices(fitfullmodvax, sort.=TRUE, maximum.number=4)
#.... now a smaller one ....
smmodelvax <- c("
               VaxFac  ~ FMFac + intFac + RelFac 
               VaxFac =~ Vax
               intFac =~ int
               RelFac =~ Rel
               FMFac =~ FM
                humFac =~ hum

               int ~~",  intSI$eSImod, "*int", 
                "Rel ~~ ", RelSI$eSImod, "*Rel",
                "hum ~~ ",  humSI$eSImod, "*hum",  
                "FM ~~ ",  FMSI$eSImod,  "*FM", 
                "Vax ~~ ", VaxSI$eSImod, "*Vax"
)
fitsmmodvax <- sem(smmodelvax, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fitsmmodvax,standardized=TRUE, fit.measures=TRUE)
smvaxgof <- fitMeasures(fitsmmodvax)
anova(fitfullmodvax,fitsmmodvax)


#**--------------- Now fit all scientific constructs together (two-factor IQ)
allfull2fmod <- c("
                ClimFac  ~ FMFac + intFac + RelFac + humFac
                VaxFac ~ FMFac + intFac + RelFac + humFac
                IQ.two ~  FMFac + intFac + RelFac + humFac
                
                humFac =~ hum
                intFac =~ int
                RelFac =~ Rel
                FMFac =~ FM
                ClimFac =~ Clim
                VaxFac =~ Vax
                IQ.her =~ Q3.13 + Q3.17 + Q3.18 
                IQ.env =~ Q3.19 + Q3.14 + Q3.15 + Q3.16
                IQ.two =~ a*IQ.her + a*IQ.env
                
                hum ~~ ",  humSI$eSImod, "*hum",
                  "int ~~",  intSI$eSImod, "*int", 
                  "Rel ~~ ", RelSI$eSImod, "*Rel",
                  "FM ~~ ",  FMSI$eSImod,  "*FM", 
                  "Clim ~~ ",ClimSI$eSImod,"*Clim",
                  "Vax ~~ ", VaxSI$eSImod, "*Vax" )

fitall2ffull  <- sem(allfull2fmod, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fitall2ffull,standardized=TRUE, fit.measures=TRUE)
#parameterEstimates(fitall2ffull, standardized=TRUE)
modindices(fitall2ffull, sort.=TRUE, maximum.number=4)
#.... and now a constrained model ....
modelall2f <- c("
               ClimFac  ~ FMFac + intFac 
               VaxFac ~ FMFac + intFac + RelFac + humFac
               IQ.two ~  intFac
               
               humFac =~ hum
               intFac =~ int
               RelFac =~ Rel
               FMFac =~ FM
               ClimFac =~ Clim
               VaxFac =~ Vax
               IQ.her =~ Q3.13 + Q3.17 + Q3.18 
               IQ.env =~ Q3.19 + Q3.14 + Q3.15 + Q3.16
               IQ.two =~ a*IQ.her + a*IQ.env
               
               hum ~~ ",  humSI$eSImod, "*hum",
                "int ~~",  intSI$eSImod, "*int", 
                "Rel ~~ ", RelSI$eSImod, "*Rel",
                "FM ~~ ",  FMSI$eSImod,  "*FM", 
                "Clim ~~ ",ClimSI$eSImod,"*Clim",
                "Vax ~~ ", VaxSI$eSImod, "*Vax" )

modelall2ffit <- sem(modelall2f, compositeRelig, std.lv=TRUE, estimator="ML")
summary(modelall2ffit,standardized=TRUE, fit.measures=TRUE)
anova(fitall2ffull,modelall2ffit)
finalmodelgof <- fitMeasures(modelall2ffit)
semPaths(modelall2ffit, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE,
         structural=TRUE, layout="tree2",rotation=2)

#** ---- now apportion variance -------------------------------------
# call to getr2(criterion,predictors,vars), 
#   e.g.: getr2("Clim", c("int", "FM"), c("Clim", "int", "Rel", "hum", "FM"))
# this decomposes the r^2 embedded in model fitall2ffull for the 3 latent criteria
# Note: this is 
climcomps <- getr2comps("Clim")
inspect(fitall2ffull,'r2')["ClimFac"]
str(climcomps)


vaxcomps  <- getr2comps("Vax")
inspect(fitall2ffull,'r2')["VaxFac"]
str(vaxcomps)
sum(unlist(vaxcomps))

IQcomps   <- getr2comps("IQ")
inspect(fitall2ffull,'r2')["IQ.two"]
str(IQcomps)

#rearrange data into a more suitable format
compnts <- cbind(unlist(climcomps),unlist(vaxcomps),unlist(IQcomps))[substr(names(climcomps),1,1)=="u",]
sharvar <- cbind(unlist(climcomps),unlist(vaxcomps),unlist(IQcomps))[substr(names(climcomps),1,1)!="u",]
compnts <- rbind(compnts,colSums(sharvar),cbind(sum(unlist(climcomps)),sum(unlist(vaxcomps)),sum(unlist(IQcomps))))
row.names(compnts) <- c("Nationalism","Religiosity","Exceptionalism","Free market","All shared","Total")

compmat <- stargazer(compnts, title="Decomposition of variance explained",digits=2,digits.extra=0) 
#compmat <- str_replace(compmat,"\\$-\\$","")
write.table(compmat[c(12,10,11,13,14,15)],file="_t.compnts.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)

#--- all constructs, predict all criteria from single 2nd-order conservatism factor --------------- 
all2ndorderfac <- c("
                ClimFac ~   C.two
                VaxFac  ~   C.two
                IQ.two  ~   C.two
                
                C.two =~ NA*FMFac + RelFac + intFac + humFac
                C.two ~~ 1*C.two 

                humFac =~ hum
                intFac =~ int
                RelFac =~ Rel
                FMFac =~ FM
                ClimFac =~ Clim
                VaxFac =~ Vax
                IQ.her =~ Q3.13 + Q3.17 + Q3.18 
                IQ.env =~ Q3.19 + Q3.14 + Q3.15 + Q3.16
                IQ.two =~ a*IQ.her + a*IQ.env
                humFac ~~ RelFac
                
                hum ~~ ",  humSI$eSImod, "*hum",
                    "int ~~",  intSI$eSImod, "*int", 
                    "Rel ~~ ", RelSI$eSImod, "*Rel",
                    "FM ~~ ",  FMSI$eSImod,  "*FM", 
                    "Clim ~~ ",ClimSI$eSImod,"*Clim",
                    "Vax ~~ ", VaxSI$eSImod, "*Vax" )

fit2ndorderfac  <- sem(all2ndorderfac, compositeRelig, std.lv=TRUE, estimator="ML")
summary(fit2ndorderfac,standardized=TRUE, fit.measures=TRUE)
modindices(fit2ndorderfac, sort.=TRUE, maximum.number=4)
