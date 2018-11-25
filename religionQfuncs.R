#**--- package of functions for religionQualtrics model ----------------

#function to obtain r^2 from predictive model
getr2 <- function(criterion,predictors,vars) {
  criterion <- paste(criterion,"Fac",sep="")
  modelcode <- paste(criterion, " ~ ", paste(paste(predictors,"Fac",sep=""), collapse=" + "), " \n ", sep="")
  modelcode <- paste(modelcode, paste(paste(vars,"Fac"," =~ ",vars," \n ",sep=""),collapse=" "))
  if (criterion=="IQFac") {
    modelcode <- paste(modelcode, "IQ.her =~ Q3.13 + Q3.17 + Q3.18 \n ",
                    "IQ.env =~ Q3.19 + Q3.14 + Q3.15 + Q3.16 \n ",
                    "IQFac =~ a*IQ.her + a*IQ.env \n ",sep="")
  }
  for (v in vars){
    modelcode <- paste(modelcode, v," ~~ ",eSImods[v], "*",v," \n ",sep="")
  }
  modelfit <- sem(modelcode, compositeRelig, std.lv=TRUE, estimator="ML")
  #https://groups.google.com/forum/#!topic/lavaan/W5rIa2eo3uQ <-- explains r-squared for latent vars
  return( inspect(modelfit, 'r2')[criterion])
}

#function to fit and print a measurement model
fitMM <- function (humspecmod,thisrelig) {
  fithumspec <- sem(humspecmod, data=thisrelig)
  summary(fithumspec, standardized=TRUE, fit.measures=TRUE)
  #parameterEstimates(fithumspec, standardized=TRUE)
  mod_ind <- modificationindices(fithumspec)
  print(head(mod_ind[order(mod_ind$mi, decreasing=TRUE), ], 4))
  return(fitMeasures(fithumspec))
}

#function to create and run a single indicator model and return omega and so on
#  indicators are provided in character vector
#  pairwise correlations are optionally provided as a list of pairwise character vectors
singleindmodel <- function(indicators,paircorrs,dat) {
  pcs <- paste(unlist(lapply(paircorrs,FUN=function(x) paste(x[1], "~~", x[2], "\n"))),collapse=" ")
  SImod <- paste("factor =~", paste(indicators, collapse=" + "), "\n",
                 pcs, 
                 "phantfac <~", paste(paste("1*",indicators,sep=""), collapse=" + "), "\n",
                 "factor ~~ 0*phantfac")
  
  fitSImod <- cfa(SImod, dat[,indicators], estimator="ML")
  ParSImod <- parameterEstimates(fitSImod, standardized=TRUE)
  LoadingsSImod <- ParSImod[1:length(indicators), "std.all"]
  ErrorvarSImod <- 1 - LoadingsSImod^2
  ImpliedCorrSImod <- lavTech(fitSImod, what='cor.lv')
  OmegaSImod <- ImpliedCorrSImod[[1]][1,2]^2            #Squared correlation between factor and phantom variable
  varSImod <- var(apply(dat[,indicators], MARGIN=1, FUN=mean)) #variance of composite
  SDSImod <- sqrt(varSImod) #SD of composite, to pass back for analysis
  eSImod <- (1-OmegaSImod)*varSImod       #error term of single indicator
  return(listN(OmegaSImod,eSImod,SDSImod))
}

#http://stackoverflow.com/questions/21011672/automatically-add-variable-names-to-elements-of-a-list
listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}
revscore <- function (x,mm) {  #this reverse scores a scale for items of reversed polarity
  return ((mm+1)-x)            #same as fixscore but to be conceptually clear it's a separate function
}
fixscore <- function (x,mm) {  #this fixes scale so SD=1 and SA=7 irrespective of polarity
  return ((mm+1)-x)
}


#--- function to compute variance share ----------------------------
getr2comps <- function(criterion) {
  #define all potential predictors and set aside space
  ijkl <- c("int", "Rel", "hum", "FM")
  names(ijkl) <- c("i","j","k","l")
  comps <- r2 <- vector("list",0)
  if (criterion=="IQ") {
    criterion2p <- NULL
  } else {
    criterion2p <- criterion
  }
  
  #unique contributions of each factor
  r2$ijkl <- getr2(criterion, ijkl,      c(criterion2p,ijkl))
  r2$jkl  <- getr2(criterion, ijkl[c("j","k","l")], c(criterion2p,ijkl[c("j","k","l")]))
  r2$ikl  <- getr2(criterion, ijkl[c("i","k","l")], c(criterion2p,ijkl[c("i","k","l")]))
  r2$ijl  <- getr2(criterion, ijkl[c("i","j","l")], c(criterion2p,ijkl[c("i","j","l")]))
  r2$ijk  <- getr2(criterion, ijkl[c("i","j","k")], c(criterion2p,ijkl[c("i","j","k")]))
  comps$uint <- r2$ijkl-r2$jkl  #i  int
  comps$uRel <- r2$ijkl-r2$ikl  #j  Rel
  comps$uhum <- r2$ijkl-r2$ijl  #k  hum
  comps$uFM  <- r2$ijkl-r2$ijk  #l  FM
  
  #pairwise shared contributions
  r2$kl  <- getr2(criterion, ijkl[c("k","l")], c(criterion2p,ijkl[c("k","l")]))
  r2$jl  <- getr2(criterion, ijkl[c("j","l")], c(criterion2p,ijkl[c("j","l")]))
  r2$jk  <- getr2(criterion, ijkl[c("j","k")], c(criterion2p,ijkl[c("j","k")]))
  r2$il  <- getr2(criterion, ijkl[c("i","l")], c(criterion2p,ijkl[c("i","l")]))
  r2$ik  <- getr2(criterion, ijkl[c("i","k")], c(criterion2p,ijkl[c("i","k")]))
  r2$ij  <- getr2(criterion, ijkl[c("i","j")], c(criterion2p,ijkl[c("i","j")]))
  
  comps$intRel <- -r2$kl + r2$ikl + r2$jkl - r2$ijkl  #ij
  comps$inthum <- -r2$jl + r2$ijl + r2$jkl - r2$ijkl  #ik
  comps$intFM  <- -r2$jk + r2$ijk + r2$jkl - r2$ijkl  #il
  comps$Relhum <- -r2$il + r2$ijl + r2$ikl - r2$ijkl  #jk
  comps$RelFM  <- -r2$ik + r2$ijk + r2$ikl - r2$ijkl  #jl
  comps$humFM  <- -r2$ij + r2$ijk + r2$ijl - r2$ijkl  #kl
  
  #triplet shared
  r2$i  <- getr2(criterion, ijkl["i"], c(criterion2p,ijkl["i"]))
  r2$j  <- getr2(criterion, ijkl["j"], c(criterion2p,ijkl["j"]))
  r2$k  <- getr2(criterion, ijkl["k"], c(criterion2p,ijkl["k"]))
  r2$l  <- getr2(criterion, ijkl["l"], c(criterion2p,ijkl["l"]))
  
  comps$intRelhum <- -r2$l + r2$il + r2$jl + r2$kl - r2$ijl - r2$ikl - r2$jkl + r2$ijkl
  comps$intRelFM  <- -r2$k + r2$ik + r2$jk + r2$kl - r2$ijk - r2$ikl - r2$jkl + r2$ijkl
  comps$inthumFM  <- -r2$j + r2$ij + r2$jk + r2$jl - r2$ijk - r2$ijl - r2$jkl + r2$ijkl
  comps$RelhumFM  <- -r2$i + r2$ij + r2$ik + r2$il - r2$ijk - r2$ijl - r2$ikl + r2$ijkl
  
  #all shared
  comps$intRelhumFM <- r2$i + r2$j + r2$k + r2$l - r2$ij - r2$ik - r2$il - r2$jk  - r2$jl -
    r2$kl + r2$ijk + r2$ijl + r2$ikl + r2$jkl - r2$ijkl
  return(comps)
}