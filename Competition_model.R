require(gdata)
library(dplyr)
library(tidyr)
library(ggplot2)
library(minpack.lm)

## Read in data and select relevant columns #

## 2015 greenhouse germinants
gh2015 <- read.xls("Carrizo 2016 Hordeum & Bromus data.xlsx", sheet=4, header=T, na.strings="#N/A!") %>%
  tbl_df() %>%
  filter(N.trt == 0) %>%
  select(Pasture, Site, Plot.ID, Sample.ID, Trophic.trt, Eng.trt, Precip.trt,  W.trt, tGrass, tForb)

gh2016 <- read.xls("Carrizo 2016 Hordeum & Bromus data.xlsx", sheet=3, header=T, na.strings="#N/A!") %>%
  tbl_df() %>%
  filter(Season == 0) %>%
  select(Pasture, Site, Plot, Sample, Troph, Eng, Precip, Water, Grass.all, Forb.all)

## standardize names
mynames <- names(gh2016)
names(gh2015) <- mynames
names(gh2015)[9:10] = c(  "Ho2015", "Er2015")
names(gh2016)[9:10] = c(  "Ho2016", "Er2016")

# put it all together
dat <- left_join(gh2015, gh2016) %>%
  filter(Precip != 1) %>%
  filter(Water == 1) %>%
  mutate(trt = interaction(Precip, Troph)) 

ggplot(dat, aes(x=log(Er2015 + 1), y= log(Er2016 +1),color = as.factor(Precip))) + geom_point() + facet_grid(Troph~Eng)
ggplot(dat, aes(x=log(Ho2015 + 1), y= log(Ho2016 +1), color = as.factor(Precip))) + geom_point() + facet_grid(Troph~Eng)
ggplot(dat, aes(x=(Er2015 + 1))) + geom_histogram() + facet_grid(interaction(Troph, Eng)~Precip)
ggplot(dat, aes(x=(Ho2015 + 1))) + geom_histogram() + facet_grid(interaction(Troph, Eng)~Precip)


##### ERODIUM MODELS ####
## Added in precip as a covariate here and for Hordeum; if keep that need to put it in the projections ##
m1 <- as.formula(log(Er2016 +1) ~  log(Er2015 +1)*((lambda  )/(1+aiE*log(Er2015 + 1) + aiH*log(Ho2015 + 1))))

treatments <- unique(dat$trt)

ERoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(ERoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  m1out <- nlsLM(m1, start=list(lambda=1, aiE = .01,  aiH=.01), 
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(dat, trt == treatments[i] ))
  summary(m1out)
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Erodium"
  ERoutput <- rbind(ERoutput, outreport)
}


### HORDEUM MODEL ###
m1 <- as.formula(log(Ho2016 +1) ~  log(Ho2015 +1)*((lambda )/(1+aiE*log(Er2015 + 1) + aiH*log(Ho2015 + 1))))

treatments <- unique(dat$trt)

HOoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(HOoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  m1out <- nlsLM(m1, start=list(lambda=1,  aiE = .01, aiH=.01), 
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(dat, trt == treatments[i]))
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Hordeum"
  HOoutput <- rbind(HOoutput, outreport)
}


## PUT THE TWO OUTPUTS TOGETHER
model.dat <- rbind(ERoutput, HOoutput) %>%
  tbl_df() %>%
  select(estimate, params, treatment, species) %>%
  spread(params, estimate)


### CREATE A FUNCTION THAT RUNS THE MODEL
growth = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Nh[i+1] = exp(log(N$Nh[i]+1)*par.dat$Hlambda[i]/(1 + par.dat$HaiE[i]*log(N$Ne[i] + 1) +  par.dat$HaiH[i]*log(N$Nh[i] + 1))) -1
      N$Na[i+1] =ifelse(N$Na[i +1] == 1, 0, N$Na[i + 1])
    N$Ne[i+1] = exp(log(N$Ne[i] + 1)*par.dat$Elambda[i]/(1 + par.dat$EaiE[i]*log(N$Ne[i] + 1) + par.dat$EaiH[i]*log(N$Nh[i] + 1))) -1
     N$Ne[i+1] =ifelse(N$Ne[i +1] == 1, 0, N$Ne[i + 1])
  }
  return(N)
}



### CREATE A FUNCTION THAT CALCULATES POPULATION CHANGE OVER TIME AND GRWR ####

grwr = function(par.dat, t.num) {
  N1 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N1) = c("Nh", "Ne")
  N1[1,] = c(1,0)
  Hequil <- growth(N1, par.dat, t.num)
  
  N2 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N2) = c("Nh", "Ne")
  N2[1,] = c(0,1)
  Eequil <- growth(N2, par.dat, t.num)
  
  N3 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N3) = c("Nh", "Ne")
  N3[1,] = c(Hequil[t.num, 1],1)
  Einvade <- growth(N3, par.dat, t.num)
  Egrwr <- Einvade[2,2]/Einvade[1,2]
  Einvade$invader <- "Erodium"
  Einvade$grwr <- Egrwr
  Einvade$time <- as.numeric(row.names(Einvade))
  
  N4 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N4) = c("Nh", "Ne")
  N4[1,] = c(1, Eequil[t.num, 2])
  Hinvade <- growth(N4, par.dat, t.num)
  
  Hgrwr <- Ainvade[2,1]/Ainvade[1,1]
  
  Hinvade$invader <- "Hordeum"
  Hinvade$grwr <- Hgrwr
  Hinvade$time <- as.numeric(row.names(Hinvade))
  
  
  out <- rbind(Hinvade, Einvade) 
  
  return(out)
}  


### CREATE A FUNCTION THAT SETS CONSISTENT PARAMETERS ACROSS TIMESTEPS ####
consistent.par <- function(model.dat, j.num, t.num){
  par = as.data.frame(matrix(NA, nrow=t, ncol=6))
  colnames(par) = c("Hlambda", "HaiE", "HaiH", "Elambda", "EaiE", "EaiH")
  trtselect = rep(j.num, t.num)
  
  for (i in 1:(t.num)){
    myparams <- subset(model.dat, treatment == treatments[trtselect[i]])
    par$Hlambda[i] <- as.numeric(subset(myparams, species == "Hordeum")[5])
    par$HaiE[i] <- as.numeric(subset(myparams, species == "Hordeum")[4])
    par$HaiH[i] <- as.numeric(subset(myparams, species == "Hordeum")[3])
    par$Elambda[i] <- as.numeric(subset(myparams, species == "Erodium")[5])
    par$EaiE[i] <- as.numeric(subset(myparams, species == "Erodium")[4])
    par$EaiH[i] <- as.numeric(subset(myparams, species == "Erodium")[3])
  }
  return(par)
}


t = 100

# troph, eng
treatments
consistentPar <- consistent.par(model.dat, 3, t)

N = as.data.frame(matrix(NA, nrow=t, ncol=2))
colnames(N) = c("Nh", "Ne")
N[1,] =c(1, 1)


growth(N, consistentPar, t)
