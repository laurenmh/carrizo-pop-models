require(gdata)
library(tidyverse)
library(minpack.lm)


# Data import and tidying ------------------------------------------------------

# long term cover
covall <-  read.xls("~/Dropbox/Carrizo-pop-models/Data/Carrizo 2009_2014 cover data.xlsx", sheet=1, header=T, na.strings="#N/A!") %>%
  tbl_df() 

covall2 <- covall %>%
  select(year, Precip.cm, Precip.lag, Site, GKR, Pasture, block, Site.num, Plot,
         Quad.ID, Precinct, erocic.c, hormur.c, grass.c, forb.c) 

names(covall2) = c("year", "Precip", "prevPrecip", "site", "Troph", "pasture",
                   "block", "sitenum", "plot", 'quadID', "Eng",
                   "erodium", "hordium", "grass", "forb")

# create previous year precip
covall3 <- covall2 %>%
  mutate(year = year + 1) %>%
  select(-Precip, -prevPrecip)

names(covall3) = c("year", "site", "Troph", "pasture",
                   "block", "sitenum", "plot", 'quadID', "Eng",
                   "preverodium", "prevhordium", "prevgrass", "prevforb")

covtog <- left_join(covall2, covall3) %>%
  filter(!is.na(prevforb)) %>%
  mutate(Eng = ifelse(Eng == "P", "burrow", "noburrow")) %>%
  mutate(Troph = ifelse(Troph == "present", "GKR", "noGKR")) %>%
  mutate(trt = as.factor(Precip)) %>%
  mutate(Trophbin = ifelse(Troph == "present", 1, 0))



# Preliminary visuals -----------------------------------------------------

ggplot(covtog, aes(x=log(preverodium + 1), y= log(erodium +1),color = as.factor(Precip))) + geom_point() + facet_grid(Troph~Eng)# + geom_smooth(se=F, method = "lm")
ggplot(covtog, aes(x=(preverodium ), y= (erodium),color = as.factor(prevPrecip))) + geom_point() + facet_grid(Troph~Eng) + geom_smooth(se=F, method = "lm")

ggplot(covtog, aes(x=log(prevhordium + 1), y= log(hordium +1), color = as.factor(prevPrecip))) + geom_point() + facet_grid(Troph~Eng) + geom_smooth(se=F, method = "lm")
ggplot(covtog, aes(x=(prevhordium ), y= (hordium), color = as.factor(prevPrecip))) + geom_point() + facet_grid(Troph~Eng) #+ geom_smooth(se=F, method = "lm")



# Erodium models ----------------------------------------------------------
# Might try to add precip as a covariate here and for Hordeum; if do that need to put it in the projections 
# Currently just running the model for each scenario (rat/no rat, engineer/no enginner, low/med/high precip)

m1 <- as.formula(log(erodium +1) ~  log(preverodium +1)*((lambda  )/(1+aiE*log(preverodium + 1) + aiH*log(prevhordium + 1))))

treatments <- unique(covtog$trt)

ERoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(ERoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  m1out <- nlsLM(m1, start=list(lambda=1, aiE = .01, aiH=.01), 
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(covtog, trt == treatments[i]  & Eng == "burrow"))
  summary(m1out)
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Erodium"
  ERoutput <- rbind(ERoutput, outreport)
}


# Hordeum models ----------------------------------------------------------
# Might try to add precip as a covariate here and for Hordeum; if do that need to put it in the projections 
# Currently just running the model for each scenario (rat/no rat, engineer/no enginner, low/med/high precip)

m1 <- as.formula(log(hordium +1) ~  log(prevhordium +1)*((lambda )/(1+aiE*log(preverodium + 1) + aiH*log(prevhordium + 1))))

treatments <- unique(covtog$trt)

HOoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(HOoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  m1out <- nlsLM(m1, start=list(lambda=1,  aiE = .01, aiH=.01), 
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(covtog, trt == treatments[i]  & Eng == "burrow"))
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



# Use parameters ----------------------------------------------------------
# Functions below to run the models under different scenarios; didn't get into running them yet

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





