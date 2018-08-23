### Relate Hordeum per capita seed production to treatments ##

require(gdata)
library(tidyverse)
library(minpack.lm)

## Read in data and select relevant columns #

## 2015 greenhouse germinants
stem2016 <- read.xls("~/Dropbox/Carrizo-pop-models/Data/Carrizo 2016 Hordeum & Bromus data.xlsx", sheet=2, header=T, na.strings="#N/A!") %>%
  tbl_df() %>%
  select(Rep, Site, Troph, Eng, Quadrat, Grass, DTiller, GTiller, ITiller, Height_mean) %>%
  group_by(Site, Troph, Eng, Quadrat, Grass, Rep) %>%
  mutate(stemCount = as.numeric(n())) %>%
  mutate(stemCount = ifelse(is.na(DTiller) & is.na(GTiller) & is.na(ITiller), 0, stemCount))

seedallom2016 <- read.xls("~/Dropbox/Carrizo-pop-models/Data/Carrizo 2016 Hordeum & Bromus data.xlsx", sheet=1, header=T, na.strings="#N/A!") %>%
  tbl_df() %>%
  select(Hordeum,	Tiller.Old,	Tiller.New.Rep,	Height,	Seed,	Maturity,	Tiller.Immature) 

ggplot(subset(seedallom2016, Maturity != "im"), aes(x=Height, y=Seed)) + geom_point(aes(color = Maturity)) + 
  geom_smooth(method = "lm") + geom_smooth(aes(color = Maturity), method = "lm")

l <- lm(Height ~ Seed, data = subset(seedallom2016, Maturity != "im"))
summary(l)


stem2016_2 <- stem2016 %>%
  filter(Grass == "H") %>%
  mutate(seedpertiller = 1.42026 + 0.79031*Height_mean,
         seeds = seedpertiller*GTiller, 
         seeds = ifelse(is.na(seeds), 0, seeds))

stem2016_3 <- stem2016_2 %>%
  group_by(Site, Troph, Eng, Quadrat) %>%
  summarize(stemCount = mean(stemCount), seedCount = sum(seeds)) %>%
  mutate(perCap = seedCount/stemCount)


ggplot(stem2016_3, aes(x = stemCount, y=seedCount)) + geom_point() + facet_grid(Troph~Eng) + geom_abline(slope = 1)

