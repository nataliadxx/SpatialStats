library(sp)
library(rgdal)
library(spdep)
library(tidyverse)
library(ggplot2)
library(mapview)

# Read the raw dataset provided by the author of the paper
LAhex <- readOGR("./Beyond405_Data/LAhex_ACS_MOVED_703_UTM11N.shp")

# Data preprocessing --------------------------------------------------------
# Exclude hexagons without COVID data or less than 1000 people
lahex.df <- data.frame(LAhex)
excl.id <- which(lahex.df$DP05_0001E < 1000 | is.na(lahex.df$crt630_mea) == T)
LAhex <- LAhex[-excl.id,]
# see the Readme file attached with the data for the field meaning

# Exclude hexagons without contiguous neighbors
QN <- poly2nb(LAhex, queen = T) # queen neighbors for each polygon
cards <- card(QN) # number of neighbors
no.neighbor.id <- which(cards == 0)
LAhex <- LAhex[-no.neighbor.id,] # 184 hexagons left

# Create the dataframe with response and independent variables
lahex.df <- data.frame(age18 = LAhex$DP05_0019P, age65 = LAhex$DP05_0024P,
                       latino = LAhex$DP05_0071P, white = LAhex$DP05_0037P,
                       black = LAhex$DP05_0038P, asian = LAhex$DP05_0044P,
                       poverty = LAhex$DP03_0119P, uninsured = LAhex$DP03_0099P,
                       bachelor = LAhex$DP02_0067P, pop.tot = LAhex$DP05_0001E,
                       hh_tot = LAhex$DP05_0086E)
lahex.df <- cbind(lahex.df, fid = LAhex$fid, area = LAhex$areasqkm, 
                  crt = LAhex$crt630_mea, adjcrt = LAhex$adjcrt630_, 
                  trt = LAhex$trt630_mea, adjtrt = LAhex$adjtrt630_,
                  westside = LAhex$Westside)
lahex.df <- lahex.df %>%
  mutate(pop.dens = pop.tot/10) %>% # calculate population density
  mutate(hh.dens = hh_tot) %>% # calculate household density
  mutate(prt = 100*crt/trt) # calculate crude positivity rates
write.csv(lahex.df, "variables.csv", row.names = F) 
LAhex$prt <- lahex.df$prt # add the calculated field to the shp file

# Check if the data ranges align with the results presented in the paper
summary(lahex.df$adjtrt)
summary(lahex.df$adjcrt)
summary(lahex.df$prt)
# the descriptive characteristics match Fig.1 top row

# Plotting Fig.1 top row ---------------------------------------------------
# Plot the quintile maps of testing, diagnosis and positivity rates (?)
ratecolor <- rev(RColorBrewer::brewer.pal(7, "PRGn"))
LAhex.sf <- st_as_sf(LAhex)
ggplot(LAhex.sf, aes(fill=adjtrt630_)) +
  geom_sf() +
  scale_fill_manual(values = ratecolor) +
  theme_minimal() +
  labs(title = "Age-adjusted testing rates (per 100000 population)")

# Producing Table 1 ----------------------------------------------------------
lahex.df <- lahex.df %>%
  mutate(prt.round = round(prt, 1)) %>%
  mutate(prt_level = NA)
# divide hexagons into 3 groups according to the positivity rate
for (i in 1:nrow(lahex.df)) {
  if (lahex.df$prt[i] < 5) lahex.df$prt_level[i] <- "Low"
  if (lahex.df$prt[i] >= 5 & lahex.df$prt[i] < 10) lahex.df$prt_level[i] <- "Med"
  if (lahex.df$prt[i] >= 10) lahex.df$prt_level[i] <- "High"
}

# generate the descriptive statistics of the independent variables
data.sum <- lahex.df %>%
  group_by(factor(prt_level)) %>%
  summarise(no.hex = length(prt_level),
            mean_age18 = mean(age18), sd_age18 = sd(age18),
            mean_age65 = mean(age65), sd_age65 = sd(age65),
            mean_white = mean(white), sd_white = sd(white),
            mean_black = mean(black), sd_black = sd(black),
            mean_asian = mean(asian), sd_asian = sd(asian),
            mean_latino = mean(latino), sd_latino = sd(latino),
            mean_poverty = mean(poverty), sd_poverty = sd(poverty),
            mean_uninsured = mean(uninsured), sd_uninsured = sd(uninsured),
            mean_bachelor = mean(bachelor), sd_bachelor = sd(bachelor),
            mean_pop.dens = mean(pop.dens), sd_pop.dens = sd(pop.dens),
            mean_hh.dens = mean(hh.dens), sd_hh.dens = sd(hh.dens))
write.csv(data.sum, "Table1_datasummary.csv", row.names = F)
# not identical with Table 1, either rounding positivity rate or not
# household density?
# correlation analysis?

# Geostatistical analysis, LISA ----------------------------------------------------
# Create spatial weights matrix
QN <- poly2nb(LAhex, queen = T)
QN1.lw <- nb2listw(QN, style = "B")

# Calculate Local Moran's I (LISA) for the three variables
I.local.test <- localmoran(LAhex$adjtrt630_, QN1.lw)
lahex.df$test.I <- I.local.test[,1]
lahex.df$test.p <- I.local.test[,5]
I.local.case <- localmoran(LAhex$adjcrt630_, QN1.lw)
lahex.df$case.I <- I.local.case[,1]
lahex.df$case.p <- I.local.case[,5]
I.local.positiv <- localmoran(LAhex$prt, QN1.lw)
lahex.df$positiv.I <- I.local.positiv[,1]
lahex.df$positiv.p <- I.local.positiv[,5]
#which(lahex.df$test.p<0.05)

# Calculate Local G* for the three variables
G.local.test <- localG(LAhex$adjtrt630_, QN1.lw)
lahex.df$G.test <- as.numeric(G.local.test)
G.local.case <- localG(LAhex$adjcrt630_, QN1.lw)
lahex.df$G.case <- G.local.case
G.local.positiv <- localG(LAhex$prt, QN1.lw)
lahex.df$G.positiv <- G.local.positiv
#which(lahex.df$G.positiv < -1.96)

# Identify clusters and write in the shp. file
LAhex$test.c <- NA
LAhex$case.c <- NA
LAhex$positiv.c <- NA
for (i in 1:nrow(lahex.df)) {
  if (lahex.df$test.I[i] >0 & lahex.df$test.p[i]<0.05) LAhex$test.c[i] <- "High"
  if (lahex.df$case.I[i] >0 & lahex.df$case.p[i]<0.05) LAhex$case.c[i] <- "High"
  if (lahex.df$positiv.I[i] >0 & lahex.df$positiv.p[i]<0.05) LAhex$positiv.c[i] <- "High"
}
summary(LAhex$test.c)
# for local G, we ought to have G>1.96 or G<-1.96 to ensure significance, right?

# Plotting Fig.1 bottom row ---------------------------------------------------
LAhex.sf <- st_as_sf(LAhex)
ggplot(LAhex.sf, aes(fill = factor(test.c))) +
  geom_sf() +
  scale_fill_manual(values = c("Tomato", "SteelBlue")) +
  theme_minimal() +
  labs(title="Testing rate clusters")
summary(LAhex.sf$test.c)

# Spatially lagged model ---------------------------------------------------
# Center the variables to a mean of 0 and scale to an sd of 1
lahex.center <- lahex.df %>%
  mutate_at(c("adjtrt", "adjcrt", "prt", "age18", "age65", "latino", "white", 
              "black", "asian", "poverty", "uninsured", "bachelor", "pop.dens", 
              "hh.dens"), ~(scale(.) %>% as.vector))

# Name the regression equation
reg.eq1 <- prt ~ age18 + age65 + latino + white + black + asian + 
  poverty + uninsured + bachelor + pop.dens + hh.dens

# Y lagged -  Spatial Autoregressive Model (SAR) 
SReg.SAR = lagsarlm(reg.eq1, data = lahex.center, QN1.lw)
summary(SReg.SAR)
summary(impacts(SReg.SAR, listw = QN1.lw, R = 500), zstats = TRUE)





