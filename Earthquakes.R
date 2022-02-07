# DSO110 Final Project
# The Science Explorers
# Borly Green, Jessenia Lorenzo, and Melissa Gonzalez
# Earthquake Season
# ------------------------------------------------------------------------------

# Evaluation Question
# Does season have any impact on either the frequency of earthquakes or the
#   magnitude of earthquakes?

# IV - Season (Categorical)
#       Calculated from latitude/hemisphere and date

# DV1 - Magnitude (Continuous)
# DV2 - Count of events/Frequency (Categorical?)

# ------------------------------------------------------------------------------

# Setup

library('tidyverse')
library('data.table')
library('dplyr')
library('IDPmisc')
library("tidyr")
library("mvnormtest") #for MANOVA
library("car") # for test of homogeneity of variance
library('rcompanion') # for plotting normality of histograms
library("caret") # linear regressions
library("gvlma")
library("predictmeans")
library("e1071")
library("lmtest")

library(readxl)

library('geosphere')
library('rgeos')
library('tmap')

library(sf) #to read xml/kml/xsd files


# Import Data

database <- read.csv("~/Desktop/Desktop Entity/DSO110 Final Group Project/Earthquakes/database 2.csv")
View(database)

# Data Wrangling

# Subset Data
earthquakes <- database[c("Date", "Time", "Latitude", "Longitude", "Type", "Depth", "Magnitude", "ID")]
View(earthquakes)

# To get an overview of columns and data types
str(earthquakes)

# Remove Missing data - not needed
NaRV.omit(earthquakes)

# Convert Date from character to date
earthquakes$Date <- as.Date(earthquakes$Date, format = "%m/%d/%Y")
str(earthquakes)

earthquakes$MonthDay <- format(earthquakes$Date, format = "%m/%d" )
view(earthquakes)

# Assign Northern/Southern Hemisphere

# Hemispheres
#   Northern Hemisphere >  23.5
#   Southern Hemisphere < -23.5

earthquakes$Hemisphere <- ifelse(earthquakes$Latitude >= 23.5, "Northern",
                     ifelse(earthquakes$Latitude <= -23.5, "Southern", "Tropics"))
View(earthquakes)

# Assign New Column to show more detail in relation to equator and tropics of capricorn and cancer.

earthquakes$TropicZone <- ifelse(earthquakes$Latitude >= 23.5, "Northern",
                                 ifelse(earthquakes$Latitude <= -23.5, "Southern",
                                        ifelse(earthquakes$Latitude > 0 & earthquakes$Latitude < 23.5, "Tropic of Cancer",
                                               ifelse(earthquakes$Latitude < 0 & earthquakes$Latitude > -23.5, "Tropic of Capricorn", 
                                                      ifelse(earthquakes$Latitude == 0 , "Equator", "Error")))))

# Just to make sure that all data points were calculated correctly
unique(earthquakes$TropicZone)

# Assign MonthDay to Season

# Seasons
#   Fall(NH)/Spring(SH):    Sept 21  - Dec 20
#   Winter(NH)/Summer(SH):  Dec 21   - March 20
#   Spring(NH)/Fall(SH):    March 21 - June 20
#   Summer(NH)/Winter:      June 21 - Sept 20

earthquakes$Season <- ifelse(
  earthquakes$Hemisphere == "Northern" & earthquakes$MonthDay >= "09/21" & earthquakes$MonthDay <= "12/20", "Fall",
  ifelse(
    earthquakes$Hemisphere == "Northern" & earthquakes$MonthDay >= "12/21" & earthquakes$MonthDay <= "12/31", "Winter", 
    ifelse(
      earthquakes$Hemisphere == "Northern" & earthquakes$MonthDay >= "01/01" & earthquakes$MonthDay <= "03/20", "Winter", 
      ifelse(
        earthquakes$Hemisphere == "Northern" & earthquakes$MonthDay >= "03/21" & earthquakes$MonthDay <= "06/20", "Spring", 
        ifelse(
          earthquakes$Hemisphere == "Northern" & earthquakes$MonthDay >= "06/21" & earthquakes$MonthDay <= "09/20", "Summer", 
          ifelse(
            earthquakes$Hemisphere == "Southern" & earthquakes$MonthDay >= "09/21" & earthquakes$MonthDay <= "12/20", "Spring",
            ifelse(
              earthquakes$Hemisphere == "Southern" & earthquakes$MonthDay >= "12/21" & earthquakes$MonthDay <= "12/31", "Summer",
              ifelse(
                earthquakes$Hemisphere == "Southern" & earthquakes$MonthDay >= "01/01" & earthquakes$MonthDay <= "03/20", "Summer",
                ifelse(
                  earthquakes$Hemisphere == "Southern" & earthquakes$MonthDay >= "03/21" & earthquakes$MonthDay <= "06/20", "Fall",
                  ifelse(
                    earthquakes$Hemisphere == "Southern" & earthquakes$MonthDay >= "06/21" & earthquakes$MonthDay <= "09/20", "Winter",
                    "Tropical"))))))))))

View(earthquakes)

# Subset without Tropical Data
NSearthquakes <- earthquakes[earthquakes$Hemisphere != "Tropical",]
View(NSearthquakes)

# Dummy code "Season"
NSearthquakes$SeasonR[NSearthquakes$Season == "Winter"] <- 0
NSearthquakes$SeasonR[NSearthquakes$Season == "Spring"] <- 1
NSearthquakes$SeasonR[NSearthquakes$Season == "Summer"] <- 2
NSearthquakes$SeasonR[NSearthquakes$Season == "Fall"] <- 3

# Total Count of Earthquake Data by ID
n_distinct(NSearthquakes$ID)


# Count of Earthquakes by Season 
NSearthquakes %>% group_by(Season) %>% summarize(count = n())

#1 Fall    2551
#2 Spring  2640
#3 Summer  2554
#4 Winter  2622
#5 NA         2

# Remove NAs
NSearthquakes <-na.omit(NSearthquakes, c("Season"))
View(NSearthquakes)

# Count of Earthquakes by Season 
NSearthquakes %>% group_by(Season) %>% summarize(count = n())

#1 Fall    2551
#2 Spring  2640
#3 Summer  2554
#4 Winter  2622

# Set new column for counts for future analysis
NSearthquakes$Count <- 1
earthquakes$Count <- 1

# Export wrangled dataframes to use in Tableau
write.csv(NSearthquakes, "~/Desktop/Desktop Entity/NSearthquakes.csv")
write.csv(earthquakes, "~/Desktop/Desktop Entity/earthquakes.csv")

#-----------------------------------------------
# Analysis of Number of Earthquakes by season
# Poisson Regression

# Y = number of earthquakes
# X = season category

PoissonModel <- glm(count~Season, dat = EarthquakeSeasonMANOVA, family = "poisson")
summary(PoissonModel)

# Deviance Residuals: 
#   [1]  0  0  0  0

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  7.844241   0.019799 396.192   <2e-16 ***
# SeasonSpring 0.034293   0.027763   1.235    0.217    
# SeasonSummer 0.001175   0.027992   0.042    0.967    
# SeasonWinter 0.027452   0.027810   0.987    0.324    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# (Dispersion parameter for poisson family taken to be 1)

# Null deviance:  2.4411e+00  on 3  degrees of freedom
# Residual deviance: -3.3107e-13  on 0  degrees of freedom
# AIC: 46.792

# Number of Fisher Scoring iterations: 2

# Exponentiate the coefficients to get ratios
exp(coef(PoissonModel))
# (Intercept) SeasonSpring SeasonSummer SeasonWinter 
# 2551.000000     1.034888     1.001176     1.027832

# However, does not pass assumption of linear relationship
# -------------------------------------------------------------------
# ChiSquare to see if frequency is different than expected
# H0 = Frequency of each season is the same
# Ha = Frequency of each season is not the same

#1 Fall    2551
#2 Spring  2640
#3 Summer  2554
#4 Winter  2622

ObservedBySeason = c(2551, 2640, 2554, 2622)
HOBySeason = c(.25, .25, .25, .25)
chisq.test(x = ObservedBySeason, p = HOBySeason)

# Chi-squared test for given probabilities

# data:  ObservedBySeason
# X-squared = 2.4419, df = 3, p-value = 0.4859
# p-value is > 0.05, so we cannot reject our null hypothesis. There is 
# not a differece of frequency than equal frequency across seasons.

# ------------------------------------------------------------------------------
# Data Wrangling for ANOVA
# Make new dataframe for ANOVA
EarthquakeSeasonEvents <- data.frame(NSearthquakes$Season, NSearthquakes$Count)
View(EarthquakeSeasonEvents)

# Rename Columns
names(EarthquakeSeasonEvents)[names(EarthquakeSeasonEvents) == "NSearthquakes.Season"] <- "Season"
names(EarthquakeSeasonEvents)[names(EarthquakeSeasonEvents) == "NSearthquakes.Count"] <- "Count"

# Make tibble for count
EarthquakeSeasonCount <- NSearthquakes %>% group_by(Season) %>% summarize(count = n())
# Turn count tibble into dataframe
as.data.frame(EarthquakeSeasonCount)
View(EarthquakeSeasonCount)

# ------------------------------------------------------------------------------
# ANOVA test for normality

# Analysis
# Does the season impact either the frequency or magnitude of earthquakes?
#   IV - Season - Categorical
#   DV1 - Count of earthquakes - Continuous
#   ANOVA

# Uses packages loaded above
#   library("mvnormtest")
#   library("car")

# Test Assumptions
# Normality of EarthquakeSeasonMANOVA$count
plotNormalHistogram(EarthquakeSeasonCount$count)

# Count skewed slightly positive, so let's square root it it
EarthquakeSeasonCount$countSQRT <- sqrt(EarthquakeSeasonCount$count)
plotNormalHistogram(EarthquakeSeasonCount$countSQRT)

# Not much better, let's log it
EarthquakeSeasonCount$countLOG <- log(EarthquakeSeasonCount$count)
plotNormalHistogram(EarthquakeSeasonCount$countLOG)

# We will use the log of EarthquakeSeasonCount$count

# ------------------------------------------------------------------------------
# Test for Homogeneity of Variance
# We corrected for normality, so we will use Bartlett's Test
View(EarthquakeSeasonCount)
bartlett.test(Count ~ Season, data = EarthquakeSeasonEvents)
#Bartlett test of homogeneity of variances

#data:  Count by Season
#Bartlett's K-squared = NaN, df = 3, p-value = NA


# ------------------------------------------------------------------------------


#Combine dataframes for MANOVA
EarthquakeSeasonMANOVA <- cbind(EarthquakeSeasonCount,EarthquakeSeasonMagnitude,EarthquakeSeasonMagnitudeMedianR)
View(EarthquakeSeasonMANOVA)
str(EarthquakeSeasonMANOVA)

# Check if count and Magnitude are numeric
str(EarthquakeSeasonMANOVA$count)
str(EarthquakeSeasonMANOVA$Magnitude)
str(EarthquakeSeasonMANOVA$Median)

# Keep only two dependent variables, Count and Magnitude
keeps <- c("count", "Magnitude")
EarthquakeSeasonMANOVADVs <- EarthquakeSeasonMANOVA[keeps]
View(EarthquakeSeasonMANOVADVs)

# Format EarthquakeSeasonMANOVADVs as a Matrix
EarthquakeSeasonMatrix <- as.matrix(EarthquakeSeasonMANOVADVs)

# Analysis
# Does the season impact either the frequency or magnitude of earthquakes?
#   IV - Season - Categorical
#   DV1 - Count of earthquakes - Continuous
#   DV2 - Magnitude of earthquakes - Continuous
#     MANOVA

# Uses packages loaded above
#   library("mvnormtest")
#   library("car")

# Test Assumptions
# Sample Size
# More than 20 cases in each section

# Multivariate normality
mshapiro.test(t(EarthquakeSeasonMatrix))

# Shapiro-Wilk normality test

# data:  Z
# W = 0.94587, p-value = 0.6904

# p-value is not significant at p > 0.05

## Look at individual DVs for normality
# Normality of EarthquakeSeasonMANOVA$count
plotNormalHistogram(EarthquakeSeasonMANOVA$count)

# Count skewed slightly positive, so let's square root it
EarthquakeSeasonMANOVA$countSQRT <- sqrt(EarthquakeSeasonMANOVA$count)
plotNormalHistogram(EarthquakeSeasonMANOVA$countSQRT)

# Not much better, let's log it
EarthquakeSeasonMANOVA$countLOG <- log(EarthquakeSeasonMANOVA$count)
plotNormalHistogram(EarthquakeSeasonMANOVA$countLOG)

# We will use the log of EarthquakeSeasonMANOVA$count

# Test Normality of EarthquakeSeasonMANOVA$Magnitude
plotNormalHistogram(EarthquakeSeasonMANOVA$Magnitude)

# EarthquakeSeasonMANOVA$Magnitude looks fairly normal! We will use as is.

# Test Normality of EarthquakeSeasonMANOVA$Median
plotNormalHistogram(EarthquakeSeasonMANOVA$Median)

# EarthquakeSeasonMANOVA$Median also looks normal, so we will use as is

# We will use the Log of count and Magnitude as it is
# Make new matrix
# Keep only two dependent variables, Count and Magnitude
keeps <- c("countLOG", "Magnitude")
EarthquakeSeasonMANOVADVs2 <- EarthquakeSeasonMANOVA[keeps]
View(EarthquakeSeasonMANOVADVs2)

# Format EarthquakeSeasonMANOVADVs2 as a Matrix
EarthquakeSeasonMatrix2 <- as.matrix(EarthquakeSeasonMANOVADVs2)

# Multivariate normality
mshapiro.test(t(EarthquakeSeasonMatrix2))

# ----------------------------------------------------------------------------------------------------------
# Mean of magnitude during each season
EarthquakeSeasonMagnitude <- aggregate(Magnitude ~ Season, NSearthquakes, mean)

#1 Fall      5.875806
#2 Spring    5.883212
#3 Summer    5.886911
#4 Winter    5.871922

# Median of magnitude during each season
EarthquakeSeasonMagnitudeMedian <- aggregate(Magnitude ~ Season, NSearthquakes, median)

#1 Fall                           5.7
#2 Spring                         5.8
#3 Summer                         5.8
#4 Winter                         5.7

# Rename Magnitude Column to Median (to avoid confusion later)
EarthquakeSeasonMagnitudeMedianR <- rename(EarthquakeSeasonMagnitudeMedian, Median = Magnitude)
View(EarthquakeSeasonMagnitudeMedianR)
#----------------------------------------------------------------------
# ChiSquare Test for Magnitude
ObservedMagnitudeSeason <- c(5.88, 5.88, 5.89, 5.87)
ExpectedNullMagnitudeSeason <- c(.25, .25, .25, .25)
chisq.test( x = ObservedMagnitudeSeason, p = ExpectedNullMagnitudeSeason)
# 	Chi-squared test for given probabilities

#data:  ObservedMagnitudeSeason
#X-squared = 3.4014e-05, df = 3, p-value = 1
# p-value is >0.05, so there is no difference than null.
# ------------------------------------------------------------------------------

# Categorize and Quantify points over land or water
library(sp)
library(rgdal)
require(maptools)

# Function from https://stackoverflow.com/questions/23459900/rmetrics-test-if-longitude-and-latitude-coordinates-are-on-land-or-sea
point.in.SpatialPolygons = function(point.x, point.y, SpP){
  X = slot(SpP,"polygons")
  is.inside = F
  for(i in 1:length(X)){
    PS   = slot(X[[i]],"Polygons")
    for(j in 1:length(PS)){
      pol.x = slot(PS[[j]],"coords")[,1]
      pol.y = slot(PS[[j]],"coords")[,2]
      pointsPosition = point.in.polygon(point.x, point.y, pol.x, pol.y)
      if(!slot(PS[[j]],"hole")) {
        is.inside = is.inside | pointsPosition != 0
      }
    }
  }
  is.outsideHole = T
  for(i in 1:length(X)){
    PS   = slot(X[[i]],"Polygons")
    for(j in 1:length(PS)){
      pol.x = slot(PS[[j]],"coords")[,1]
      pol.y = slot(PS[[j]],"coords")[,2]
      pointsPosition = point.in.polygon(point.x, point.y, pol.x, pol.y)
      if(slot(PS[[j]],"hole")) {
        is.outsideHole = is.outsideHole & (pointsPosition == 0 |  pointsPosition == 3)
      }
    }
  }
  is.inside & is.outsideHole
}

# Oceans' Shapefile (from https://www.naturalearthdata.com/downloads/110m-physical-vectors/  Oceans)
setwd("/Users/gonzalezfamily/Desktop/Desktop Entity/DSO110 Final Group Project/Earthquakes/ne_110m_ocean")

oceans <- readShapePoly("ne_110m_ocean.shp", repair=TRUE)
View(oceans)

# New Column to Catorgorize plot point as True for Ocean location and False for land location
earthquakes$LandSea <- point.in.SpatialPolygons(earthquakes$Longitude, earthquakes$Latitude, oceans)

# Count of earthquakes over ocean/true or land/false
EarthquakesOceanLand <- NA
EarthquakesOceanLand$Ocean <- sum(earthquakes$LandSea == TRUE)
EarthquakesOceanLand$Land <- sum(earthquakes$LandSea == FALSE)
View(EarthquakesOceanLand)

# Ocean locations - 18,523
# Land locations - 4,889
# Total locations - 23,412

# One proportion testing in R
# H0 : Ocean locations = Land locations
# Hz : Ocean locations != Land locations

prop.test(x = 18523, n = 23412, p = 0.71, alternative = "two.sided", correct = FALSE)

# 	1-sample proportions test without continuity correction

# data:  18523 out of 23412, null probability 0.5
# X-squared = 7939.8, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.5
# 95 percent confidence interval:
#   0.7859213 0.7963341
# sample estimates:
#   p 
# 0.7911755 

# p-value of 2.2e-16 is < 0.05, so there is a statistical difference
# There is a statistically significant difference between the number of ocean
#    and land earthquake locations, with ocean locations having more earthquakes.

# ------------------------------------------------------------------------------
# Does the depth of the earthquake affect the magnitude of the earthquake?
# Linear Regression

# Uses
# library("car")
# library("caret")
# library("gvlma")
# library("predictmeans")
# library("e1071")
# library("lmtest")

# Test Assumptions
# Test for Linearity

scatter.smooth(x = earthquakes$Depth, y = earthquakes$Magnitude, main = "Earthquake Magnitude by Depth", xlab = "Depth", ylab = "Magnitude")

# Does not pass assumption of linearity, so the depth of earthquakes do not affect magnitude.

# ------------------------------------------------------------------------------
# Separate by land/water and run magnitude/depth for each
EarthquakesInWater <- earthquakes %>% filter(LandSea == TRUE) 
View(EarthquakesInWater)

EarthquakesOnLand <- earthquakes %>% filter(LandSea == FALSE)
View(EarthquakesOnLand)

scatter.smooth(x = EarthquakesInWater$Magnitude, y = EarthquakesInWater$Depth, ylim = rev(range(EarthquakesInWater$Depth)), main = "Ocean Earthquakes Magnitude by Depth", xlab = "Magnitude", ylab = "Depth")

scatter.smooth(x = EarthquakesOnLand$Magnitude, y = EarthquakesOnLand$Depth, ylim = rev(range(EarthquakesOnLand$Depth)), main = "Land Earthquakes Magnitude by Depth", xlab = "Magnitude", ylab = "Depth")

# ---------------------------------------------------------------------
# Does the location on the earth impact the frequency and magnitude of earthquakes?
# How many earthquakes per section above?
EarthquakesPerZone <- earthquakes %>% count(TropicZone)
View(EarthquakesPerZone)

# Remove the 1 occurance of an earthquake right on equator
EarthquakesPerZone1 <- subset(EarthquakesPerZone[2:5,])
View(EarthquakesPerZone1)
# ------------------------------------------------------------------------------
# Import KML data from Google Earth Plate Boundary file from USGS
library(sf)
library(maptools)

temp <- tempfile()
temp2 <- tempfile()
download.file("https://prd-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/atoms/files/plate-boundaries.kmz", temp)
unzip(zipfile = temp, exdir = temp2)
GoogleEarthPlateBoundaries <- st_read(file.path(temp2, "doc.kml"))
unlink(c(temp,temp2))

View(GoogleEarthPlateBoundaries)

# Split Description Column into Multiple Columns to extract data
GoogleEarthPlateBoundaries2 <- separate(GoogleEarthPlateBoundaries, Description, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33"), sep = "<")
View(GoogleEarthPlateBoundaries2)

# Keep only Name, V25, and Geometry
keeps <- c("Name", "V25", "geometry")
GoogleEarthPlateBoundaries3 <- GoogleEarthPlateBoundaries2[keeps]
View(GoogleEarthPlateBoundaries3)

# Remove unnecessary XML info
GoogleEarthPlateBoundaries4 <- separate(GoogleEarthPlateBoundaries3, V25, c("delete", "PlateBoundaryType"), sep = ">")
View(GoogleEarthPlateBoundaries4)

keepsagain <- c("Name", "PlateBoundaryType", "geometry")
GoogleEarthPlateBoundaries5 <- GoogleEarthPlateBoundaries4[keepsagain]
View(GoogleEarthPlateBoundaries5)

#Export as csv
write_csv(GoogleEarthPlateBoundaries5, "GoogleEarthPlateBoundaries5.csv")

# Unnest list in 'geometry'
# Boundary Coordinates into Dataframe
BoundaryCoordinates <- as.data.frame(do.call(rbind,GoogleEarthPlateBoundaries5$geometry))
View(BoundaryCoordinates)

# Four identical columns of vectors, drop all but V1
keeps3 <- c("V1")
BoundaryCoordinates1 <- BoundaryCoordinates[keeps3]
View(BoundaryCoordinates1)

#Remove Zeros
BoundaryCoordinates2 <- cbind(lapply(BoundaryCoordinates1$V1, function(x) {x[x!=0]}))
View(BoundaryCoordinates2)

# Split Lists in Half to separate lats from longs
BoundaryCoordinates3 <- lapply(BoundaryCoordinates2, function(x) list(
  Longitudes = I(list(x[1:(length(x) / 2)])), 
  Latitudes = I(list(x[(length(x) / 2 + 1):length(x)]))))
BoundaryCoordinates3 <- do.call('rbind',BoundaryCoordinates3)   
View(BoundaryCoordinates3)

# Convert to dataframe
BoundaryCoordinates4 <- as.data.frame(BoundaryCoordinates3)
View(BoundaryCoordinates4)

Longitudes <- as.data.frame(do.call(rbind,BoundaryCoordinates4$Longitudes))
View(Longitudes)

Latitudes <- as.data.frame(do.call(rbind,BoundaryCoordinates4$Latitudes))
View(Latitudes)

# Separate out again
Longitudes1 <- unnest_wider(Longitudes, col = V1, names_sep = " Long ")
View(Longitudes1)

Latitudes1 <- unnest_wider(Latitudes, col = V1, names_sep = " Lat ")
View(Latitudes1)

# Rename Column names
names(Longitudes1)[1:98] <- c("Long1","Long2","Long3","Long4","Long5","Long6","Long7","Long8","Long9","Long10","Long11","Long12","Long13","Long14","Long15","Long16","Long17","Long18","Long19","Long20","Long21","Long22","Long23","Long24","Long25","Long26","Long27","Long28","Long29","Long30","Long31","Long32","Long33","Long34","Long35","Long36","Long37","Long38","Long39","Long40","Long41","Long42","Long43","Long44","Long45","Long46","Long47","Long48","Long49","Long50")
names(Latitudes1)[1:98] <- c("Lat1","Lat2","Lat3","Lat4","Lat5","Lat6","Lat7","Lat8","Lat9","Lat10","Lat11","Lat12","Lat13","Lat14","Lat15","Lat16","Lat17","Lat18","Lat19","Lat20","Lat21","Lat22","Lat23","Lat24","Lat25","Lat26","Lat27","Lat28","Lat29","Lat30","Lat31","Lat32","Lat33","Lat34","Lat35","Lat36","Lat37","Lat38","Lat39","Lat40","Lat41","Lat42","Lat43","Lat44","Lat45","Lat46","Lat47","Lat48","Lat49","Lat50")

# Remove V1 and space
names(Longitudes1) <- gsub("V1", "", names(Longitudes1))
names(Longitudes1) <- gsub(" ", "", names(Longitudes1))

names(Latitudes1) <- gsub("V1", "", names(Latitudes1))
names(Latitudes1) <- gsub(" ", "", names(Latitudes1))

# Change BACK to dataframes
Latitudes2 <- as.data.frame(Latitudes1)
Longitudes2 <- as.data.frame(Longitudes1)

# Combine datasets
BoundaryCoordinates5 <- cbind(Longitudes2[1:98],Latitudes2[1:98])
View(BoundaryCoordinates5)

# Combine coordinates to make geopoints
BoundaryCoordinates5$Point1 <- paste(BoundaryCoordinates5$Lat1,',',BoundaryCoordinates5$Long1)
BoundaryCoordinates5$Point2 <- paste(BoundaryCoordinates5$Lat2,',',BoundaryCoordinates5$Long2)
BoundaryCoordinates5$Point3 <- paste(BoundaryCoordinates5$Lat3,',',BoundaryCoordinates5$Long3)
BoundaryCoordinates5$Point4 <- paste(BoundaryCoordinates5$Lat4,',',BoundaryCoordinates5$Long4)
BoundaryCoordinates5$Point5 <- paste(BoundaryCoordinates5$Lat5,',',BoundaryCoordinates5$Long5)
BoundaryCoordinates5$Point6 <- paste(BoundaryCoordinates5$Lat6,',',BoundaryCoordinates5$Long6)
BoundaryCoordinates5$Point7 <- paste(BoundaryCoordinates5$Lat7,',',BoundaryCoordinates5$Long7)
BoundaryCoordinates5$Point8 <- paste(BoundaryCoordinates5$Lat8,',',BoundaryCoordinates5$Long8)
BoundaryCoordinates5$Point9 <- paste(BoundaryCoordinates5$Lat9,',',BoundaryCoordinates5$Long9)
BoundaryCoordinates5$Point10 <- paste(BoundaryCoordinates5$Lat10,',',BoundaryCoordinates5$Long10)
BoundaryCoordinates5$Point11 <- paste(BoundaryCoordinates5$Lat11,',',BoundaryCoordinates5$Long11)
BoundaryCoordinates5$Point12 <- paste(BoundaryCoordinates5$Lat12,',',BoundaryCoordinates5$Long12)
BoundaryCoordinates5$Point13 <- paste(BoundaryCoordinates5$Lat13,',',BoundaryCoordinates5$Long13)
BoundaryCoordinates5$Point14 <- paste(BoundaryCoordinates5$Lat14,',',BoundaryCoordinates5$Long14)
BoundaryCoordinates5$Point15 <- paste(BoundaryCoordinates5$Lat15,',',BoundaryCoordinates5$Long15)
BoundaryCoordinates5$Point16 <- paste(BoundaryCoordinates5$Lat16,',',BoundaryCoordinates5$Long16)
BoundaryCoordinates5$Point17 <- paste(BoundaryCoordinates5$Lat17,',',BoundaryCoordinates5$Long17)
BoundaryCoordinates5$Point18 <- paste(BoundaryCoordinates5$Lat18,',',BoundaryCoordinates5$Long18)
BoundaryCoordinates5$Point19 <- paste(BoundaryCoordinates5$Lat19,',',BoundaryCoordinates5$Long19)
BoundaryCoordinates5$Point20 <- paste(BoundaryCoordinates5$Lat20,',',BoundaryCoordinates5$Long20)
BoundaryCoordinates5$Point21 <- paste(BoundaryCoordinates5$Lat21,',',BoundaryCoordinates5$Long21)
BoundaryCoordinates5$Point22 <- paste(BoundaryCoordinates5$Lat22,',',BoundaryCoordinates5$Long22)
BoundaryCoordinates5$Point23 <- paste(BoundaryCoordinates5$Lat23,',',BoundaryCoordinates5$Long23)
BoundaryCoordinates5$Point24 <- paste(BoundaryCoordinates5$Lat24,',',BoundaryCoordinates5$Long24)
BoundaryCoordinates5$Point25 <- paste(BoundaryCoordinates5$Lat25,',',BoundaryCoordinates5$Long25)
BoundaryCoordinates5$Point26 <- paste(BoundaryCoordinates5$Lat26,',',BoundaryCoordinates5$Long26)
BoundaryCoordinates5$Point27 <- paste(BoundaryCoordinates5$Lat27,',',BoundaryCoordinates5$Long27)
BoundaryCoordinates5$Point28 <- paste(BoundaryCoordinates5$Lat28,',',BoundaryCoordinates5$Long28)
BoundaryCoordinates5$Point29 <- paste(BoundaryCoordinates5$Lat29,',',BoundaryCoordinates5$Long29)
BoundaryCoordinates5$Point30 <- paste(BoundaryCoordinates5$Lat30,',',BoundaryCoordinates5$Long30)
BoundaryCoordinates5$Point31 <- paste(BoundaryCoordinates5$Lat31,',',BoundaryCoordinates5$Long31)
BoundaryCoordinates5$Point32 <- paste(BoundaryCoordinates5$Lat32,',',BoundaryCoordinates5$Long32)
BoundaryCoordinates5$Point33 <- paste(BoundaryCoordinates5$Lat33,',',BoundaryCoordinates5$Long33)
BoundaryCoordinates5$Point34 <- paste(BoundaryCoordinates5$Lat34,',',BoundaryCoordinates5$Long34)
BoundaryCoordinates5$Point35 <- paste(BoundaryCoordinates5$Lat35,',',BoundaryCoordinates5$Long35)
BoundaryCoordinates5$Point36 <- paste(BoundaryCoordinates5$Lat36,',',BoundaryCoordinates5$Long36)
BoundaryCoordinates5$Point37 <- paste(BoundaryCoordinates5$Lat37,',',BoundaryCoordinates5$Long37)
BoundaryCoordinates5$Point38 <- paste(BoundaryCoordinates5$Lat38,',',BoundaryCoordinates5$Long38)
BoundaryCoordinates5$Point39 <- paste(BoundaryCoordinates5$Lat39,',',BoundaryCoordinates5$Long39)
BoundaryCoordinates5$Point40 <- paste(BoundaryCoordinates5$Lat40,',',BoundaryCoordinates5$Long40)
BoundaryCoordinates5$Point41 <- paste(BoundaryCoordinates5$Lat41,',',BoundaryCoordinates5$Long41)
BoundaryCoordinates5$Point42 <- paste(BoundaryCoordinates5$Lat42,',',BoundaryCoordinates5$Long42)
BoundaryCoordinates5$Point43 <- paste(BoundaryCoordinates5$Lat43,',',BoundaryCoordinates5$Long43)
BoundaryCoordinates5$Point44 <- paste(BoundaryCoordinates5$Lat44,',',BoundaryCoordinates5$Long44)
BoundaryCoordinates5$Point45 <- paste(BoundaryCoordinates5$Lat45,',',BoundaryCoordinates5$Long45)
BoundaryCoordinates5$Point46 <- paste(BoundaryCoordinates5$Lat46,',',BoundaryCoordinates5$Long46)
BoundaryCoordinates5$Point47 <- paste(BoundaryCoordinates5$Lat47,',',BoundaryCoordinates5$Long47)
BoundaryCoordinates5$Point48 <- paste(BoundaryCoordinates5$Lat48,',',BoundaryCoordinates5$Long48)
BoundaryCoordinates5$Point49 <- paste(BoundaryCoordinates5$Lat49,',',BoundaryCoordinates5$Long49)
BoundaryCoordinates5$Point50 <- paste(BoundaryCoordinates5$Lat50,',',BoundaryCoordinates5$Long50)
BoundaryCoordinates5$Point51 <- paste(BoundaryCoordinates5$Lat51,',',BoundaryCoordinates5$Long51)
BoundaryCoordinates5$Point52 <- paste(BoundaryCoordinates5$Lat52,',',BoundaryCoordinates5$Long52)
BoundaryCoordinates5$Point53 <- paste(BoundaryCoordinates5$Lat53,',',BoundaryCoordinates5$Long53)
BoundaryCoordinates5$Point54 <- paste(BoundaryCoordinates5$Lat54,',',BoundaryCoordinates5$Long54)
BoundaryCoordinates5$Point55 <- paste(BoundaryCoordinates5$Lat55,',',BoundaryCoordinates5$Long55)
BoundaryCoordinates5$Point56 <- paste(BoundaryCoordinates5$Lat56,',',BoundaryCoordinates5$Long56)
BoundaryCoordinates5$Point57 <- paste(BoundaryCoordinates5$Lat57,',',BoundaryCoordinates5$Long57)
BoundaryCoordinates5$Point58 <- paste(BoundaryCoordinates5$Lat58,',',BoundaryCoordinates5$Long58)
BoundaryCoordinates5$Point59 <- paste(BoundaryCoordinates5$Lat59,',',BoundaryCoordinates5$Long59)
BoundaryCoordinates5$Point60 <- paste(BoundaryCoordinates5$Lat60,',',BoundaryCoordinates5$Long60)
BoundaryCoordinates5$Point61 <- paste(BoundaryCoordinates5$Lat61,',',BoundaryCoordinates5$Long61)
BoundaryCoordinates5$Point62 <- paste(BoundaryCoordinates5$Lat62,',',BoundaryCoordinates5$Long62)
BoundaryCoordinates5$Point63 <- paste(BoundaryCoordinates5$Lat63,',',BoundaryCoordinates5$Long63)
BoundaryCoordinates5$Point64 <- paste(BoundaryCoordinates5$Lat64,',',BoundaryCoordinates5$Long64)
BoundaryCoordinates5$Point65 <- paste(BoundaryCoordinates5$Lat65,',',BoundaryCoordinates5$Long65)
BoundaryCoordinates5$Point66 <- paste(BoundaryCoordinates5$Lat66,',',BoundaryCoordinates5$Long66)
BoundaryCoordinates5$Point67 <- paste(BoundaryCoordinates5$Lat67,',',BoundaryCoordinates5$Long67)
BoundaryCoordinates5$Point68 <- paste(BoundaryCoordinates5$Lat68,',',BoundaryCoordinates5$Long68)
BoundaryCoordinates5$Point69 <- paste(BoundaryCoordinates5$Lat69,',',BoundaryCoordinates5$Long69)
BoundaryCoordinates5$Point70 <- paste(BoundaryCoordinates5$Lat70,',',BoundaryCoordinates5$Long70)
BoundaryCoordinates5$Point71 <- paste(BoundaryCoordinates5$Lat71,',',BoundaryCoordinates5$Long71)
BoundaryCoordinates5$Point72 <- paste(BoundaryCoordinates5$Lat72,',',BoundaryCoordinates5$Long72)
BoundaryCoordinates5$Point73 <- paste(BoundaryCoordinates5$Lat73,',',BoundaryCoordinates5$Long73)
BoundaryCoordinates5$Point74 <- paste(BoundaryCoordinates5$Lat74,',',BoundaryCoordinates5$Long74)
BoundaryCoordinates5$Point75 <- paste(BoundaryCoordinates5$Lat75,',',BoundaryCoordinates5$Long75)
BoundaryCoordinates5$Point76 <- paste(BoundaryCoordinates5$Lat76,',',BoundaryCoordinates5$Long76)
BoundaryCoordinates5$Point77 <- paste(BoundaryCoordinates5$Lat77,',',BoundaryCoordinates5$Long77)
BoundaryCoordinates5$Point78 <- paste(BoundaryCoordinates5$Lat78,',',BoundaryCoordinates5$Long78)
BoundaryCoordinates5$Point79 <- paste(BoundaryCoordinates5$Lat79,',',BoundaryCoordinates5$Long79)
BoundaryCoordinates5$Point80 <- paste(BoundaryCoordinates5$Lat80,',',BoundaryCoordinates5$Long80)
BoundaryCoordinates5$Point81 <- paste(BoundaryCoordinates5$Lat81,',',BoundaryCoordinates5$Long81)
BoundaryCoordinates5$Point82 <- paste(BoundaryCoordinates5$Lat82,',',BoundaryCoordinates5$Long82)
BoundaryCoordinates5$Point83 <- paste(BoundaryCoordinates5$Lat83,',',BoundaryCoordinates5$Long83)
BoundaryCoordinates5$Point84 <- paste(BoundaryCoordinates5$Lat84,',',BoundaryCoordinates5$Long84)
BoundaryCoordinates5$Point85 <- paste(BoundaryCoordinates5$Lat85,',',BoundaryCoordinates5$Long85)
BoundaryCoordinates5$Point86 <- paste(BoundaryCoordinates5$Lat86,',',BoundaryCoordinates5$Long86)
BoundaryCoordinates5$Point87 <- paste(BoundaryCoordinates5$Lat87,',',BoundaryCoordinates5$Long87)
BoundaryCoordinates5$Point88 <- paste(BoundaryCoordinates5$Lat88,',',BoundaryCoordinates5$Long88)
BoundaryCoordinates5$Point89 <- paste(BoundaryCoordinates5$Lat89,',',BoundaryCoordinates5$Long89)
BoundaryCoordinates5$Point90 <- paste(BoundaryCoordinates5$Lat90,',',BoundaryCoordinates5$Long90)
BoundaryCoordinates5$Point91 <- paste(BoundaryCoordinates5$Lat91,',',BoundaryCoordinates5$Long91)
BoundaryCoordinates5$Point92 <- paste(BoundaryCoordinates5$Lat92,',',BoundaryCoordinates5$Long92)
BoundaryCoordinates5$Point93 <- paste(BoundaryCoordinates5$Lat93,',',BoundaryCoordinates5$Long93)
BoundaryCoordinates5$Point94 <- paste(BoundaryCoordinates5$Lat94,',',BoundaryCoordinates5$Long94)
BoundaryCoordinates5$Point95 <- paste(BoundaryCoordinates5$Lat95,',',BoundaryCoordinates5$Long95)
BoundaryCoordinates5$Point96 <- paste(BoundaryCoordinates5$Lat96,',',BoundaryCoordinates5$Long96)
BoundaryCoordinates5$Point97 <- paste(BoundaryCoordinates5$Lat97,',',BoundaryCoordinates5$Long97)
BoundaryCoordinates5$Point98 <- paste(BoundaryCoordinates5$Lat98,',',BoundaryCoordinates5$Long98)


BoundaryCoordinates7 <- BoundaryCoordinates5[197:294]
View(BoundaryCoordinates7)

# Rename columns
names(Longitudes) <- gsub("V1", "Longitudes", names(Longitudes))
names(Latitudes) <- gsub("V1", "Latitudes", names(Latitudes))

# Add columns back in to original
GoogleEarthPlateBoundaries6 <- bind_cols(GoogleEarthPlateBoundaries5, BoundaryCoordinates7)
View(GoogleEarthPlateBoundaries6)

# Drop geometry lists
drop <- c("geometry")
GoogleEarthPlateBoundaries7 <- GoogleEarthPlateBoundaries6[,!(names(GoogleEarthPlateBoundaries6) %in% drop)]

#Export to Excel
write_xlsx(GoogleEarthPlateBoundaries6, "GoogleEarthPlateBoundaries6.xlsx")
write_csv(GoogleEarthPlateBoundaries6, "GoogleEarthPlateBoundaries6.csv")

# ------------------------------------------------------------------------------

# Add Distance from Distance tool at "http://www.zonums.com/online/kmlArea/"
library(readxl)
BoundaryDistance <- read_excel("~/Desktop/Desktop Entity/DSO110 Final Group Project/Earthquakes/GoogleEarthPlateBoundaries7.xlsx")
View(BoundaryDistance)

# Distances of Plate Barriers
PlateBoundaryMeasurements <- BoundaryDistance %>% group_by(Name) %>%
  summarize(SumDist = sum(DistanceKM))
View(PlateBoundaryMeasurements)

# Distances by Section of Earth
DistancebySection <- BoundaryDistance %>% group_by(Section) %>%
  summarize(Distance = sum(DistanceKM))
View(DistancebySection)

#1 Northern              55869.
#2 Southern              59932.
#3 Tropic of Cancer      35934.
#4 Tropic of Capricorn   37873.

#------------------------------------------------------------------------------
# Combine EarthquakesPerZone1 and DistancebySection
EarthquakesPerZone2 <- merge(x = EarthquakesPerZone1, y = DistancebySection, by.x = "TropicZone", by.y = "Section", all = TRUE)
View(EarthquakesPerZone2)

# Rename Column
names(EarthquakesPerZone2)[names(EarthquakesPerZone2) == "n"] <- "NumberOfEarthquakes"
names(EarthquakesPerZone2)[names(EarthquakesPerZone2) == "Distance"] <- "PlateBoundaryDistance"

print(EarthquakesPerZone2)
write_csv(EarthquakesPerZone2, "EarthquakesPerZone.csv")

# Distance by Plate Type
#      TropicZone             NumberOfEarthquakes    PlateBoundaryDistance
# 1            Northern                6248              55869.41
# 2            Southern                4121              59931.75
# 3    Tropic of Cancer                4369              35933.56
# 4 Tropic of Capricorn                8673              37873.44

# Turn Plate Boundary Distance into percentage of whole earth to run a Goodness of Fit Chi-Square
EarthquakesPerZone2$PercentofWhole <- EarthquakesPerZone2$PlateBoundaryDistance/sum(EarthquakesPerZone2$PlateBoundaryDistance)

# Set up Observed Values of Number of Earthquakes
observed = c(6248, 4121, 4369, 8673)
expected = c(0.29, 0.32, 0.19, 0.20)

# ChiSquare
chisq.test( x = observed, p = expected)

#Chi-squared test for given probabilities

#data:  observed
#X-squared = 4962.5, df = 3, p-value < 2.2e-16

#p-value is < -.05, so rates of earthquakes do not line up with the total amount of distance of plate boundaries


# --------------------------------------------------------------
# depth v. magnitude

# Tectonic Summaries for M7+ Earthquakes 2000-2015

temp <- tempfile()
temp2 <- tempfile()
download.file("https://prd-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/atoms/files/tectonic-summ2000-2015.kmz", temp)
unzip(zipfile = temp, exdir = temp2)
TectonicSummary <- st_read(file.path(temp2, "doc.kml"))
unlink(c(temp,temp2))

View(TectonicSummary)

# Wrangle "Description"
# Split Description Column into Multiple Columns to extract data
TectonicSummary2 <- separate(TectonicSummary, Description, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33"), sep = "<")
View(TectonicSummary2)

TectonicSummary3 = subset(TectonicSummary2, select = -c(V1, V5, V7, V23, V25, V27, V28, V29, V30, V31, V32, V33))
View(TectonicSummary3)

TectonicSummary4 <- separate(TectonicSummary3, Name, c("V1", "Name"), sep = ":")
View(TectonicSummary4)

TectonicSummary5 <- separate(TectonicSummary4, V2, c("Delete", "Date"), sep = ">")
View(TectonicSummary5)

TectonicSummary6 <- separate(TectonicSummary5, Date, c("Date", "Time", "Delete2", "Delete3"), sep = " ")
View(TectonicSummary6)

TectonicSummary7 <- separate(TectonicSummary6, V3, c("D4", "Coordinates"), sep = ">")
View(TectonicSummary7)

TectonicSummary8 <- separate(TectonicSummary7, V4, c("D5", "Magnitude"), sep = ">")
View(TectonicSummary8)

TectonicSummary9 <- separate(TectonicSummary8, V6, c("D6", "Depth"), sep = ">")
View(TectonicSummary9)

TectonicSummary10 <- separate(TectonicSummary9, V8, c("D7", "Description"), sep = ">")
View(TectonicSummary10)

TectonicSummary11 = subset(TectonicSummary10, select = -c(V1, Delete, Delete2, Delete3, D4, D5, D6, D7))
View(TectonicSummary11)

str(TectonicSummary11$Magnitude)

TectonicSummary12 <- strsplit(TectonicSummary11$Magnitude, split = " - ")
TectonicSummary12 <- separate(TectonicSummary11, Magnitude, c("Magnitude", "Region"), sep = "-")
View(TectonicSummary12)

keeps <- c("Name", "V25", "geometry")
GoogleEarthPlateBoundaries3 <- GoogleEarthPlateBoundaries2[keeps]
View(GoogleEarthPlateBoundaries3)

GoogleEarthPlateBoundaries4 <- separate(GoogleEarthPlateBoundaries3, V25, c("delete", "PlateBoundaryType"), sep = ">")
View(GoogleEarthPlateBoundaries4)

keepsagain <- c("Name", "PlateBoundaryType", "geometry")
GoogleEarthPlateBoundaries5 <- GoogleEarthPlateBoundaries4[keepsagain]
View(GoogleEarthPlateBoundaries5)
# ---------------------------------------------------------------------
# Has the number of earthquakes per year changed over time?
earthquakesByYear <- separate(earthquakes, Date, c("Year", "Month", "Day"))
View(earthquakesByYear)

earthquakeCountYears <- earthquakesByYear %>% group_by(Year) %>% summarize(count = n())
View(earthquakeCountYears)

# drop NAs
earthquakeCountYears <- drop_na(earthquakeCountYears)

ggplot(earthquakeCountYears, aes( x = Year, y = count)) +
  ggtitle("Numbers of Earthquakes over the Years") +
  scale_y_continuous(name = "Number of Earthquakes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(hjust = 0.5)) +
  geom_col()

# Linear Model
EarthquakesOverTheYears <- lm(count ~ Year, data = earthquakeCountYears)
par(mfrow=c(2,2))
plot(EarthquakesOverTheYears)


# Breush-Pagan test
lmtest::bptest(EarthquakesOverTheYears)
#studentized Breusch-Pagan test
#data:  EarthquakesOverTheYears
#BP = 3.9961, df = 1, p-value = 0.0456

# literally right at cutoff for homoscedasticity

# NCV test
car::ncvTest(EarthquakesOverTheYears)
#Non-constant Variance Score Test 
#Variance formula: ~ fitted.values 
#Chisquare = 3.236142, Df = 1, p = 0.07203
# p-value > 0.05, so this barely passes test for homoscedasticity

summary(EarthquakesOverTheYears)
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -8033.9657  1192.9835  -6.734 1.56e-08 ***
#  Year            4.2623     0.5993   7.112 4.01e-09 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 64.86 on 50 degrees of freedom
#Multiple R-squared:  0.5029,	Adjusted R-squared:  0.4929 
#F-statistic: 50.58 on 1 and 50 DF,  p-value: 4.01e-09

# If the trend continues, we can expect to see about 4 more earthquakes every year.

# Plot Trend line

plot(earthquakeCountYears$Year, earthquakeCountYears$count)
lines(earthquakeCountYears$Year, predict(EarthquakesOverTheYears), col = 2, lwd = 2)

earthquakeCountYears$Year <- as.numeric(earthquakeCountYears$Year)
ggplot(earthquakeCountYears, aes(x = Year, y = count)) +
  ggtitle("Earthquakes Over the Years") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_y_continuous("Number of Earthquakes") +
  theme(plot.title = element_text(hjust = 0.5))

#------------------------------------------------------------------------------------
# Is the increase in earthquakes shown due to increase number of sensors/monitoring?

#Download station txt file (http://www.iris.washington.edu/gmap/#network=_REALTIME&plates=on&planet=earth)
library(readr)

StationXML <- read_delim("~/Desktop/Desktop Entity/DSO110 Final Group Project/Earthquakes/gmap-stations.csv", 
                            "|", escape_double = FALSE, trim_ws = TRUE)
View(StationXML)

#Export to Excel to wrangle DataCenter Info
install.packages('xlsx')
library('xlsx')
write.xlsx(StationXML, "StationXML.xlsx")

# Import back in after fixing DataCenter info and also duplicate station Dates occurring in the future as checked at http://www.iris.washington.edu/gmap/#network=_REALTIME&station=MWC&plates=on&planet=earth
StationXML1 <- read_excel("~/Desktop/Desktop Entity/DSO110 Final Group Project/Earthquakes/StationXML.xlsx")
View(StationXML1) 
str(StationXML1)

# Separate Start Date Year from Date and Time
StationXML2 <- separate(StationXML1, StartTime, c("Date", "Time"), sep = " ")
View(StationXML2)
StationXML3 <- separate(StationXML2, Date, c("Year Deployed", "Month", "Day"), sep = "-")
View(StationXML3)

#Change Column Name
StationXML3$YearDeployed <- StationXML3$`Year Deployed`
StationYears <- StationXML3 %>% group_by(YearDeployed) %>% summarize(count = n())
View(StationYears)

#Export StationYears and earthquakeCountYears to Excel
write.xlsx(StationYears, "StationYears.xlsx")
write.xlsx(earthquakeCountYears, "earthquakeCountYears.xlsx")

#After combining info read back earthquakeCountYears with Station number and Increase for each year of data
earthquakeCountYearsStations <- read_excel("~/Desktop/Desktop Entity/DSO110 Final Group Project/Earthquakes/earthquakeCountYears.xlsx")
View(earthquakeCountYearsStations)

# Plot Stations v. Earthquake Count
scatter.smooth(x = earthquakeCountYearsStations$count, y = earthquakeCountYearsStations$Stations)

ggplot(earthquakeCountYearsStations) +
  ggtitle("Comparison of Number of Recorded Earthquakes and Number of Active Seismometers per Year") +
  geom_bar(aes( x = Year, y = count), stat = "identity") +
  geom_line(aes(x = Year, y = Stations), stat = "identity") +
  scale_y_continuous(name = "Number of Earthquakes Each Year", 
                     sec.axis = sec_axis(~., name = "Number of Seismometers") )


# ------------------------------------------------------------------------------
# Are certain types of plate boundaries more likely to have earthquakes?
# Import joined data from QGIS - PlateBoundaries and earthqakes
# Closest PlateBoundary was joined to each earthquake data.

EarthquakesToPlateInterface <- read.csv("~/Desktop/Desktop Entity/DSO110 Final Group Project/Earthquakes/EarthquakesToPlateInterface.csv")
View(EarthquakesToPlateInterface)
str(EarthquakesToPlateInterface)

#Edit out plate boundary type from "join_description"
# Separate out plate boundary by "<"
EarthquakesToPlateInterface2 <- separate(EarthquakesToPlateInterface, join_description, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33"), sep = "<")
View(EarthquakesToPlateInterface2)

#Remove out data I don't need
EarthquakesToPlateInterface3 = subset(EarthquakesToPlateInterface2, select = -c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21, V22, V23, V24, V26, V27, V28, V29, V30, V31, V32, V33))
View(EarthquakesToPlateInterface3)

#Separate out V25 by ">"
EarthquakesToPlateInterface4 <- separate(EarthquakesToPlateInterface3, V25, c("delete", "PlateBoundaryType"), sep = ">")
View(EarthquakesToPlateInterface4)
str(EarthquakesToPlateInterface4)

# Group by PlateBoundaryType
EarthquakesOnEachPlateBoundaryType <- EarthquakesToPlateInterface4 %>% group_by(PlateBoundaryType) %>% summarize(count = n())
View(EarthquakesOnEachPlateBoundaryType)

names(EarthquakesOnEachPlateBoundaryType)[names(EarthquakesOnEachPlateBoundaryType) == "count"] <- "NumberOfEarthquakes"

# Compare to distance of boundaries of each type
BoundaryDistanceByType <- BoundaryDistance %>% group_by(PlateBoundaryType) %>% summarize(TotalPlateBoundaryDistanceKM = sum(DistanceKM))
View(BoundaryDistanceByType)

PlateBoundaryComparison <- merge(EarthquakesOnEachPlateBoundaryType, BoundaryDistanceByType, by = "PlateBoundaryType")
View(PlateBoundaryComparison)
write_csv(PlateBoundaryComparison, "PlateBoundaryComparison.csv")

# PlateBoundaryType NumberOfEarthquakes  SumDist
#1 Convergent Boundary               17356 62698.63
#2  Divergent Boundary                1461 61859.67
#3               Other                 981 13163.88
#4  Transform Boundary                3614 51885.98

# Turn Plate Boundary Distance into percentage of whole earth to run a Goodness of Fit Chi-Square
PlateBoundaryComparison$PercentofWhole <- PlateBoundaryComparison$TotalPlateBoundaryDistanceKM/sum(PlateBoundaryComparison$TotalPlateBoundaryDistanceKM)

# Set up Observed Values of Number of Earthquakes
observedPlateType = c(17356, 1461, 981, 3614)
expectedPlateType = c(0.33, 0.33, 0.07, 0.27)

# ChiSquare
chisq.test( x = observedPlateType, p = expectedPlateType)

#Chi-squared test for given probabilities

#data:  observedPlateType
#X-squared = 18507, df = 3, p-value < 2.2e-16


#p-value is < -.05, so rates of earthquakes on each boundary type do not line up with what is expected from distance of each plate boundary type

# Look at Percentage of earthquakes on each type of boundary
PlateBoundaryComparison$earthquakePercentage <- PlateBoundaryComparison$NumberOfEarthquakes/sum(PlateBoundaryComparison$NumberOfEarthquakes)
print(PlateBoundaryComparison)
