# DSO110 Final Project
# The Science Explorers
# Borly Green, Jessenia Lorenzo, and Melissa Gonzalez

#How do Earthquakes impact Tsunamis afterward? What about Location, Depth, Season, Hemisphere and other factors?

# IV (categorical) - Date, Time, Location Name, Country, Region, 
# DV (continuous) - Latitude and Longitude, Tsunami Intensity

#Import Packages
library('tidyverse') #clean and tidying data
library('data.table') #aggregation of dataset (remove, add columns, etc.)
library('dplyr') #data manipulation, wrangling
library('IDPmisc') #manipulation (ex. omitting irregular values/selecting data by logical vectors, including NAs)
library("tidyr")
library("mvnormtest") #for MANOVA
library("car") # for test of homogeneity of variance
library('rcompanion') # for plotting normality of histograms
library("caret") # linear regressions
library("gvlma")
library("predictmeans")
library("e1071")
library("lmtest")
library("readxl")
library('geosphere')
library('rgeos')
library('tmap')

#Import the data
tsunami <- read.csv("~/R Files/tsunami")
View(tsunami)

##### WRANGLE TIME!!!
#Subset to what we need
Tsunami <- tsunami[c("ID", "YEAR", "MONTH", "DAY", "HOUR", "MINUTE", "LATITUDE", "LONGITUDE", "LOCATION_NAME", "COUNTRY", "REGION", "CAUSE", "TS_INTENSITY")]
View(Tsunami)

# Remove Missing data
na.omit(Tsunami)
View(Tsunami)

#Export Data to visualize and story tell in Tableau
write.csv(Tsunami, "~/R Files/Tsunami.csv")