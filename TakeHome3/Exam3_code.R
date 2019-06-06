####################################################################################################################################
#################################################### Preliminary Data Gathering ####################################################
####################################################################################################################################
# Read In Data #
library(RCurl)
gitstring = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome1/BirdDat.csv"
RTH_Dat <-read.csv(text=getURL(gitstring))[,-4]
gitstr2 = "https://raw.githubusercontent.com/msilva00/AMS207/master/TakeHome3/RedShouldered.csv"
RSH_dat = na.omit(read.csv(text=getURL(gitstr2)))

####################################################################################################################################
head(RTH_Dat)
yi1 = RTH_Dat$RedtailedHawk
i1_years = length(yi1)
sum_x = sum(yi1)
c_i2 = RTH_Dat$RouteCount

yi2 = RSH_dat$RedShouldered
i2_years = RSH_dat$Year
c_i2 = RSH_dat$RouteCount
