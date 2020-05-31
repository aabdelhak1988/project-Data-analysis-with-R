#read data to g
library(haven)
g <- read_sav("Data/Endfile.sav")
View(g)

#change mslauf to ms-Verlauf with only PPMS(1) and SPMS(2)
g$ms_verlauf <- length(g$mslauf) 
g$ms_verlauf[g$mslauf == 4] <- "PPMS"
g$ms_verlauf[g$mslauf == 3] <- "SPMS"

g$ms_verlauf <- as.factor(g$ms_verlauf) # change to factor

#activate psych to perform summary analysis 
library(psych)
basic_variables_cont. <- c("edss", "age", "GFAP_singlePlex", "NfL_singlePlex", "Disease_duration", "Progression_Duration", "ms_verlauf")
basic_variables_data <- g[basic_variables_cont.]
describeBy(basic_variables_data, basic_variables_data$ms_verlauf)

#histograms 
#change width through binwidth = 
# dached line = mean, can be changed to median 
# can change the black lining through color = "desired color"
# combine = false -> law 3ndi akter mi parameter fi x= c("age","edss") -> combine/ merge can be changed to TRUE
library(ggpubr)
gghistogram(g, x = "age", y= "..count..", combine = FALSE, ylab = "Count", add = "mean", fill= "dark grey", rug = TRUE, bins = 25, binwidth = 4)
gghistogram(g, x = "edss", y= "..count..", combine = FALSE, ylab = "Count", add = "mean", fill= "dark grey", rug = TRUE, bins = 25, binwidth = 0.5)
gghistogram(g, x = "GFAP_singlePlex", y= "..count..", combine = FALSE, ylab = "Count", add = "mean", fill= "dark grey", rug = TRUE, bins = 25, binwidth = 15)
gghistogram(g, x = "NfL_singlePlex", y= "..count..", combine = FALSE, ylab = "Count", add = "mean", fill= "dark grey", rug = TRUE, bins = 25, binwidth = 2)


