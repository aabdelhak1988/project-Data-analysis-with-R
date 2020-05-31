#read data to g
library(haven)
g <- read_sav("Data/Endfile.sav")

#change mslauf to ms-Verlauf with only PPMS(1) and SPMS(2)
g$ms_verlauf <- length(g$mslauf) 
g$ms_verlauf[g$mslauf == 4] <- "PPMS"
g$ms_verlauf[g$mslauf == 3] <- "SPMS"

g$ms_verlauf <- as.factor(g$ms_verlauf) # change to factor

#remove labels 
library(labelled)
remove_labels(g$ongoing)
remove_labels(g$sex)

#calculate ratio
g$GFAP_NfL <- g$GFAP_singlePlex/g$NfL_singlePlex

#find out if normal 
shapiro.test(g$edss)
shapiro.test(g$age)
shapiro.test(g$GFAP_singlePlex)
shapiro.test(g$NfL_singlePlex)
shapiro.test(g$Disease_duration)
shapiro.test(g$GFAP_NfL)

#creat log values
g$edss_log <- log(g$edss)
g$gfap_log <- log(g$GFAP_singlePlex)
g$Nfl_log <- log(g$NfL_singlePlex)
g$dis_duration_log <- log(g$Disease_duration)
g$GFAP_NfL_log <- log(g$GFAP_NfL)

#activate psych to perform summary analysis 
library(psych)
basic_variables_cont <- c("edss", "age", "GFAP_singlePlex", "NfL_singlePlex", "GFAP_NfL", "Disease_duration", "Progression_Duration", "ms_verlauf")
basic_variables_data <- g[basic_variables_cont]
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
gghistogram(g, x = "GFAP_NfL", y= "..count..", combine = FALSE, ylab = "Count", add = "mean", fill= "dark grey", rug = TRUE, bins = 25, binwidth = 4)

library(xtable)
library(Hmisc)
library(htmlTable)
library(lattice)
library(Formula)
library(survival)

#add corstars
# x is a matrix containing the data
# method : correlation method. "pearson"" or "spearman"" is supported
# removeTriangle : remove upper or lower triangle
# results :  if "html" or "latex"
# the results will be displayed in html or latex format
corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex")){
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
  
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex") 
  }
} 


#correlation matrix with p values in ***
corr_table_edss_bio <- corstars(g[,c("edss","age","GFAP_singlePlex", "NfL_singlePlex", "GFAP_NfL","Disease_duration")], method = "spearman", result = "html")
htmlTable(corr_table_edss_bio)
my_cols <- c("red", "black") 
pairs(g[,c("age","edss", "GFAP_singlePlex", "Nfl_log", "GFAP_NfL", "Disease_duration")], pch = 19, lower.panel = NULL, col = my_cols[g$ms_verlauf], main = "Correlation matrix, red -> PPMS, black -> SPMS")


#calculate linear regression EDSS GFAP/NFL/Ratio
lm_edss_gfap <- lm(g$edss_log ~ g$gfap_log + g$age + g$dis_duration_log + factor(g$sex) + factor(g$ongoing) + factor(g$ms_verlauf))
summary(lm_edss_gfap)
lm_edss_nfl <- lm(g$edss_log ~ g$Nfl_log + g$age + g$dis_duration_log + factor(g$sex) + factor(g$ongoing) + factor(g$ms_verlauf))
summary(lm_edss_nfl)
lm_edss_GFAP_nfl <- lm(g$edss_log ~ g$GFAP_NfL_log + g$age + g$dis_duration_log + factor(g$sex) + factor(g$ongoing) + factor(g$ms_verlauf))
summary(lm_edss_GFAP_nfl)
str

#comparison_Progression vs activity 
