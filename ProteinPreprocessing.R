# Script created by Frank D'Agostino 2020

rm(list = ls())
# Clear console
# install.packages("Interpol")
library(Interpol)

#-------------------------------------------------------------------------------
# Sequence Preprocessing
#-------------------------------------------------------------------------------

# We must clean the CovAbDab data from Oxford
covAbDab <- read.csv("CoV-AbDab_120720.csv"); head(covAbDab)

temp <- subset(covAbDab, CDRH3 != "ND"); head(temp)

datTemp <- temp$CDRH3; head(datTemp)
dat <- as.vector(datTemp); head(dat)

# Seems like it works!
test <- AAdescriptor(dat[3]); test

desiredLength <- 35

newDat <- matrix(nrow=length(dat), 
                 ncol=desiredLength)

for (i in 1:length(dat)) {
  # AAdescriptor encodes sequences
  # Interpol normalizes to desiredLength
  
  newDat[i,] <- Interpol(AAdescriptor(dat[i]), dims = desiredLength)
}

# We now have a list of normalized lengths
# of amino acid sequences that are numerical as well
print(newDat)

##################################################################

