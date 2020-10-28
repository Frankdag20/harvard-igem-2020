# Script created by Frank D'Agostino 2020

rm(list = ls())
# Clear console
install.packages("Interpol")
library(Interpol)
library(tidyverse)
data(AAindex)

#-------------------------------------------------------------------------------
# Sequence Preprocessing
#-------------------------------------------------------------------------------

# We must clean the CovAbDab data from Oxford
covAbDab <- read.csv("CoV-AbDab_120720.csv"); head(covAbDab)
temp <- subset(covAbDab, CDRH3 != "ND"); head(temp)

# Clean the data to change bind, neutralizing data into numerical data
bindDat <- numeric(nrow(temp))
k <- 0
for (i in temp$Binds.to) {
  if (grepl("SARS-CoV2 (weak)", i, fixed=TRUE)) {
    bindDat[k] <- 0.5
  } else if (grepl("SARS-CoV2", i, fixed=TRUE)) {
    bindDat[k] <- 1
  } else {
    bindDat[k] <- 0
  }
  
  k = k + 1
}

# Repeat for neutralization
neutrDat <- numeric(nrow(temp))
k <- 0
for (i in temp$Neutralising.Vs) {
  if (grepl("SARS-CoV2 (weak)", i, fixed=TRUE)) {
    neutrDat[k] <- 0.5
  } else if (grepl("SARS-CoV2", i, fixed=TRUE)) {
    neutrDat[k] <- 1
  } else {
    neutrDat[k] <- 0
  }
  
  k = k + 1
}

datTemp <- temp$CDRH3; head(datTemp)
dat <- as.vector(datTemp); head(dat)

# Seems like it works!
test <- AAdescriptor(dat[3], descriptor=151); test

desiredLength <- 15

newDat <- matrix(nrow=length(dat), 
                 ncol=desiredLength)

for (i in 1:length(dat)) {
  # AAdescriptor encodes sequences
  # Interpol normalizes to desiredLength
  
  newDat[i,] <- Interpol(AAdescriptor(dat[i]), dims = desiredLength)
}

# Create a function to pass into Python Code
interp <- function(desiredLength, data) {
  # dat <- as.vector(data)
  newDat <- matrix(nrow=length(dat), 
                   ncol=desiredLength)
  
  for (i in 1:length(dat)) {
    # AAdescriptor encodes sequences
    # Interpol normalizes to desiredLength
    
    newDat[i,] <- Interpol(AAdescriptor(dat[i]), dims = desiredLength)
  }
  
  # Extract rows with NA values
  X <- which(is.na(newDat), arr.ind = TRUE)
  Y <- X[1:nrow(X), 1]
  # Remove repeats and remove those in corresponding lists
  Z <- unique(Y)
  bindDatRem <- bindDat[-Z]
  neutrDatRem <- neutrDat[-Z]
  frameDat <- na.omit(newDat)
  writeOutput <- data.frame(frameDat, bindingData = bindDatRem, neutralizationData = neutrDatRem)
  
  return(writeOuput)
}

# Extract rows with NA values
X <- which(is.na(newDat), arr.ind = TRUE)
Y <- X[1:nrow(X), 1]

# Remove repeats and remove those in corresponding lists
Z <- unique(Y)
bindDatRem <- bindDat[-Z]
neutrDatRem <- neutrDat[-Z]

# We now have a list of normalized lengths
# of amino acid sequences that are numerical as well
frameDat <- na.omit(newDat)
#bindDat <- c(bindDat[1:715], bindDat[720:length(bindDat)])
#neutrDat <- c(neutrDat[1:715], neutrDat[720:length(neutrDat)])
print(na.omit(newDat))
newDat[1,]

writeOutput <- data.frame(frameDat, bindingData = bindDatRem, neutralizationData = neutrDatRem)
head(writeOutput)
write.csv(writeOutput, "interpolOutput.csv", row.names=FALSE)

##################################################################

