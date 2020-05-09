## Create known Ancestral Assignment data from SAEGUS output for use in calculating Ancestral Assignment Accuracy

# Required Packages
library("tidyr")
library("gdata")
library("data.table")
library("dplyr")

# User inputs
# Working directory
wd <- "~/working/directory"
# Name of simulated data
sd <- "simdata_n1000_test_set.csv"
# Name of founder key data
fd <- "test_founder_data_key.txt"
# Number of samples
sn <- 1000

# Set working directory, which contains the simulated output from SAEGUS and the founder key data
setwd(wd)

# SimuPop Output
sim <- read.csv(sd,head=T, stringsAsFactors = FALSE) #read in simulated data
sim <- as.data.frame(t(sim[,3:ncol(sim)])) #Transpose so markers are rows and samples are columns

# Founder Data
key <- read.table(fd, head=T, stringsAsFactors = FALSE, sep="\t")
key <- rename.vars(key,"chr","CHROM") #rename chromosome column
founders <- names(key[,6:ncol(key)]) #extract founder names for later

# Split simdata for formatting, each file corresponds to each allele
sim2_allele1 <- sim[seq(1,nrow(sim),2),]
sim2_allele2 <- sim[seq(2,nrow(sim),2),]

# Replace index with name of founder
# 0=ParentA, 1=ParentB, 2=ParentC, 3=ParentD, 4=ParentE, 5=ParentF, 6=ParentG, key for test data
for (i in 1:length(founders)){
  sim2_allele1 <- apply(sim2_allele1,2,function(x) ifelse(x==(i-1),founders[i],x))
}

for (i in 1:length(founders)){
  sim2_allele2 <- apply(sim2_allele2,2,function(x) ifelse(x==(i-1),founders[i],x))
}

sim2_allele1 <- as.data.frame(sim2_allele1)
sim2_allele2 <- as.data.frame(sim2_allele2)

## Concatenate alleles to represent genotype in format founder/founder...
all_data <- data.frame(snpID=key$snpID)
for (i in 1:ncol(sim2_allele1)){
  all_data[,i+1] <- paste(sim2_allele1[,i],sim2_allele2[,i],sep="/")
}

# Rename Columns to match Sample Names in RABBIT output
samples <-seq(1,sn,1)
samples <- paste("sim.sample",samples,sep="")
names(all_data)[2:ncol(all_data)] = samples

# Output File for AAA
write.csv(all_data,"known_AAA_simdata.csv",row.names = F)

