### Format RABBIT output from the origPosteriorDecoding Algorithm for plotting
### Outputs the difference in maximum likelihood of the top two founders at each marker for each sample

### Required Packages
library("dplyr")
library("tidyr")
library("data.table")
library("gdata")
library("ggplot2")

# User inputs
# Working directory (contains RABBIT folder)
wd <- "/working/directory"
# Name of founder key data
fd <- "founder_key.txt"
# Number of chromosomes
chrom <- 10
# Number of samples
sn <- 1000
#Rsq threshold
rsq <- 0.5

# Set working directory
setwd(wd)

#Founder Data
key <- read.table(fd, head=T, stringsAsFactors = FALSE, sep="\t")
key <- rename.vars(key,"chr","CHROM")

#Pull mlinfo data from MACH that contains Rsq data
mlinfo <- list()
for (i in 1:chrom){
  
  setwd(paste0(wd,'/MACH/chrom',i))
  mlinfo[[i]]<- read.table(paste0("step2_chrom_",i,".mlinfo",sep=""), header=T, sep="\t")
  
}

mlinfo2 <- do.call(rbind,mlinfo)
hist(mlinfo2$Rsq) #Check R-square distribution from MACH output

mlinfo3 <- subset(mlinfo2, mlinfo2$Rsq>=rsq) #Subset based on R-square threshold
mlinfo3$SNP <- gsub("SNP",replacement = "",mlinfo3$SNP) #Format SNP column

key <- key[which(key$snpID %in% mlinfo3$SNP),] #subset key based on R-square threshold
key2 <- separate(key, names(key[6]), c(paste(names(key[6]),".1",sep=""), paste(names(key[6]),".2",sep="")), sep="/") #Separate alleles in GT for parent1, used for fomatting later

setwd(wd)

# Format RABBIT OPD Output
fn <- length(names(key))-5 # Number of founders
gtnum <- sum(seq(1,fn,1)) # Number of possible genotype combinations
test_acc <- list()
for (c in 1:chrom){
  ## Probabilities from RABBIT for each chromosome
  haplo_prob <- read.table(paste("./RABBIT/SimData_Rsq",rsq*100,"pct_chrom1_RABBIT_jointModel_OPD_output_magicReconstruct_Summary.csv", sep=""),head=T,
                           sep=",",skip=(9 + gtnum + sn + fn),nrows = ((gtnum*sn)+2), stringsAsFactors = F)
  haplo_prob <- haplo_prob[-c(1:2),]
  haplo_prob <- separate(haplo_prob, SNP, c("Sample", "GT"), sep="_")
  
  ## Founder key from RABBIT output
  haplo_prob_key <- read.table(paste("./RABBIT/SimData_Rsq",rsq*100,"pct_chrom1_RABBIT_jointModel_OPD_output_magicReconstruct_Summary.csv", sep=""),head=T,
                               sep=",",skip=(7 + sn + fn),nrows = gtnum, stringsAsFactors = F)

  ## Assign ancester for each sample at each marker based on RABBIT output
  prob_t2 <- list()
  for (i in 1:sn) {
      sample1 <- subset(haplo_prob,haplo_prob$Sample==paste("sim.sample",i,sep="")) #subset sample
      sample1 <- melt(sample1, id.vars = c("Sample","GT")) #format into long format
      sample1 <- rename.vars(sample1, c("variable","value", "Sample"),c("snpID","ML","sample")) #change variable names
      sample1$ML <- as.numeric(sample1$ML) #make sure probabilities are numeric
    
    sample1 <- sample1 %>% 
      arrange(desc(ML)) %>% 
      group_by(snpID) %>% slice(1:2) #slices data (which is the probability for all founders, and outputs the top two founders
    
    sample1_m1 <- sample1[seq(1,nrow(sample1),2),] #top founder
    sample1_m2 <- sample1[seq(2,nrow(sample1),2),] #second founder
    
    sample1_m1 <- rename.vars(sample1_m1,c("ML","GT"),c("ml1","GT_ml1")) #rename
    sample1_m2 <- rename.vars(sample1_m2,c("ML","GT"),c("ml2","GT_ml2")) #rename
    
    sample1 <- cbind(as.data.frame(sample1_m1),as.data.frame(sample1_m2)) #combines data so there is a separate column for top founder and second founder
    
    sample1$snpID <- as.numeric(sub("X","", sample1$snpID)) #reformat snpID
    sample1$diff <- as.numeric(sample1$ml1)-as.numeric(sample1$ml2) #difference between top two founders
    prob_t2[[i]] <- sample1 #output sample data
    
  }
  
  test_acc[[c]] <- as.data.frame(do.call("rbind",prob_t2)) #combine all sample data
  
}

bymarker <- list()
for (i in 1:chrom){
  
  bymarker[[i]] <- aggregate(test_acc[[i]]$diff,list(test_acc[[i]]$snpID),mean) #take average probabilities across all samples for each marker in each chromosome
  
}

ml_bymarker <- do.call("rbind",bymarker)
ml_bymarker <- rename.vars(ml_bymarker,c("Group.1","x"),c("snpID","ml_diff"))
setwd(wd)
write.csv(ml_bymarker,"ml_diff_bymarker_18MAR20.csv",row.names = F) #output results

