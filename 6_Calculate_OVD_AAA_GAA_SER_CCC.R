# Calculate AAA, AAG, PAA, and CCO for RABBIT OVD output
#Need to run lines 1:157 first, all calculations based on datalist3 (formatted output from RABBIT and known truth)
#Running entire script will output 3 tables to wd: 1) SPEARS metrics by sample (AAA, GAA, PAA, CCC),
#2) SPEARS metrics by marker (AAA, GAA) that includes missing and error distributions, 3) output from Pearson's cor.test on CO counts
#Tables on SPEARS metrics are used for genome plots

### Required Packages
library("dplyr")
library("tidyr")
library("data.table")
library("gdata")

#Designate working directory
wd <- "~/working/directory/"
setwd(wd)

#User inputs (read from file created in 2_SAEGUS_to_MACH_format.R
user <- read.table("user_input.txt", header = F, sep = "\t", stringsAsFactors = F)

# Name of parent key data
#fd <- user[user[[1]]=="fd",2]
# Number of chromosomes
chrom <- as.numeric(user[user[[1]]=="chrom",2])
# Name of Population (for MACH output)
popID <- user[user[[1]]=="popID",2]
#Sample Number
sn <- as.numeric(user[user[[1]]=="sn",2])
#Rsq Threshold
rsq <- as.numeric(user[user[[1]]=="rsq",2])
#Known parent-of-origin
ld <- user[user[[1]]=="ld",2]
#Known GT Data
kd <- user[user[[1]]=="kd",2]
#Number of parents
fn <- as.numeric(user[user[[1]]=="fn",2])

# parent key, known GT, and known parent-of-origin
#Known GT (subset parent key from this data)
GT_known <- read.csv(kd, head=T, stringsAsFactors = FALSE) #Known GT data created from script 1_SAEGUS_to_MACH_format.R
err_miss_dist <- GT_known[,c(1,5,6)] #distribution of genotyping error and missing data
key <- GT_known[,1:(8+fn)] #parent key
GT_known <- melt(as.data.table(GT_known), id.vars = c("snpID","CHROM"), measure.vars = patterns("^sim.sample")) #Reformat known data into long format
GT_known <- rename.vars(GT_known,c("variable","value"),c("sample","sim_GT")) #Assign variable names to melted columns

#Format parent key
key[,9:ncol(key)] <- lapply(key[,9:ncol(key)], gsub, pattern='/[0-9]', replacement='') #Remove second allele from GT, (0/0 to 0) to be used for inferring GT based on parent of origin

#Known parent-of-origin from SAEGUS output
origin_known <- read.csv(ld, head=T, stringsAsFactors = FALSE)
# Replace index with name of parent
# 1=ParentA, 2=ParentB, 3=ParentC, 4=ParentD, 5=ParentE, 6=ParentF, 7=ParentG, key for test data
parents <- names(key[,9:ncol(key)]) #parent names taken from parent key
for (i in 1:length(parents)){
  origin_known[which(origin_known$parent_of_origin==i),4] <- parents[i]
}

origin_known <- dcast(as.data.table(origin_known), sample + snpID ~ allele) #expand origin data so there is a column for each homolog
origin_known <- rename.vars(origin_known, c("snpID","a1","a2"), c("snpID","sim_parent1","sim_parent2"))
origin_known <- merge(origin_known,key[,1:2],all.x = T,by="snpID", sort=F) #Add chrom variable 
origin_known$sample <- paste0("sim.sample",origin_known$sample) #Rename samples to match RABBIT output
#####
# Format RABBIT OVD Output

# Diplotype key from RABBIT Output
#Based on output name from RABBIT output script
#Number of rows to skip when reading in data changes depending on number of parents and number of samples, 
#the formula 9 + sn +2*fn accounts for this. 9 is the number of headers (doesn't change in RABBIT output, 
#sn is sample number and accounts for rows that have likelihood output for each samples, and 2*fn accounts 
#for the rows that have the original parent GT data and haplotype columns
#
d.key <- read.table(paste("./RABBIT/reconstructed_chrom_1_jointModel_OVD_summary.csv", sep=""),head=T,
                             sep=",",skip=(9 + sn + 2*fn),nrows = (fn*fn), stringsAsFactors = F)
d.key$Diplotype <- as.numeric(sub("diplotype","", d.key$Diplotype))

## SNP index information for each chromosome
#RABBIT has the original snpIDs that were assigned in 4_MACH_to_RABBIT_format.R, which are read in here
#Viterbi output, however, is based on the index of the SNP in the file, 
#so if snpID skips a number it wouldn't match what the viterbi output is. 
#This index "SNPkey[[i]]$index" accounts for this difference.
SNPkey <- list()
for (c in 1:chrom){
SNPkey[[c]] <- read.table(paste("./RABBIT/reconstructed_chrom_",c,"_jointModel_OVD_summary.csv",sep=""),head=F,
                    sep=",", skip=1, nrows=2, stringsAsFactors = F)
SNPkey[[c]] <- as.data.frame(t(SNPkey[[c]][,2:ncol(SNPkey[[c]])]))
SNPkey[[c]] <- rename.vars(SNPkey[[c]],c("V1","V2"),c("SNP","CHROM"))
SNPkey[[c]]$index <- seq(1,nrow(SNPkey[[c]]),1)
}


### Format Viterbi Paths from RABBIT output for calculating AAA,AAG

datalist3 <- list() #create empty list to store formatted output from RABBIT and SAEGUS

#Iterate through each sample and read in the Viterbi Path and reformat the path into ranges and diplotypes
for (c in 1:chrom){
  
  datalist <- list() #create empty list to store each sample for each chromosome
  key_sub2 <- subset(key, key$CHROM == c) #subset parent data needed for chromosome c
  key_sub2 <- key_sub2[which(key_sub2$snpID %in% SNPkey[[c]]$SNP),] #subset parent data to only include SNPs in RABBIT output (post rsq filtering)

  
for(i in 1:sn){
  
  vit_path <- read.table(paste("./RABBIT/reconstructed_chrom_",c,"_jointModel_OVD_summary.csv",sep=""),head=T,
                              sep=",", skip=(11 + sn + 2*fn + (fn*fn)), nrows=sn, stringsAsFactors = F)

  #Subset sample i and Separate viterbi path
  
  dat <- as.data.frame(unlist(strsplit(vit_path[i,2],split = "-")))
  #Viterbi output is like this: (1-4-50-3-60-2…). Which translates to Diplotype 4 is present for markers with the index 1through 50, diplotype 3 for index 51-60, etc… When we run the line above, this unlists and separates data into index, diplotype, index, etc...along rows
  even <- seq(2,nrow(dat)-1,2) #This will be used to extract the diplotype from dat
  odd <- seq(1,nrow(dat)-1,2) #This will be used to extract start index of each diplotype range from dat
  #Extract diplotype and start index and format into dataframe
  dat <- data.frame(index_start=(dat$`unlist(strsplit(vit_path[i, 2], split = "-"))`[odd]), 
                    diplotype=(dat$`unlist(strsplit(vit_path[i, 2], split = "-"))`[even]))
  dat$index_start <- as.numeric(as.character(dat$index_start)) #make sure index is numeric
  dat$diplotype <- as.numeric(as.character(dat$diplotype)) #make sure diplotype is numeric
  dat <- as.data.table(dat) #needs to be data.table for next lines of code
  #Assign stop index based on start index of next haplotype (start index - 1)
  dat[,index_stop := shift(index_start, fill = first(index_start), type="lead")] 
  dat$index_stop <- dat$index_stop-1
  dat[dat==0] <- nrow(SNPkey[[c]]) #last stop position is assigned "0", so replace with end index (total number of SNPs for chromosome)
  dat <- merge(dat, d.key, by.x="diplotype", by.y="Diplotype") #merge so snpIDs can be referenced
  dat <- merge(dat, SNPkey[[c]], by.x="index_start",by.y="index") #replace start index with snpID
  dat <- rename.vars(dat,"SNP","SNP_start")
  dat <- merge(dat, SNPkey[[c]][,c(1,3)], by.x="index_stop",by.y="index", all.x=T) #replace stop index with snpID
  dat <- rename.vars(dat,"SNP","SNP_stop")
  dat$sample <- vit_path[i,1] #rename samples
  datalist[[i]] <- dat #assign sample to list
  
  #create new datatable that has a row for each snpID instead of start and stop on same row
  datalist[[i]] <- do.call(rbind, 
                       lapply(seq(nrow(datalist[[i]])), 
                              function(x){data.frame(
                                snpID=seq(as.numeric(datalist[[i]][x,"SNP_start"]), as.numeric(datalist[[i]][x,"SNP_stop"]) ,by=1),
                                sample = as.character(datalist[[i]][x,"sample"]),
                                CHROM = as.numeric(datalist[[i]][x,"CHROM"]),
                                OVD_parent = as.character(datalist[[i]][x,"founder"])
                              ) } ))
  
  #Separate diplotypes for easy comparison
  datalist[[i]] <- separate(datalist[[i]], OVD_parent, c("OVD_parent1", "OVD_parent2"), sep="\\|")
  datalist[[i]]$OVD_allele1 <- datalist[[i]]$OVD_parent1 #create new column to store parent assignment, will be reformatting OVD_parent1 to GT in next steps
  datalist[[i]]$OVD_allele2 <- datalist[[i]]$OVD_parent2 #create new column to store parent assignment
  datalist[[i]] <- datalist[[i]][which(datalist[[i]]$snpID %in% key_sub2$snpID),] #subset based on snps actually analyzed (datframe is created by filling in gaps between start and stop, some snps don't exist)
  
  for (j in 1:length(parents)){
  
    datalist[[i]][,6:7] <- apply(as.data.frame(datalist[[i]][,6:7]),2,function(x) 
    ifelse(x==names(key_sub2)[j+8],as.character(key_sub2[,(j+8)]),x)) #reformat parent assignment into genotypes
  }
  
}
  #now that samples are formatted in a way we can use, append all samples within a chromosome and reformat chromosome data into long format
  datalist3[[c]] <- as.data.frame(do.call("rbind",datalist))
  datalist3[[c]] <- datalist3[[c]][which(datalist3[[c]]$snpID %in% SNPkey[[c]]$SNP),] #subset only markers that met the original rsq threshold
  
  #Now merge known data with inferred data for calculating statistics, this will have known parent-of-origin (sim_parent1, sim_parent2), inferred parent-of-origin (OVD_parent1,OVD_parent2), known GT (sim_allele1,sim_allele2), and inferred GT (OVD_allele1,OVD_allele2)
  datalist3[[c]] <- merge(datalist3[[c]],subset(origin_known, origin_known$CHROM==c), sort = F, all.x=T)
  datalist3[[c]] <- merge(datalist3[[c]],subset(GT_known, GT_known$CHROM==c), sort = F, all.x=T)
  datalist3[[c]] <- separate(datalist3[[c]], sim_GT, c("sim_allele1", "sim_allele2"), sep="/")
  
  
}

########
## Calculate AAA
########

for_AAA <- datalist3
AAA_bymarker <- list()
AAA_bysample <- list()
for (c in 1:chrom){
  
  #calculate the proportion of matching sites, because this is phase independent it takes into account both possible orientations for het sites. Also, this includes sites that may only match 1 of 2 parents (assigned 0.5 as a half match)
  for_AAA[[c]]$parent_match <- ifelse(paste(for_AAA[[c]]$OVD_parent1,"|", for_AAA[[c]]$OVD_parent2,sep="") == paste(for_AAA[[c]]$sim_parent1,"|",for_AAA[[c]]$sim_parent2,sep=""), 1,
                                     ifelse(paste(for_AAA[[c]]$OVD_parent1,"|", for_AAA[[c]]$OVD_parent2,sep="") == paste(for_AAA[[c]]$sim_parent2,"|",for_AAA[[c]]$sim_parent1,sep=""),1,
                                            ifelse((for_AAA[[c]]$OVD_parent1 == for_AAA[[c]]$sim_parent1 | for_AAA[[c]]$OVD_parent1 == for_AAA[[c]]$sim_parent2),0.5,
                                                   ifelse((for_AAA[[c]]$OVD_parent2 == for_AAA[[c]]$sim_parent1 | for_AAA[[c]]$OVD_parent2 == for_AAA[[c]]$sim_parent2),0.5,0))))

  #Aggregate match data by sample and by marker and calculate mean AAA for each chromosome
  AAA_bymarker[[c]] <- aggregate(for_AAA[[c]][,12], list(for_AAA[[c]]$snpID), mean)
  AAA_bymarker[[c]] <- rename.vars(AAA_bymarker[[c]],"Group.1","snpID")
  AAA_bysample[[c]] <- aggregate(for_AAA[[c]][,12], list(for_AAA[[c]]$sample), mean)
  AAA_bysample[[c]] <- rename.vars(AAA_bysample[[c]],"Group.1","sample")


}

# AAA by marker
#Append all chromosomes
AAA_bymarker <- do.call("rbind",AAA_bymarker)
AAA_bymarker <- aggregate(AAA_bymarker[,2], list(AAA_bymarker$snpID), mean)
AAA_bymarker <- rename.vars(AAA_bymarker,c("Group.1","x"),c("snpID","AAA"))

# AAA by sample
#Append all chromosomes and calculate mean per sample
AAA_bysample <- do.call("rbind",AAA_bysample)
AAA_bysample <- aggregate(AAA_bysample[,2], list(AAA_bysample$sample), mean)
AAA_bysample <- rename.vars(AAA_bysample,c("Group.1","x"),c("sample","AAA"))

# Genomewide mean and sd of AAA for sample data
mean(AAA_bysample$AAA)
sd(AAA_bysample$AAA)

#####
## Calculate GAA
#####

#Create GAA from datalist for calculations and create empty lists to store GAA calculations
GAA <- datalist3
GAA_bymarker <- list()
GAA_bysample <- list()
for (c in 1:chrom){

    
    #Calculate matches between known and inferred genotype, same criteria as AAA
    GAA[[c]]$GT_match <- ifelse(GAA[[c]]$OVD_allele1 == GAA[[c]]$sim_allele1 & GAA[[c]]$OVD_allele2 == GAA[[c]]$sim_allele2,1,
                                 ifelse(GAA[[c]]$OVD_allele1 == GAA[[c]]$sim_allele2 & GAA[[c]]$OVD_allele2 == GAA[[c]]$sim_allele1,1,
                                        ifelse(GAA[[c]]$OVD_allele1 == GAA[[c]]$sim_allele1 | GAA[[c]]$OVD_allele2 == GAA[[c]]$sim_allele2 | GAA[[c]]$OVD_allele1 == GAA[[c]]$sim_allele2 | GAA[[c]]$OVD_allele2 == GAA[[c]]$sim_allele1,0.5,0)))
  
    #Aggregate match data by sample and by marker and calculate mean GAA for each chromosome
    GAA_bymarker[[c]] <- aggregate(GAA[[c]][,12], list(GAA[[c]]$snpID), mean)
    GAA_bymarker[[c]] <- rename.vars(GAA_bymarker[[c]],"Group.1","snpID")
    GAA_bysample[[c]] <- aggregate(GAA[[c]][,12], list(GAA[[c]]$sample), mean)
    GAA_bysample[[c]] <- rename.vars(GAA_bysample[[c]],"Group.1","sample")  

  }


#Append all chromosomes
GAA_bymarker <- do.call("rbind",GAA_bymarker)
GAA_bymarker <- aggregate(GAA_bymarker[,2], list(GAA_bymarker$snpID), mean)
GAA_bymarker <- rename.vars(GAA_bymarker,c("Group.1","x"),c("snpID","GAA"))

# GAA by sample
#Append all chromosomes and calculate mean per sample
GAA_bysample <- do.call("rbind",GAA_bysample)
GAA_bysample <- aggregate(GAA_bysample[,2], list(GAA_bysample$sample), mean)
GAA_bysample <- rename.vars(GAA_bysample,c("Group.1","x"),c("sample","GAA"))

# Genomewide mean and sd of GAA for sample data
mean(GAA_bysample$GAA)
sd(GAA_bysample$GAA)

#####
## Calculate PAA
#####

###Switch Accuracy (Lin et al. 2002): (n-1-sw)/(n-1); n=#het sites, sw=# switches to get correct phase
###Switch error rate (Stevens and Donnely, 2003): 1-switch accuracy
###SER is based on heterozygous sites, so does not perform well for inbred/homozygous samples
###Subset heterozygous sites that are accurate (AAG == 1)  and mark switch sites (those where phase switches
list_df <- list()
for (i in 1:chrom){
  list_df[[i]] <- subset(GAA[[i]],GAA[[i]]$GT_match==1)
  #list_df[[i]] <- AAG[[i]]
  list_df[[i]] <- subset(list_df[[i]],list_df[[i]]$sim_allele1!=list_df[[i]]$sim_allele2)
  list_df[[i]]$switch <- ifelse(list_df[[i]]$sim_parent1==list_df[[i]]$OVD_parent1 & list_df[[i]]$sim_parent2==list_df[[i]]$OVD_parent2,
                                0,1)
}

#This will calculate runs of markers that have the same switch set as the markers that precede it, so switches are only counted if it switches related to preceding het sites
#calculated for each chromosome (j) in each sample (i)
samples <- as.data.frame(as.character(unique(datalist3[[1]]$sample)),stringsAsFactors = FALSE)
allrunsum = list()
for (j in 1:chrom){
  runsum <- list()
  for (i in samples[,1]){
    counter <- vector(length=0)
    counter <- sequence(rle(as.character(subset(list_df[[j]],sample==i)$switch))$lengths)
    runsum[[i]] <- (nrow(subset(list_df[[j]],sample==i))-1-sum(counter==1))/(nrow(subset(list_df[[j]],sample==i))-1)
    runsum2 <- t(as.data.frame(runsum))
  }
  allrunsum[[j]] <- runsum2
}

#Combine all chromosomes
allrunsum2 <- do.call("cbind",allrunsum)

#Average for each sample
PAA <- as.data.frame(rowMeans(allrunsum2))
#Add sample names to dataframe
PAA$sample <- row.names(allrunsum2)
names(PAA) <- c("PAA","sample")
#mean and standard deviation
mean(PAA$PAA) #mean SER
sd(PAA$PAA) #sd of SER

#####
## Calculate CO counts and CCO for RABBIT output
#####

###Subset chromosomes, reformat for calculating CO
chrom_CO <- list()
for (i in 1:chrom) {
  CO <- datalist3[[i]]
  CO <- as.data.table(CO)
  CO <- melt(CO, id.vars = c("snpID","CHROM","sample"), measure.vars = c("OVD_parent1","OVD_parent2"))
  CO <- CO[order(CO$sample,as.numeric(CO$snpID)),]
  CO <- rename.vars(CO, c("variable","value"),c("homolog","parent"))
  CO$homolog <- ifelse(CO$homolog=="OVD_parent1",1,2)
  chrom_CO[[i]] <- CO
}

###For each homolog (b) in chromosome (j) for each sample (i), calculate the number of times the parent changes (a crossover) and sum total crossovers
samples <- as.data.frame(unique(chrom_CO[[1]]$sample),stringsAsFactors = FALSE)
look_run = NULL
for (j in 1:chrom){
  runsum = list()
  for (i in samples[,1]){
    runs <- vector(length=0)
    for (b in 1:2)
    {
      runs= c(runs, length(rle(subset(chrom_CO[[j]], sample==i & homolog==b)$parent)$lengths))
    }
    runsum[[i]] <- sum(runs)
  }
  look_run = rbind(look_run,data.frame(CHROM=j, runsum))
}

# Format output
OVD_CO_counts <- look_run
OVD_CO_sums <- colSums(OVD_CO_counts[2:(sn+1)]) #Calculate total for each chromosome
OVD_CO_sums <- as.data.frame(OVD_CO_sums)
OVD_CO_sums$sample <- row.names(OVD_CO_sums)
OVD_CO_sums$sample2 <- row.names(OVD_CO_sums)
OVD_CO_sums$sample2 <- sub("sim.sample","", OVD_CO_sums$sample2)
OVD_CO_sums <- OVD_CO_sums[order(as.numeric(OVD_CO_sums$sample2)),]
OVD_CO_sums <- OVD_CO_sums[,1:2]
names(OVD_CO_sums) <- c("CO_RABBIT","sample")
 
# Mean and sd
mean(OVD_CO_sums$CO_RABBIT) 
sd(OVD_CO_sums$CO_RABBIT)

#####
## Calculate CO counts and CCO for known data for comparison
#####

###Same process as simulated data
# Known CO counts for comparison
chrom_CO <- list()
for (i in 1:chrom) {
  CO <- datalist3[[i]]
  CO <- as.data.table(CO)
  CO <- melt(CO, id.vars = c("snpID","CHROM","sample"), measure.vars = c("sim_parent1","sim_parent2"))
  CO <- CO[order(CO$sample,as.numeric(CO$snpID)),]
  CO <- rename.vars(CO, c("variable","value"),c("homolog","parent"))
  CO$homolog <- ifelse(CO$homolog=="sim_parent1",1,2)
  chrom_CO[[i]] <- CO
}

samples <- as.data.frame(unique(chrom_CO[[1]]$sample),stringsAsFactors = FALSE)
look_run = NULL
for (j in 1:chrom){
  runsum = list()
  for (i in samples[,1]){
    runs <- vector(length=0)
    for (b in 1:2)
    {
      runs= c(runs, length(rle(subset(chrom_CO[[j]], sample==i & homolog==b)$parent)$lengths))
    }
    runsum[[i]] <- sum(runs)
  }
  look_run = rbind(look_run,data.frame(CHROM=j, runsum))
}

# Format output
known_CO_counts <- look_run
known_CO_sums <- colSums(known_CO_counts[2:(sn+1)])
known_CO_sums <- as.data.frame(known_CO_sums)
known_CO_sums$sample <- row.names(known_CO_sums)
known_CO_sums$sample2 <- row.names(known_CO_sums)
known_CO_sums$sample2 <- sub("sim.sample","", known_CO_sums$sample2)
known_CO_sums <- known_CO_sums[order(as.numeric(known_CO_sums$sample2)),]
known_CO_sums <- known_CO_sums[,1:2]
names(known_CO_sums) <- c("CO_known","sample")

# Mean and sd
mean(known_CO_sums$CO_known)
sd(known_CO_sums$CO_known)

# Merge known and RABBIT inferred
all_CO <- merge(known_CO_sums,OVD_CO_sums, sort=F)

# Pearson's correlation coefficient
cor(all_CO$CO_known,all_CO$CO_RABBIT)


##Output
#Table of per-sample metrics (AAA, GAA, PAA, CO_known, CO_RABBIT)
sample_table <- Reduce(merge, list(AAA_bysample,GAA_bysample, PAA, known_CO_sums, OVD_CO_sums))
write.table(sample_table,"SPEARS_by_sample_Metrics.csv",sep=",", row.names = F)
write.table(cbind("sample","SPEARS_by_sample_Metrics.csv"),"user_input.txt", col.names = F, sep="\t", append = TRUE, row.names = F)

#Table of per-marker metrics (AAA, GAA)
marker_table <- Reduce(merge, list(AAA_bymarker,GAA_bymarker,err_miss_dist))
write.table(marker_table,"SPEARS_by_marker_Metrics.csv",sep=",", row.names = F)
write.table(cbind("marker","SPEARS_by_marker_Metrics.csv"),"user_input.txt", col.names = F, sep="\t", append = TRUE, row.names = F)

#Pearson's Correlation Coefficient Output
corrFunc <- function(var1, var2, data) {
  result = cor.test(data[,var1], data[,var2])
  data.frame(var1, var2, result[c("estimate","p.value","statistic","method")], 
             stringsAsFactors=FALSE)
}
P_cor <- corrFunc("CO_known","CO_RABBIT",sample_table)
write.table(P_cor,"SPEARS_CO_Pearson_results.csv",sep=",", row.names = F)







