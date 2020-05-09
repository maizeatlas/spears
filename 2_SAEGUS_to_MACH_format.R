# This takes output from SAEGUS and formats for MACH based on missingness data and founder genotypes
# Requires sim output from SAEGUS and founder key data (see file formats in README)

# Required Packages
library("tidyr")
library("gdata")
library("data.table")

# User inputs
# Working directory
wd <- "~/working/directory/"
# Name of simulated data
sd <- "simdata_n1000_test_set.csv"
# Name of founder key data
fd <- "test_founder_data_key.txt"
# Number of chromosomes
chrom <- 10
# Name of Population (for MACH output)
popID <- "Pop"
# Sample number
sn <- 1000

# Set working directory, which contains the simulated output from SAEGUS and the founder key data
setwd(wd)

# SimuPop Output
sim <- read.csv(sd,head=T, stringsAsFactors = FALSE) #Reads output from 1_SAEGUS.py
sim <- as.data.frame(t(sim[,3:ncol(sim)])) #Transposes data so markers are rows and samples are columns

# Founder Data
key <- read.table(fd, head=T, stringsAsFactors = FALSE,sep='\t') #Founder data key for GT info and marker info
key <- rename.vars(key,"chr","CHROM") #Rename because SAEGUS uses "chr" as the variable name

# Subset CHROM, POS, snpID, and GT for MACH founder.chrom.haplos files later on
founders <- key[,c(1:3,6:ncol(key))] #Subset "snpID", "CHROM", "POS", and parent columns
founders2 <- names(key[,6:ncol(key)]) #Extract founder ids

#Separate the genotypes into two columns for each founder
#These need to be completely homozygous
for (i in founders2) {
  key <- separate(key,i,c(paste(i,".1",sep=""),paste(i,".2",sep="")),sep="/")
}
key <- key[,c(1:5,seq(6,ncol(key),2))] #Only subset the first allele for each GT in each founder
key <- key[order(as.numeric(key$snpID)),] #Make sure markers are in numerical order by snpID


# Split simdata for formatting into separate alleles
sim2_allele1 <- sim[seq(1,nrow(sim),2),]
sim2_allele2 <- sim[seq(2,nrow(sim),2),]


# Format alleles to represent actual nucleotide and not generic number for founder
# 0=ParentA, 1=ParentB, 2=ParentC, 3=ParentD, 4=ParentE, 5=ParentF, 6=ParentG, key for test data

for (i in 1:length(founders2)){
  sim2_allele1 <- apply(sim2_allele1,2,function(x) ifelse(x==(i-1),key[,i+5],x))
}

for (i in 1:length(founders2)){
  sim2_allele2 <- apply(sim2_allele2,2,function(x) ifelse(x==(i-1),key[,i+5],x))
}

# Concatenate alleles to represent genotype in format G/G, A/G, etc...
all_data <- key[,1:5] #Subset snpID, CHROM, POS, cM, F_MISS to append to in loop
for (i in 1:ncol(sim2_allele1)){
  all_data[,i+5] <- paste(sim2_allele1[,i],sim2_allele2[,i],sep="/")
}

# Name samples sim.sample# and remove allele1 column
samples <-seq(1,sn,1) #Creates sequence of numbers for total number of samples
samples <- paste("sim.sample",samples,sep="") #Concatenate "sim.sample" to numbers in samples
names(all_data)[6:ncol(all_data)] = samples #Put sample names in GT matrix

# Format with missing info and final output for running MACH based on column in key data
# Removes rows with 100% missing data based on real data and then randomly samples samples to mask 
# for each marker
# This data is rounded based on 1000 samples, may need to change depending on your sample number

all_data$F_MISS_r <- (round(all_data$F_MISS,3)*sn) #Round the proportion missing to 3 decimal places
all_data <- subset(all_data,all_data$F_MISS_r<sn) #Remove markers with completely missing data, will be imputed bakc in later by MACH

#Randomly pick n number of samples based on the rounded missing proportion for each marker and set to missing
for (i in 1:nrow(all_data)){
  
  index_sub <- sample(1:sn,all_data[i,]$F_MISS_r)
  all_data[i,6:(ncol(all_data)-1)][index_sub] <- "./."
  
}

# Create MACH data, based on required format for MACH, see MACH help page
# Create MACH directory and folder for each chromosome within working directory
dir.create("MACH")
setwd("./MACH") #navigate to directory
# Create directory for each chromosome
for (i in seq(1,chrom)){
  dir.create(paste0('chrom',i))
}

# Create founder_map, which lists chrom and pos of founder data, population.dat, which lists the overlapping markers between founders and population
founders_map <- founders[,1:3]
founders_map <- founders_map[order(as.numeric(founders_map$snpID)),]

# Subset all_data to keep snpID, Chrom, Pos, and samples
all_data2 <- all_data[,-c(4:5,ncol(all_data))] #Removes cM, F_MISS, and F_MISS_r
founders_map[,1] <- sub("^","M SNP",founders_map[,1]) #Add "M_SNP" in front of marker ids in population (required format for MACH)
#This is the set of overlapping markers between the founder set and the population
population.dat <- founders_map[which(interaction(founders_map$CHROM, founders_map$POS) %in% 
                                       interaction(all_data2$CHROM,all_data2$POS)),1:2]

# Create .dat files for each chromosome from population.dat
for(i in seq(1,chrom)){
  setwd(paste0(wd,"/MACH/chrom",i))
  write.table(population.dat[population.dat$CHROM == i,][1], file = paste0("simdata_chrom_",i,".dat"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# Create .PED file (genotype data on specific chromosome for population)
all_data2 <- all_data2[with(all_data2,order(as.numeric(snpID))),] #Make sure order of markers is maintained
for (i in seq(1,chrom)){
  population.ped <- as.data.frame(t(all_data2[all_data2$CHROM == i,][,-c(1:3)]))
  population.ped <- data.frame(fID=popID, iID=rownames(population.ped), pID=0, mID=0, sex=0, population.ped) #Required IDs for MACH
  setwd(paste0(wd,"/MACH/chrom",i))
  write.table(population.ped, file = paste0("simdata_chrom_",i,".PED"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# Create founders_chrom_X.snps file (list of founder SNP#'s on specfic chromosome)
for(i in seq(1,chrom)){
  founders.snps <- data.frame(snpID=substring(founders_map[founders_map$CHROM == i,][,1],3))
  setwd(paste0(wd,"/MACH/chrom",i))
  write.table(founders.snps, file = paste0("founders_chrom_",i,".snps"), quote = FALSE,
              sep = "\t", row.names = FALSE,col.names = FALSE)
}

# Create founders_chrom_X.haplos file (founder haplotype data on specific chromosome)
for(i in seq(1,chrom)){
    founders.chrom.haps <- as.data.frame(t(founders[founders$CHROM == i,c(-1:-3)])) #transpose data and remove snpID, CHROM, and POS columns
  
  rownames(founders.chrom.haps) <- paste0(' ',rownames(founders.chrom.haps),' (',c(1:length(founders2)),')') #Set row names as founder IDs
  founders.chrom.haps <- as.data.frame(founders.chrom.haps[rep(1:nrow(founders.chrom.haps), each=2),]) #Creates two rows for each founder (each homolog)
  
  founders.chrom.haps <- data.frame("fID" = rownames(founders.chrom.haps),founders.chrom.haps, row.names = NULL)
  founders.chrom.haps <- cbind.data.frame('index' = rownames(founders.chrom.haps),
                                          "fID" =gsub('\\.1',"", founders.chrom.haps$fID),"HapID" = c("HAP1","HAP2"),
                                          founders.chrom.haps[2:ncol(founders.chrom.haps)]) #creates and formats dataframe
  
  # Removing the second allele from first haplotype and first allele from second haplotype
  founders.chrom.haps <- rbind.data.frame(lapply(founders.chrom.haps[c(TRUE,FALSE),], gsub, pattern='/[A-Z]', replacement=''),
                                          lapply(founders.chrom.haps[c(FALSE,TRUE),], gsub, pattern='[A-Z]/', replacement=''))
  # Reordering into correct order
  founders.chrom.haps <- founders.chrom.haps[order(as.numeric(levels(founders.chrom.haps$index)[founders.chrom.haps$index])),]
  setwd(paste0(wd,"/MACH/chrom",i))
  write.table(founders.chrom.haps[,-1],file = paste0('founders_chrom_',i,".haplos"),quote = FALSE, sep = '\t',
              row.names = FALSE,col.names = FALSE)
}




