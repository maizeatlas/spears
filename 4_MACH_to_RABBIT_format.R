# This takes mlinfo and mlgeno output from MACH and formats for use in RABBIT
# Filters markers that have an Rsq<0.80

#Designate working directory
wd <- "~/Dropbox/Maize_ATLAS_share/ParallelSelection/GBS/Manuscripts/RABBIT_Bio_App/scripts/test_run/"

# Required Packages
library("tidyr")
library("gdata")

setwd(wd)

#User inputs (read from file created in 2_SAEGUS_to_MACH_format.R
user <- read.table("user_input.txt", header = F, sep = "\t")

# Name of founder key data
fd <- user[user[[1]]=="fd",2]
# Number of chromosomes
chrom <- as.numeric(user[user[[1]]=="chrom",2])
# Name of Population (for MACH output)
popID <- user[user[[1]]=="popID",2]
#Sample Number
sn <- as.numeric(user[user[[1]]=="sn",2])
#Rsq Threshold
rsq <- as.numeric(user[user[[1]]=="rsq",2])


# MACH Output
# mlinfo files: Info on markers, SNP, each allele, Freq of allele1, MAF, Quality, Rsq
# If MACH is run in Linux, output will need to be unzipped prior to run
mlinfo <- list()
#Read in mlinfo files for each chromosome
for (i in 1:chrom){
  
  setwd(paste0(wd,'/MACH/chrom',i))
  mlinfo[[i]]<- read.table(paste0("step2_chrom_",i,".mlinfo",sep=""), header=T, sep="\t")
  
}

mlinfo2 <- do.call(rbind,mlinfo) #Append all chromosomes together
hist(mlinfo2$Rsq) #Check R-square distribution from MACH output

mlinfo3 <- subset(mlinfo2, mlinfo2$Rsq>=rsq) #Subset mlinfo based on R-square threshold
mlinfo3$SNP <- gsub("SNP",replacement = "",mlinfo3$SNP) #Format SNP IDs for merging later

### Imputed genotypes from MACH
#Read in each mlgeno file for each chromosome and reformat
mlgeno <- list()
for (i in 1:chrom){
  
  setwd(paste0(wd,'/MACH/chrom',i))
  mlgeno[[i]]<- read.table(paste0("step2_chrom_",i,".mlgeno",sep=""), header=F, sep=" ", stringsAsFactors = FALSE)
  mlgeno[[i]]$V1 <- sub(paste(popID,"->",sep=""), "", mlgeno[[i]]$V1)
  mlgeno[[i]] <- mlgeno[[i]][,-2]
  
}

#Get sample ids from output
mlgeno_col <- mlgeno[[1]]
mlgeno_col <- cbind.data.frame(sample = mlgeno_col[,1])

#Combine sample IDs with corresponding genotype
for (i in 1:chrom){
  
  mlgeno_col <- cbind.data.frame(mlgeno_col,mlgeno[[i]][2:ncol(mlgeno[[i]])])  
  
}

#Assign sample IDs to corresponding columns
mlgeno2 <- mlgeno_col
colnames(mlgeno2) <- c("sample", as.character(mlinfo2$SNP))
orig_mlgeno <- mlgeno2


#working directory: the one that contains founder data
setwd(wd)

#Founder Data
key <- read.table(fd, head=T, stringsAsFactors = FALSE, sep="\t") #Upload founder data
key <- rename.vars(key,"chr","CHROM") #Rename chromosome column
key <- key[which(key$snpID %in% mlinfo3$SNP),] #Subset based on markers kept after R-square filter
#key2 <- separate(key, names(key[6]), c(paste(names(key[6]),".1",sep=""), paste(names(key[6]),".2",sep="")), sep="/") #Separate the GT data for the first founder into separate alleles. The first allele from the first founder will be used as allele "1" for RABBIT formatting

##Formatting mlgeno data
names <- as.character(mlgeno2$sample) #make sure sample ids are in character format
look2 <- as.data.frame(t(mlgeno2[,2:ncol(mlgeno2)])) #transpose data into markers as rows and samples as columns
names(look2) <- names #Set sample names
look2$snpID <- row.names(look2) #set marker names
look2$snpID <- gsub("SNP",replacement = "",look2$snpID) #Reformat snpID for RABBIT
look2 <- look2[which(look2$snpID %in% mlinfo3$SNP),] #Make sure same markers in each set are being used
test <- merge(key[,c(1,2,4)],look2,by="snpID") #Merge snpID, CHROM, and cM to genotype data
test <- test[order(as.numeric(test$snpID)),] #Check marker order

#Run for sample columns, reformat GT data to exclude "/"
for (a in 4:(sn+3)) {
  test[,a] <- sub("/", "", test[,a])
}

#Create RABBIT folder
dir.create("RABBIT")

#Create header, founder, population data for RABBIT input, append, and write csv files for each chromosome
for (i in 1:chrom){
    mlgeno <- subset(test,test$CHROM==i) #subset chromosome
    mlgeno <- mlgeno[order(as.numeric(mlgeno$snpID)),] #check marker order
    mlgeno <- as.data.frame(t(mlgeno[,4:(sn+3)])) #transpose data so samples are rows and markers are columns
    key3 <- subset(key,key$CHROM==i) #subset founder key for each chromosome

#Formats genotypes to 11,12,22 (1 is the allele present in first founder (ParentA in this case) and 2 is the alternate allele)
for (a in 1:ncol(mlgeno)) {
  mlgeno[,a] <- gsub(as.character(key3[,6][a]), "1", mlgeno[,a])
  mlgeno[,a] <- gsub("[A-Z]","2", mlgeno[,a])
}
look <- mlgeno
colnames(look) <- subset(test,test$CHROM==i)$snpID #Assign snpIDs

#founders
mlgeno <- subset(key,key$CHROM==i) #subset original founder key for each chromosome
mlgeno <- mlgeno[order(mlgeno$snpID),] #check marker order
mlgeno2 <- mlgeno

#Subset founder columns
mlgeno <- as.data.frame(t(mlgeno[,8:ncol(mlgeno)])) #transpose founder data to match pop data

#Formats genotypes to 11,12,22 (1 is the first allele in first founder and 2 is the alternate allele)
#for (a in 1:ncol(mlgeno)) {
# mlgeno[,a] <- gsub(as.character(key3[,7][a]), "1", mlgeno[,a])
#  mlgeno[,a] <- sub("/", "", mlgeno[,a])
#  mlgeno[,a] <- gsub("[A-Z]","2", mlgeno[,a])
#}

#Formats 11,12,22 to 1,2,N which is needed by RABBIT for founders
mlgeno <- apply(mlgeno, 2, function(x) sub("0/0", "1", x))
mlgeno <- apply(mlgeno, 2, function(x) sub("1/1", "2", x))
mlgeno <- apply(mlgeno, 2, function(x) sub("0/1", "N", x))
mlgeno <- apply(mlgeno, 2, function(x) sub("1/0", "N", x))
mlgeno <- as.data.frame(mlgeno)
colnames(mlgeno) <- colnames(look)

#Subset snpID, CHROM, and cM information
output.2 <- key[,c(1,2,4)]
output.2 <- subset(output.2,output.2$CHROM==i)
output.2 <- t(output.2)
rownames(output.2)[1:3] <- c("SNP", "CHROM", "cM")
colnames(output.2) <- colnames(look)

#Append data and write output 
output.3 <- as.matrix(rbind(output.2, mlgeno, look))
write.table(output.3, paste0(wd,"/RABBIT/SimData_Rsq",rsq*100,"pct_chrom",i,"_RABBIT_input.csv",sep=""), sep=",",
            row.names=T,
            col.names=c("#founders", "7", rep("", (ncol(output.3)-2))), quote=F)
}


