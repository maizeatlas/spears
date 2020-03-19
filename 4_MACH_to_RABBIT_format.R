# This takes mlinfo and mlgeno output from MACH and formats for use in RABBIT
# Filters markers that have an Rsq<0.80

# Required Packages
library("tidyr")
library("gdata")

#User inputs (wd,fd,chrom,popID)
# Working directory
wd <- "~/working/directory"
# Name of founder key data
fd <- "founder_key.txt"
# Number of chromosomes
chrom <- 10
# Name of Population (for MACH output)
popID <- "Pop"
#Sample Number
sn <- 1000
#Rsq Threshold
rsq <- 0.5


# MACH Output
# mlinfo files: Info on markers, SNP, each allele, Freq of allele1, MAF, Quality, Rsq
# If MACH is run in Linux, output will need to be unzipped prior to run
mlinfo <- list()
for (i in 1:chrom){
  
  setwd(paste0(wd,'/MACH/chrom',i))
  mlinfo[[i]]<- read.table(paste0("step2_chrom_",i,".mlinfo",sep=""), header=T, sep="\t")
  
}

mlinfo2 <- do.call(rbind,mlinfo)
hist(mlinfo2$Rsq) #Check R-square distribution from MACH output

mlinfo3 <- subset(mlinfo2, mlinfo2$Rsq>=rsq)
mlinfo3$SNP <- gsub("SNP",replacement = "",mlinfo3$SNP)

### Imputed genotypes
mlgeno <- list()
for (i in 1:chrom){
  
  setwd(paste0(wd,'/MACH/chrom',i))
  mlgeno[[i]]<- read.table(paste0("step2_chrom_",i,".mlgeno",sep=""), header=F, sep=" ", stringsAsFactors = FALSE)
  mlgeno[[i]]$V1 <- sub(paste(popID,"->",sep=""), "", mlgeno[[i]]$V1)
  mlgeno[[i]] <- mlgeno[[i]][,-2]
  
}

mlgeno_col <- mlgeno[[1]]
mlgeno_col <- cbind.data.frame(sample = mlgeno_col[,1])


for (i in 1:chrom){
  
  mlgeno_col <- cbind.data.frame(mlgeno_col,mlgeno[[i]][2:ncol(mlgeno[[i]])])  
  
}

mlgeno2 <- mlgeno_col
colnames(mlgeno2) <- c("sample", as.character(mlinfo2$SNP))
orig_mlgeno <- mlgeno2


#working directory: the one that contains founder data
setwd(wd)

#Founder Data
key <- read.table(fd, head=T, stringsAsFactors = FALSE, sep="\t")
key <- rename.vars(key,"chr","CHROM")
key <- key[which(key$snpID %in% mlinfo3$SNP),]
key2 <- separate(key, names(key[6]), c(paste(names(key[6]),".1",sep=""), paste(names(key[6]),".2",sep="")), sep="/")

##Formatting mlgeno data
names <- as.character(mlgeno2$sample)
look2 <- as.data.frame(t(mlgeno2[,2:ncol(mlgeno2)]))
names(look2) <- names
look2$snpID <- row.names(look2)
look2$snpID <- gsub("SNP",replacement = "",look2$snpID)
look2 <- look2[which(look2$snpID %in% mlinfo3$SNP),]
test <- merge(look2,key2[,c(1,2,4)],by="snpID")
test <- test[order(as.numeric(test$snpID)),]

#Run for sample columns
for (a in 2:(sn+1)) {
  test[,a] <- sub("/", "", test[,a])
}

#Create RABBIT folder
dir.create("RABBIT")

#Create header, founder, population data for RABBIT input, append, and write csv files for each chromosome
for (i in 1:chrom){
mlgeno <- subset(test,test$CHROM==i)
mlgeno <- mlgeno[order(as.numeric(mlgeno$snpID)),]
mlgeno <- as.data.frame(t(mlgeno[,2:(sn+1)]))
key3 <- subset(key2,key2$CHROM==i)

#Formats genotypes to 11,12,22 (1 is the allele present in first founder (ParentA in this case) and 2 is the alternate allele)
for (a in 1:ncol(mlgeno)) {
  mlgeno[,a] <- gsub(as.character(key3[,6][a]), "1", mlgeno[,a])
  mlgeno[,a] <- gsub("[A-Z]","2", mlgeno[,a])
}
look <- mlgeno
colnames(look) <- subset(test,test$CHROM==i)$snpID

#founders
mlgeno <- subset(key,key$CHROM==i)
mlgeno <- mlgeno[order(mlgeno$snpID),]
mlgeno2 <- mlgeno

#Subset founder columns
mlgeno <- as.data.frame(t(mlgeno[,6:12]))

#Formats genotypes to 11,12,22 (1 is the first allele in first founder and 2 is the alternate allele)
for (a in 1:ncol(mlgeno)) {
  mlgeno[,a] <- gsub(as.character(key3[,6][a]), "1", mlgeno[,a])
  mlgeno[,a] <- sub("/", "", mlgeno[,a])
  mlgeno[,a] <- gsub("[A-Z]","2", mlgeno[,a])
}

#Formats 11,12,22 to 1,2,N which is needed by RABBIT for founders
mlgeno <- apply(mlgeno, 2, function(x) sub("11", "1", x))
mlgeno <- apply(mlgeno, 2, function(x) sub("22", "2", x))
mlgeno <- apply(mlgeno, 2, function(x) sub("12", "N", x))
mlgeno <- apply(mlgeno, 2, function(x) sub("21", "N", x))
mlgeno <- as.data.frame(mlgeno)
colnames(mlgeno) <- colnames(look)

#Subset snpID, CHROM, and cM information
output.2 <- key[,c(1,2,4)]
output.2 <- subset(output.2,output.2$CHROM==i)
output.2 <- t(output.2)
rownames(output.2)[1:3] <- c("SNP", "CHROM", "cM")
colnames(output.2) <- colnames(look)

output.3 <- as.matrix(rbind(output.2, mlgeno, look))
write.table(output.3, paste0(wd,"/RABBIT/SimData_Rsq",rsq*100,"pct_chrom",i,"_RABBIT_input.csv",sep=""), sep=",",
            row.names=T,
            col.names=c("#founders", "7", rep("", (ncol(output.3)-2))), quote=F)
}

