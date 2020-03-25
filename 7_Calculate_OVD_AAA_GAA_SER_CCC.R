# Calculate AAA, AAG, SER, and CCO for RABBIT OVD output
# Need to create known truth data using script #6 first and have saved in working directory

### Required Packages
library("dplyr")
library("tidyr")
library("data.table")
library("gdata")

# User inputs
# Working directory
wd <- "~/working/directory"
# Name of founder key data
fd <- "founder_key.txt"
# Known Data (output from script 7)
kd <- "known_AAA_gbsonly.csv"
# Number of chromosomes
chrom <- 10
# Number of samples
sn <- 1000
# Rsq Threshold
rsq <- 0.5


# Set working directory
setwd(wd)

# Founder key and known truth from script 7
key <- read.table(fd, head=T, stringsAsFactors = FALSE, sep="\t")
key <- rename.vars(key,"chr","CHROM")
known <- read.csv(kd, head=T, stringsAsFactors = FALSE)
known <- merge(key[,1:2],known, sort = F)
known <- melt(as.data.table(known), id.vars = c("snpID","CHROM"), measure.vars = patterns("^sim.sample"))
known <- rename.vars(known,c("variable","value"),c("sample","sim_founder"))

#####
# Format RABBIT OVD Output
fn <- length(names(key))-5 # Number of founders
gtnum <- sum(seq(1,fn,1)) # Number of possible genotype combinations

# Diplotype key from RABBIT Output
d.key <- read.table(paste("./RABBIT/SimData_Rsq",rsq*100,"pct_chrom1_RABBIT_jointModel_OVD_output_magicReconstruct_Summary.csv", sep=""),head=T,
                             sep=",",skip=(9 + sn + 2*fn),nrows = (fn*fn), stringsAsFactors = F)
d.key$Diplotype <- as.numeric(sub("diplotype","", d.key$Diplotype))

## SNP index information for each chromosome
SNPkey <- list()
for (i in 1:chrom){
SNPkey[[i]] <- read.table(paste("./RABBIT/SimData_Rsq",rsq*100,"pct_chrom",i,"_RABBIT_jointModel_OVD_output_magicReconstruct_Summary.csv",sep=""),head=F,
                    sep=",", skip=1, nrows=2, stringsAsFactors = F)
SNPkey[[i]] <- as.data.frame(t(SNPkey[[i]][,2:ncol(SNPkey[[i]])]))
SNPkey[[i]] <- rename.vars(SNPkey[[i]],c("V1","V2"),c("SNP","CHROM"))
SNPkey[[i]]$index <- seq(1,nrow(SNPkey[[i]]),1)
}


### Format Viterbi Paths from RABBIT output for calculating AAA,AAG

datalist3 <- list()
for (c in 1:chrom){
  
  datalist <- list()

  
for(i in 1:sn){
  
  vit_path <- read.table(paste("./RABBIT/SimData_Rsq",rsq*100,"pct_chrom",c,"_RABBIT_jointModel_OVD_output_magicReconstruct_Summary.csv",sep=""),head=T,
                              sep=",", skip=(11 + sn + 2*fn + (fn*fn)), nrows=sn, stringsAsFactors = F)
  dat <- as.data.frame(unlist(strsplit(vit_path[i,2],split = "-")))
  even <- seq(2,nrow(dat)-1,2)
  odd <- seq(1,nrow(dat)-1,2)
  dat <- data.frame(index_start=(dat$`unlist(strsplit(vit_path[i, 2], split = "-"))`[odd]), 
                    diplotype=(dat$`unlist(strsplit(vit_path[i, 2], split = "-"))`[even]))
  dat$index_start <- as.numeric(as.character(dat$index_start))
  dat$diplotype <- as.numeric(as.character(dat$diplotype))
  dat <- as.data.table(dat)
  dat[,index_stop := shift(index_start, fill = first(index_start), type="lead")]
  dat$index_stop <- dat$index_stop-1
  dat[dat==0] <- nrow(SNPkey[[c]])
  dat <- merge(dat, d.key, by.x="diplotype", by.y="Diplotype")
  dat <- merge(dat, SNPkey[[c]], by.x="index_start",by.y="index")
  dat <- rename.vars(dat,"SNP","SNP_start")
  dat <- merge(dat, SNPkey[[c]][,c(1,3)], by.x="index_stop",by.y="index", all.x=T)
  dat <- rename.vars(dat,"SNP","SNP_stop")
  dat$sample <- paste("sim.sample",i,sep="")
  datalist[[i]] <- dat
  
}
  
  viterbi_all <- as.data.frame(do.call("rbind",datalist))
  datalist3[[c]] <- do.call(rbind, 
                            lapply(seq(nrow(viterbi_all)), 
                                   function(x){data.frame(
                                     snpID=seq(viterbi_all[x,"SNP_start"], viterbi_all[x,"SNP_stop"] ,by=1),
                                     sample = viterbi_all[x,"sample"],
                                     CHROM = viterbi_all[x,"CHROM"],
                                     OVD_founder = viterbi_all[x,"founder"]
                                   ) } ))
  
  datalist3[[c]] <- datalist3[[c]][which(datalist3[[c]]$snpID %in% SNPkey[[c]]$SNP),]
  
  known_sub <- subset(known, known$CHROM==c)
  
  datalist3[[c]] <- merge(datalist3[[c]],known_sub[,c(1,3,4)],sort = F, all.x=T)
  datalist3[[c]] <- separate(datalist3[[c]], sim_founder, c("sim_founder1", "sim_founder2"), sep="/")
  datalist3[[c]] <- separate(datalist3[[c]], OVD_founder, c("OVD_founder1", "OVD_founder2"), sep="\\|")
  
  
}

########
## Calculate AAA
########

for_AAA <- datalist3
AAA_bymarker <- list()
AAA_bysample <- list()
for (c in 1:chrom){
  
  for_AAA[[c]]$OVD_match <- ifelse(paste(for_AAA[[c]]$OVD_founder1,"|", for_AAA[[c]]$OVD_founder2,sep="") == paste(for_AAA[[c]]$sim_founder1,"|",for_AAA[[c]]$sim_founder2,sep=""), 1,
                                     ifelse(paste(for_AAA[[c]]$OVD_founder1,"|", for_AAA[[c]]$OVD_founder2,sep="") == paste(for_AAA[[c]]$sim_founder2,"|",for_AAA[[c]]$sim_founder1,sep=""),1,
                                            ifelse((for_AAA[[c]]$OVD_founder1 == for_AAA[[c]]$sim_founder1 | for_AAA[[c]]$OVD_founder1 == for_AAA[[c]]$sim_founder2),0.5,
                                                   ifelse((for_AAA[[c]]$OVD_founder2 == for_AAA[[c]]$sim_founder1 | for_AAA[[c]]$OVD_founder2 == for_AAA[[c]]$sim_founder2),0.5,0))))

  AAA_bymarker[[c]] <- aggregate(for_AAA[[c]][,8], list(for_AAA[[c]]$snpID), mean)
  AAA_bymarker[[c]] <- rename.vars(AAA_bymarker[[c]],"Group.1","snpID")
  AAA_bysample[[c]] <- aggregate(for_AAA[[c]][,8], list(for_AAA[[c]]$sample), mean)
  AAA_bysample[[c]] <- rename.vars(AAA_bysample[[c]],"Group.1","sample")


}

# AAA by marker
AAA_bymarker <- do.call("rbind",AAA_bymarker)
AAA_bymarker <- aggregate(AAA_bymarker[,2], list(AAA_bymarker$snpID), mean)
AAA_bymarker <- rename.vars(AAA_bymarker,c("Group.1","x"),c("snpID","AAA"))

# AAA by sample
AAA_bysample <- do.call("rbind",AAA_bysample)
AAA_bysample <- aggregate(AAA_bysample[,2], list(AAA_bysample$sample), mean)
AAA_bysample <- rename.vars(AAA_bysample,c("Group.1","x"),c("sample","AAA"))

# Overall mean and sd of AAA
mean(AAA_bysample$AAA)
sd(AAA_bysample$AAA)

write.csv(AAA_bysample,"AAA_by_sample_OVD.csv",row.names = F)
write.csv(AAA_bymarker,"AAA_by_marker_OVD.csv",row.names = F)

#####
## Calculate AAG
#####

AAG <- list()
key[,6:ncol(key)] <- lapply(key[,6:ncol(key)], gsub, pattern='/[A-Z]', replacement='')

founders <- names(key[,6:ncol(key)])
for (c in 1:chrom){
  key_sub2 <- key[which(key$snpID %in% SNPkey[[c]]$SNP),]
  temp <- list()
  for (i in 1:sn){
    
    temp[[i]] <- subset(datalist3[[c]],datalist3[[c]]$sample==paste("sim.sample",i,sep=""))
    temp[[i]] <- merge(temp[[i]],key_sub2[,c(1,6:ncol(key_sub2))], by="snpID", sort=F)
    temp[[i]] <- temp[[i]][order(as.numeric(temp[[i]]$snpID)),]
    
    for (j in 1:length(founders)){
      
      temp[[i]][,4:7] <- apply(as.data.frame(temp[[i]][,4:7]),2,function(x) ifelse(x==names(temp[[i]])[j+7],as.character(temp[[i]][,(j+7)]),x))
      
    }
    
    temp[[i]]$GT_match <- ifelse(temp[[i]]$OVD_founder1 == temp[[i]]$sim_founder1 & temp[[i]]$OVD_founder2 == temp[[i]]$sim_founder2,1,
                                 ifelse(temp[[i]]$OVD_founder1 == temp[[i]]$sim_founder2 & temp[[i]]$OVD_founder2 == temp[[i]]$sim_founder1,1,
                                        ifelse(temp[[i]]$OVD_founder1 == temp[[i]]$sim_founder1 | temp[[i]]$OVD_founder2 == temp[[i]]$sim_founder2 | temp[[i]]$OVD_founder1 == temp[[i]]$sim_founder2 | temp[[i]]$OVD_founder2 == temp[[i]]$sim_founder1,0.5,0)))
    

  }
  
  AAG[[c]] <- do.call("rbind",temp)
  
}

AAG2 <- do.call("rbind",AAG)

# AAG by marker
AAG_bymarker <- aggregate(AAG2$GT_match, list(AAG2$snpID), mean)
AAG_bymarker <- rename.vars(AAG_bymarker,c("Group.1","x"),c("snpID","AAG"))

# AAG by sample
AAG_bysample <- aggregate(AAG2$GT_match, list(AAG2$sample), mean)
AAG_bysample <- rename.vars(AAG_bysample,c("Group.1","x"),c("sample","AAG"))

# Overall mean and sd of AAG
mean(AAG_bysample$AAG) # mean AAG
sd(AAG_bysample$AAG) # sd of AAG by sample

write.csv(AAG_bysample,"AAG_by_sample_OVD.csv",row.names = F)
write.csv(AAG_bymarker,"AAG_by_marker_OVD.csv",row.names = F)


#####
## Calculate SER
#####

list_df <- list()
for (i in 1:chrom){
  list_df[[i]] <- subset(AAG[[i]],AAG[[i]]$GT_match==1)
  list_df[[i]] <- AAG[[i]]
  list_df[[i]] <- subset(list_df[[i]],list_df[[i]]$sim_founder1!=list_df[[i]]$sim_founder2)
  list_df[[i]]$switch <- ifelse(list_df[[i]]$sim_founder1==list_df[[i]]$OVD_founder1 & list_df[[i]]$sim_founder2==list_df[[i]]$OVD_founder2,
                                0,1)
}

samples <- as.data.frame(as.character(unique(datalist3[[1]]$sample)),stringsAsFactors = FALSE)
allrunsum = list()
for (j in 1:10){
  runsum <- list()
  for (i in samples[,1]){
    counter <- vector(length=0)
    counter <- sequence(rle(as.character(subset(list_df[[j]],sample==i)$switch))$lengths)
    runsum[[i]] <- (nrow(subset(list_df[[j]],sample==i))-1-sum(counter==1))/(nrow(subset(list_df[[j]],sample==i))-1)
    runsum2 <- t(as.data.frame(runsum))
  }
  allrunsum[[j]] <- runsum2
}

allrunsum2 <- do.call("cbind",allrunsum)
SER <- as.data.frame(rowMeans(allrunsum2))
SER$sample <- row.names(allrunsum2)
names(SER) <- c("SER","sample")
SER$SER <- 1-SER$SER
mean(SER$SER) #mean SER
sd(SER$SER) #sd of SER
write.csv(SER,"SER_OVD_bysample.csv",row.names = F)


#####
## Calculate CO counts and CCO for RABBIT output
#####

chrom_CO <- list()
for (i in 1:chrom) {
  CO <- datalist3[[i]]
  CO <- melt(CO, id.vars = c("snpID","CHROM","sample"), measure.vars = c("OVD_founder1","OVD_founder2"))
  CO <- CO[order(CO$sample,as.numeric(CO$snpID)),]
  CO <- rename.vars(CO, c("variable","value"),c("homolog","founder"))
  CO$homolog <- ifelse(CO$homolog=="OVD_founder1",1,2)
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
      runs= c(runs, length(rle(subset(chrom_CO[[j]], sample==i & homolog==b)$founder)$lengths))
    }
    runsum[[i]] <- sum(runs)
  }
  look_run = rbind(look_run,data.frame(CHROM=j, runsum))
}

# Format output
OVD_CO_counts <- look_run
OVD_CO_sums <- colSums(OVD_CO_counts[2:(sn+1)])
OVD_CO_sums <- as.data.frame(OVD_CO_sums)
OVD_CO_sums$sample <- row.names(OVD_CO_sums)
OVD_CO_sums$sample2 <- row.names(OVD_CO_sums)
OVD_CO_sums$sample2 <- sub("sim.sample","", OVD_CO_sums$sample2)
OVD_CO_sums <- OVD_CO_sums[order(as.numeric(OVD_CO_sums$sample2)),]
OVD_CO_sums <- OVD_CO_sums[,1:2]
names(OVD_CO_sums) <- c("OVD_CO","sample")
 
# Mean and sd
mean(OVD_CO_sums$OVD_CO) 
sd(OVD_CO_sums$OVD_CO)

#####
## Calculate CO counts and CCO for known data for comparison
#####

# Known CO counts for comparison
chrom_CO <- list()
for (i in 1:chrom) {
  CO <- datalist3[[i]]
  CO <- melt(CO, id.vars = c("snpID","CHROM","sample"), measure.vars = c("sim_founder1","sim_founder2"))
  CO <- CO[order(CO$sample,as.numeric(CO$snpID)),]
  CO <- rename.vars(CO, c("variable","value"),c("homolog","founder"))
  CO$homolog <- ifelse(CO$homolog=="sim_founder1",1,2)
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
      runs= c(runs, length(rle(subset(chrom_CO[[j]], sample==i & homolog==b)$founder)$lengths))
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
names(known_CO_sums) <- c("known_CO","sample")

# Mean and sd
mean(known_CO_sums$known_CO)
sd(known_CO_sums$known_CO)

# Merge known and RABBIT inferred
all_CO <- merge(known_CO_sums,OVD_CO_sums, sort=F)

# Pearson's correlation coefficient
cor(all_CO$known_CO,all_CO$OVD_CO)
write.csv(all_CO, "CO_counts_known_OVD.csv",row.names = F)




