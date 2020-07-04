# This takes output from SAEGUS and formats for MACH based on missingness data and parent genotypes
# Requires sim output from SAEGUS and parent key data (see file formats in README)

# Required Packages
library("tidyr")
library("gdata")
library("data.table")

# User inputs  
# Working directory
wd <- "~/working/directory/"
# Simulated GT matrix
sd <- "known_GT.csv"
#Simulated Parent_of_origin
ld <- "known_parent_of_origin.csv"
# Name of parent data key
fd <- "parent_key_vcf.txt"
# Number of chromosomes
chrom <- 10
# Name of Population (for MACH output)
popID <- "Pop"
# Sample number
sn <- 100
#Rsq Threshold (for imputation output)
rsq <- 0.8
#Genotyping Error
g_error <- 0.006

# Set working directory, which contains the simulated output from SAEGUS and the parent key data
setwd(wd)

#write file that contains all user inputs for downstream scripts to use
write.table(rbind(wd,sd,ld,fd,chrom,popID,sn,rsq),"user_input.txt", col.names = F, sep="\t")


# parent Data
parents <- read.table(fd, head=T, stringsAsFactors = FALSE,sep='\t') #parent data key for GT info and marker info
parents <- rename.vars(parents,"chr","CHROM") #Rename because SAEGUS uses "chr" as the variable name

# SimuPop Output
sim <- read.csv(sd,head=T, stringsAsFactors = FALSE) #Reads in GT output from 1_SAEGUS.py

#Collapse alleles into Genotypes and transpose for MACH
out <- sim[1]
out[2:(((ncol(sim)-1)/2)+1)] <- lapply(seq(2, ncol(sim), 2), 
                                  function(i) do.call(paste0, c(sim[i],"/",sim[(i+1)])))
sim <- as.data.frame(t(out[2:ncol(out)]))
#Rename samples to sim.sample#
names(sim) <- paste("sim.sample",out$ind_id,sep="")
#Concatenate snpID, CHROM, POS, cM, and F_MISS data from parent key
sim <- cbind(parents[,1:7],sim)
known <- sim

#Genotyping error
#1 Create a dataframe (num_loci x num_indv) filled with random sampling from uniform distribution
err_df <- as.data.frame(matrix(runif(nrow(parents)*sn),nrow(parents),sn))
#2 sample GTs and change
#index positions in err_df < user defined ERR at that locus
gt_change <- which(err_df<g_error, arr.ind = T)
#locate indexed positions in known data and randomly change to one of two remaining GTs (equal prob)
sim[,8:ncol(sim)][gt_change] <- ifelse(sim[,8:ncol(sim)][gt_change]=="0/0", sample(c("0/1","1/1"),1), ifelse(sim[,8:ncol(sim)][gt_change]=="1/1", sample(c("0/0","0/1"),1), sample(c("0/0","1/1"))))


#Calculate realized error at each locus
err_df$R_ERR <- apply(err_df,1,function(x) length(which(x<g_error))/sn)
#Histogram of realized error
hist(err_df$R_ERR, main = "Realized Genotyping Error Across all Loci", xlab = "Realized GT Error")
#Append realized error distribution to known GT data set
known <- cbind.data.frame(known[,1:5],R_ERR = err_df[,(sn+1)], parents[,6:ncol(parents)], known[,8:ncol(known)])

#Output known GT data for later
write.csv(known,"known_GT2.csv",row.names = F)
nparents <- length(parents)-7 #Number of parents
write.table(cbind("kd","known_GT2.csv"),"user_input.txt", col.names = F, sep="\t", append = TRUE, row.names = F)
write.table(cbind("fn",nparents),"user_input.txt", col.names = F, sep="\t", append = TRUE, row.names = F)


#Induce missing data 
sim$F_MISS_r <- (round(sim$F_MISS,(nchar(sn)-1))*sn) #Create an integer for number of genotypes to set to missing (rounded to whole number based on original sample number)
sim <- subset(sim,sim$F_MISS_r<sn) #Remove markers with completely missing data, will be imputed back in later by MACH (also reduces time to run for next step)

#Randomly pick n number of samples based on the rounded missing proportion for each marker and set to missing
look <- split(as.matrix(sim[,8:(ncol(sim))]), row(as.matrix(sim[,8:(ncol(sim))])))
for (i in 1:nrow(sim)){
  
  look[[i]][sample(sn,look[[i]][sn+1])] <- "./."
  
}
look <- as.data.frame(do.call("rbind",look))
names <- (colnames(sim[1:ncol(sim)-1]))
sim <- cbind(sim[,1:7],look[,1:sn])
colnames(sim) <- names

#Format alleles in required format for MACH (A,C,T,G) for parents and population
#parents
look <- split(as.matrix(parents[,6:(ncol(parents))]), row(as.matrix(parents[,6:(ncol(parents))])))
for (i in 1:nrow(parents)){
  look[[i]] <- gsub("0/0", (paste0(as.character(look[[i]][1]), "/", as.character(look[[i]][1]))), look[[i]])
  look[[i]] <- gsub("1/1", (paste0(as.character(look[[i]][2]), "/", as.character(look[[i]][2]))), look[[i]])
}
look <- as.data.frame(do.call("rbind",look))
names <- colnames(parents)
parents <- cbind(parents[,1:5],look)
colnames(parents) <- names

#Population

look <- split(as.matrix(sim[,6:(ncol(sim))]), row(as.matrix(sim[,6:(ncol(sim))])))
for (i in 1:nrow(sim)){
  
  look[[i]] <- gsub("0/0", (paste0(as.character(look[[i]][1]), "/", as.character(look[[i]][1]))), look[[i]])
  look[[i]] <- gsub("1/1", (paste0(as.character(look[[i]][2]), "/", as.character(look[[i]][2]))), look[[i]])
  look[[i]] <- gsub("0/1", (paste0(as.character(look[[i]][1]), "/", as.character(look[[i]][2]))), look[[i]])
  look[[i]] <- gsub("1/0", (paste0(as.character(look[[i]][1]), "/", as.character(look[[i]][2]))), look[[i]])
  
}
look <- as.data.frame(do.call("rbind",look))
names <- (colnames(sim))
sim <- cbind(sim[,1:5],look)
colnames(sim) <- names


### Create MACH data, based on required format for MACH, see MACH help page
# Create MACH directory and folder for each chromosome within working directory
dir.create("MACH")
setwd("./MACH") #navigate to directory
# Create directory for each chromosome
for (c in seq(1,chrom)){
  dir.create(paste0('chrom',c))
}

# Create parent_map, which lists chrom and pos of parent data, population.dat, which lists the overlapping markers between parents and population
parents_map <- parents[,1:3]
parents_map <- parents_map[order(as.numeric(parents_map$snpID)),]

# Subset all_data to keep snpID, Chrom, Pos, and samples
all_data2 <- sim[,-c(4:7)] #Removes cM, F_MISS, REF, ALT
parents_map[,1] <- sub("^","M SNP",parents_map[,1]) #Add "M_SNP" in front of marker ids in population (required format for MACH)
#This is the set of overlapping markers between the parent set and the population
population.dat <- parents_map[which(interaction(parents_map$CHROM, parents_map$POS) %in% 
                                       interaction(all_data2$CHROM,all_data2$POS)),1:2]

# Create .dat files for each chromosome from population.dat
for(c in seq(1,chrom)){
  setwd(paste0(wd,"/MACH/chrom",c))
  write.table(population.dat[population.dat$CHROM == c,][1], file = paste0("progeny_chrom_",c,".dat"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# Create .PED file (genotype data on specific chromosome for population)
all_data2 <- all_data2[with(all_data2,order(as.numeric(snpID))),] #Make sure order of markers is maintained
for (c in seq(1,chrom)){
  population.ped <- as.data.frame(t(all_data2[all_data2$CHROM == c,][,-c(1:3)]))
  population.ped <- data.frame(fID=popID, iID=rownames(population.ped), pID=0, mID=0, sex=0, population.ped) #Required IDs for MACH
  setwd(paste0(wd,"/MACH/chrom",c))
  write.table(population.ped, file = paste0("progeny_chrom_",c,".PED"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# Create parents_chrom_X.snps file (list of parent SNP#'s on specfic chromosome)
for(c in seq(1,chrom)){
  parents.snps <- data.frame(snpID=substring(parents_map[parents_map$CHROM == c,][,1],3))
  setwd(paste0(wd,"/MACH/chrom",c))
  write.table(parents.snps, file = paste0("parents_chrom_",c,".snps"), quote = FALSE,
              sep = "\t", row.names = FALSE,col.names = FALSE)
}

# Create parents_chrom_X.haplos file (parent haplotype data on specific chromosome)
for(c in seq(1,chrom)){
  parents.chrom.haps <- as.data.frame(t(parents[parents$CHROM == c,c(-1:-7)])) #transpose data and remove snpID, CHROM, and POS columns
  
  rownames(parents.chrom.haps) <- paste0(' ',rownames(parents.chrom.haps),' (',c(1:length(8:ncol(parents))),')') #Set row names as parent IDs
  parents.chrom.haps <- as.data.frame(parents.chrom.haps[rep(1:nrow(parents.chrom.haps), each=2),]) #Creates two rows for each parent (each homolog)
  
  parents.chrom.haps <- data.frame("fID" = rownames(parents.chrom.haps),parents.chrom.haps, row.names = NULL)
  parents.chrom.haps <- cbind.data.frame('index' = rownames(parents.chrom.haps),
                                          "fID" =gsub('\\.1',"", parents.chrom.haps$fID),"HapID" = c("HAP1","HAP2"),
                                          parents.chrom.haps[2:ncol(parents.chrom.haps)]) #creates and formats dataframe
  
  # Removing the second allele from first haplotype and first allele from second haplotype
  parents.chrom.haps <- rbind.data.frame(lapply(parents.chrom.haps[c(TRUE,FALSE),], gsub, pattern='/[A-Z]', replacement=''),
                                          lapply(parents.chrom.haps[c(FALSE,TRUE),], gsub, pattern='[A-Z]/', replacement=''))
  # Reordering into correct order
  parents.chrom.haps <- parents.chrom.haps[order(as.numeric(parents.chrom.haps$index)),]
  setwd(paste0(wd,"/MACH/chrom",c))
  write.table(parents.chrom.haps[,-1],file = paste0('parents_chrom_',c,".haplos"),quote = FALSE, sep = '\t',
              row.names = FALSE,col.names = FALSE)
}
