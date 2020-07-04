##Karyogram plot for a single sample from RABBIT OVD output

#Required Packages
library("ggplot2")
library("gdata")
library("data.table")
library("tidyr")

# User inputs
# Working directory
wd <- "~/working/directory/"
#Sample to plot (format "sim.sample#")
s_id <- "sim.sample5564"
#Color Scheme for Parent haplotype blocks (needs to be length of total number of founders)
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Set working directory
setwd(wd)

#User inputs (read from file created in 2_SAEGUS_to_MACH_format.R
user <- read.table("user_input.txt", header = F, sep = "\t", stringsAsFactors = F)

#Chromosome coordinates
agpv4Ideogram <- read.table("chrom_coord.csv",head=T,sep=",",stringsAsFactors = F)

# Number of chromosomes
chrom <- as.numeric(user[user[[1]]=="chrom",2])
#Sample Number
sn <- as.numeric(user[user[[1]]=="sn",2])

#Founder Data
key <- read.table(user[user[[1]]=="fd",2], head=T, stringsAsFactors = FALSE, sep="\t")
key <- rename.vars(key,"chr","CHROM")
founders <- names(key)[8:ncol(key)]

## Make ideogram: with p and q arms (sets chromosome framework for plot)
agpv4Ideogram2 <- rbind(agpv4Ideogram,agpv4Ideogram)
agpv4Ideogram2$homolog <- c(rep("a",chrom),rep("b",chrom))
agpv4Ideogram2$homolog_start <- c(rep(3,chrom),rep(0,chrom))
agpv4Ideogram2$homolog_stop <- c(rep(5,chrom),rep(2,chrom))

####
####Karyogram from RABBIT Output###
####

# Format RABBIT OVD Output
fn <- length(names(key))-7 # Number of founders
gtnum <- sum(seq(1,fn,1)) # Number of possible genotype combinations

# Diplotype key from RABBIT Output
d.key <- read.table(paste("./RABBIT/reconstructed_chrom_1_jointModel_OVD_summary.csv", sep=""),head=T,
                    sep=",",skip=(9 + sn + 2*fn),nrows = (fn*fn), stringsAsFactors = F)
d.key$Diplotype <- as.numeric(sub("diplotype","", d.key$Diplotype))

## SNP index information for each chromosome
SNPkey <- list()
for (i in 1:chrom){
  SNPkey[[i]] <- read.table(paste("./RABBIT/reconstructed_chrom_",i,"_jointModel_OVD_summary.csv",sep=""),head=F,
                            sep=",", skip=1, nrows=2, stringsAsFactors = F)
  SNPkey[[i]] <- as.data.frame(t(SNPkey[[i]][,2:ncol(SNPkey[[i]])]))
  SNPkey[[i]] <- rename.vars(SNPkey[[i]],c("V1","V2"),c("SNP","CHROM"))
  SNPkey[[i]]$index <- seq(1,nrow(SNPkey[[i]]),1)
}


### Format Viterbi Paths from RABBIT output for calculating AAA,AAG

datalist3 <- list()
for (c in 1:chrom){
  
    vit_path <- read.table(paste("./RABBIT/reconstructed_chrom_",c,"_jointModel_OVD_summary.csv",sep=""),head=T,
                           sep=",", skip=(11 + sn + 2*fn + (fn*fn)), nrows=sn, stringsAsFactors = F)
    dat <- subset(vit_path,vit_path==s_id)
    #dat <- as.data.frame(unlist(strsplit(vit_path[i,2],split = "-")))
    dat <- as.data.frame(unlist(strsplit(dat[,2],split = "-")))
    even <- seq(2,nrow(dat)-1,2)
    odd <- seq(1,nrow(dat)-1,2)
    dat <- data.frame(index_start=(dat$`unlist(strsplit(dat[, 2], split = "-"))`[odd]), 
                      diplotype=(dat$`unlist(strsplit(dat[, 2], split = "-"))`[even]))
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
    #dat$sample <- paste("sim.sample",i,sep="")
    dat$sample <- s_id
    datalist3[[c]] <- dat
}

df_RABBIT <- do.call("rbind",datalist3)
#Need physical positions
df_RABBIT <- merge(df_RABBIT, key[,c(1,3)], by.x=c("SNP_start"), by.y=c("snpID"), all.x=TRUE, all.y=FALSE)
df_RABBIT <- merge(df_RABBIT, key[,c(1,3)], by.x=c("SNP_stop"), by.y=c("snpID"), all.x=TRUE, all.y=FALSE)
#Split parent diplotypes
df_RABBIT <- separate(df_RABBIT, founder, c("a", "b"), sep="\\|")
#Reformat for plotting
df_RABBIT <- melt(as.data.table(df_RABBIT),id=c(1,2,9:12),measure=7:8)
df_RABBIT <- rename.vars(df_RABBIT,c("variable","value"), c("homolog","founder"))
#Reorder for plotting
df_RABBIT <- df_RABBIT[order(df_RABBIT$CHROM,df_RABBIT$POS.x,df_RABBIT$homolog),]

## a temporary solution ... filling gaps
df_RABBIT <- df_RABBIT[order(df_RABBIT$CHROM,df_RABBIT$homolog,df_RABBIT$POS.x),]
t0 <- (df_RABBIT$POS.x[2:nrow(df_RABBIT)] - df_RABBIT$POS.y[1:nrow(df_RABBIT)-1]) / 2
t0[t0<0] <- 0
df_RABBIT$POS.x[2:nrow(df_RABBIT)] <- df_RABBIT$POS.x[2:nrow(df_RABBIT)] - t0
df_RABBIT$POS.y[1:nrow(df_RABBIT)-1] <- df_RABBIT$POS.y[1:nrow(df_RABBIT)-1] + t0

## define continuous value for y-width and spacing
df_RABBIT <- df_RABBIT[order(df_RABBIT$CHROM,df_RABBIT$POS.x,df_RABBIT$homolog),]
df_RABBIT$homolog_start <- c(rep(c(3,0),(nrow(df_RABBIT)/2)))
df_RABBIT$homolog_stop <- c(rep(c(5,2),(nrow(df_RABBIT)/2)))

#Rename for plotting
colnames(df_RABBIT)[which(colnames(df_RABBIT)=="POS.x")] <- "Start"
colnames(df_RABBIT)[which(colnames(df_RABBIT)=="POS.y")] <- "End"

#Merge RABBIT output with ideogram information
df_RABBIT <- merge(df_RABBIT, 
                   agpv4Ideogram2, 
                   by=c("CHROM","homolog","Start","End","homolog_start","homolog_stop"), 
                   all = TRUE)


#Plot
tiff(paste0("../RABBIT_karyogram_",s_id,".tiff"), width=8.5 , height=11, units="in", compression="none", res=600)
ggplot(df_RABBIT) +
  geom_rect(mapping=aes(ymin=homolog_start, ymax=homolog_stop, xmin=Start, xmax=End, fill=founder)) +
  scale_fill_manual(values = cbbPalette, breaks=founders) +
  geom_rect(data=agpv4Ideogram2, mapping=aes(ymin=homolog_start, ymax=homolog_stop, xmin=Start, xmax=End), fill=NA, colour="black") +
  geom_point(data=agpv4Ideogram2, 
             aes(x=rowMeans(cbind(centromer_start,centromer_stop)), y=rowMeans(cbind(homolog_stop,homolog_start))), 
             pch=21, size=4, stroke=1, fill="grey95", colour="black") +
  facet_grid(CHROM ~ ., switch="y") +
  xlab("") + ylab("") +
  #scale_x_continuous(breaks=c(0,5e7,1e8,1.5e8,2e8,2.5e8,3e8),
  #                   labels=c("0 Mb", "50 Mb", "100 Mb", "150 Mb", "200 Mb", "250 Mb", "300 Mb")) +
  
  theme(plot.background = element_blank()
        #,legend.position = "none"
        ,legend.direction="vertical"
        ,legend.key.size = unit(0.4, "cm")
        ,legend.title = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.grid.major.x = element_blank()
        ,panel.grid.major.y = element_blank()
        ,panel.spacing = unit(1, "lines")
        ,panel.border = element_blank()
        ,panel.background = element_blank()
        ,plot.title = element_text(face="bold")
        ,axis.title.x = element_text(size=14, margin = margin(t=10))
        ,axis.text.x = element_text(size=10, margin = margin(t=2.5), vjust=0.5, hjust=1, angle=90)
        ,axis.text.y = element_blank()
        ,axis.title.y = element_text(size=14, vjust=2)
        ,axis.ticks.y = element_blank()
        ,plot.margin=unit(c(0.25,0.25,1,0.5), "cm")
  )

dev.off()
