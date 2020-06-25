### For plotting Genome-wide AAA vs Parent Certainty

### Required Packages
library("dplyr")
library("tidyr")
library("data.table")
library("gdata")
library("ggplot2")

# User inputs
# Working directory (contains RABBIT folder)
wd <- "~/working/directory/"

# Set working directory
setwd(wd)

#User inputs (read from file created in 2_SAEGUS_to_MACH_format.R
user <- read.table("user_input.txt", header = F, sep = "\t", stringsAsFactors = F)

#Founder Data
key <- read.table(user[user[[1]]=="fd",2], head=T, stringsAsFactors = FALSE, sep="\t")
key <- rename.vars(key,"chr","CHROM")

#Output from mldiff analysis
ml_bymarker <- read.csv("ml_diff_bymarker_16JUN20.csv",head=T, stringsAsFactors = F)

#Compare to Ancestral Assignement Accuracy results from Viterbi output
AAA <- read.csv(user[user[[1]]=="marker",2],head=T,stringsAsFactors = F)
AAA <- rename.vars(AAA, "AAA", "mean_AAA")

for_plot <- merge(AAA, ml_bymarker,sort=F)
for_plot <- merge(for_plot,key[,1:3],sort = F)

# All chromosomes
for_plot_all <- melt(as.data.table(for_plot), id.vars = c("snpID","CHROM","POS"), measure.vars = c("mean_AAA","ml_diff"))
#Centromere and chromosome information for rug plots
agpv4Ideogram <- read.table("chrom_coord.csv",head=T,sep=",",stringsAsFactors = F)
# agpv4Ideogram <- cbind.data.frame(CHROM=(1:10),
#                                   Start=rep(1,10), 
#                                   End=c(307041717, 244442276, 
#                                         235667834, 246994605, 
#                                         223902240, 174033170, 
#                                         182381542, 181122637, 
#                                         159769782, 150982314),
#                                   centromer_start=c(136.77e6,95.51e6, 
#                                                     85.78e6,109.07e6, 
#                                                     104.54e6,52.3e6,
#                                                     56.38e6,53.75e6,
#                                                     57.36e6,51.39e6),
#                                   centromer_stop=c(137.12e6,97.49e6,
#                                                    86.39e6,110.5e6,
#                                                    106.82e6,53.11e6,
#                                                    56.68e6,55.39e6,
#                                                    57.76e6,52.78e6))


tiff("./AAA_vs_PC.tiff", width=8 , height=10, units="in", compression="none", res=600)
ggplot(for_plot_all) +
  geom_line(data=for_plot_all, 
            aes(x=POS, y=value, color=variable)) +
  scale_colour_manual(values=c(mean_AAA="black",ml_diff="goldenrod"), labels=c("AAA","Parent Certainty")) +
  #Creates the border around marker density rug (may need to adjust ymin and ymax values)
  geom_rect(data = agpv4Ideogram, mapping=aes(ymin=(min(for_plot_all$value)-.15), ymax=(min(for_plot_all$value)), xmin=Start, xmax=End), fill=NA ,colour="black") + 
  #Creates the marker density rug (may need to adjust y and yend values)
  geom_segment(data=for_plot_all, aes(x=POS,y=(min(for_plot_all$value)-.125), xend=POS, yend=(min(for_plot_all$value)-.025)), alpha=0.1) +
  #Creates circle for centromere position (may need to adjust pch and ymin)
  geom_point(data=agpv4Ideogram, 
             aes(x=rowMeans(cbind(centromer_start,centromer_stop)), y=(min(for_plot_all$value)-.075)), 
             pch=21, size=5, stroke=1, fill="grey60", colour="black") +
  facet_grid(CHROM ~ ., switch="y") +
  xlab("Position") + ylab("Chromosome") +
  scale_x_continuous(breaks=c(0,5e7,1e8,1.5e8,2e8,2.5e8,3e8),
                     labels=c("0 Mb", "50 Mb", "100 Mb", "150 Mb", "200 Mb", "250 Mb", "300 Mb")) +
  theme(plot.background = element_blank()
        #,legend.position="none"
        #,legend.direction="vertical"
        #,legend.key.size = unit(0.4, "cm")
        #,legend.title = element_blank()
        #,legend.background = element_rect(linetype="solid", colour ="black")
        #,legend.text = element_text(size=14)
        ,panel.grid.minor = element_blank()
        ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=2)
        ,panel.grid.major.y = element_line(colour="grey", size = 0.5)
        ,panel.spacing = unit(1, "lines")
        ,panel.border = element_blank()
        ,panel.background = element_blank()
        ,plot.title = element_text(face="bold")
        ,axis.title.x = element_text(size=14, margin = margin(t=10))
        ,axis.text.x = element_text(size=10, margin = margin(t=2.5), vjust=0.5, hjust=1, angle=90)
        ,axis.title.y = element_text(size=14, vjust=2)
        ,plot.margin=unit(c(0.25,0.25,1,0.5), "cm")
  )

dev.off()


