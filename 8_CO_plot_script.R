##Script for plotting crossover count comparisons with histogram in margin

# Required Packages
library("ggplot2")
library("ggExtra")
library("ggpubr") ## used for stat_cor

# User inputs

# Working directory
wd <- "/working/directory/"
# Name of crossover count data from script 7_Calculate_OVD_AAA_GAA_SER_CCC.R
co <- "CO_counts_known_vs_OVD.csv"

# Set working directory
setwd(wd)

# Import crossover count data (from script #7)  
all_CO <- read.csv(co, head=T)

# Plot

tiff("./CO_plot.tiff", width=6 , height=6, units="in", compression="none", res=600)

### Will need to customize axes limits and breaks, currently commented out
p <- ggplot(all_CO, aes(x=known_CO, y=OVD_CO)) +
  geom_point() +
  stat_cor(method = "pearson") +
  #scale_x_continuous(limits=c(175,325), breaks=c(200,250,round(mean(all_CO$known_CO),0),300)) +
  #scale_y_continuous(limits=c(175,325), breaks=c(200,round(mean(all_CO$OVD_CO),0),250,300)) +
  geom_hline(yintercept = round(mean(all_CO$OVD_CO),0), color="grey", size=0.75) +
  geom_vline(xintercept = round(mean(all_CO$known_CO),0), color="grey", size=0.75) +
  xlab("Expected (Simulation)") + ylab("Observed (RABBIT)") + 
  theme(plot.background = element_blank()
        ,legend.direction="vertical"
        ,legend.key.size = unit(0.4, "cm")
        ,panel.grid.minor = element_blank()
        ,panel.grid.major.x = element_line(colour="grey", size = 0.5, linetype=2)
        ,panel.grid.major.y = element_line(colour="grey", size = 0.5, linetype=2)
        ,panel.spacing = unit(1, "lines")
        ,panel.border = element_rect(colour="grey", fill=NA)
        ,panel.background = element_blank()
        ,plot.title = element_text(face="bold")
        ,axis.title.x = element_text(size=14, margin = margin(t=10))
        ,axis.text.x = element_text(size=12, margin = margin(t=2.5), vjust=0.5, hjust=1, angle=90)
        ,axis.text.y = element_text(size=12)
        ,axis.title.y = element_text(size=14, vjust=2)
        ,plot.margin=unit(c(0.25,0.25,1,0.5), "cm"),
        aspect.ratio = 1
  )
ggMarginal(p, type = "histogram", color = "white", binwidth = 2)

dev.off()





