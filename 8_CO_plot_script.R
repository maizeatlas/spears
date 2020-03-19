##Script for plotting crossover count comparisons with histogram in margin

# Required Packages
library("ggplot2")
library("ggExtra")

# User inputs

# Working directory
wd <- "~/Dropbox/Maize_ATLAS_share/ParallelSelection/GBS/Manuscripts/RABBIT_Bio_App/Data/"
# Name of crossover count data
co <- "CO_counts_known_OVD_03MAR20.csv"

# Set working directory
setwd(wd)

# Import crossover count data (from script #7)  
all_CO <- read.csv(co, head=T)

# Plot

tiff("../For_submission/Figure_S5.tiff", width=6 , height=6, units="in", compression="none", res=600)

p <- ggplot(all_CO) +
  geom_point( 
    aes(x=known_CO, y=OVD_CO)) +
  xlim(175,325) + 
  ylim(175,325) +
  geom_hline(yintercept = mean(all_CO$OVD_CO), color="red") +
  geom_vline(xintercept = mean(all_CO$known_CO), color="red") +
  xlab("Known") + ylab("Inferred by RABBIT") + 
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
ggMarginal(p, type = "histogram")

dev.off()





