#29-Sep-2022 Diogo Ribeiro @ UNIL
# CH results summary

library(data.table)
library(ggplot2)

png("main_a.png",width = 1000, height = 900)

data1 = fread("data/lof.nosingle.results")
data2 = fread("data/lof.nosingle.cutoff0.99.results")
data3 = fread("data/lof.nosingle.rand.results")

labels = c("Full data","High confidence","Random phasing")
colors = c("#c6dbef","#2171b5","#fc9272")

#############
# Plot # CH genes per category
#############

## Basic plot
n1CH = nrow(data1[numCH > 0])
n2CH = nrow(data2[numCH > 0])
n3CH = nrow(data3[numCH > 0])

dt = data.table(labels = labels, CH = c(n1CH,n2CH,n3CH), cols = colors)

g1 = ggplot(dt, aes(x = labels, y = CH, fill = cols)) + 
  geom_bar(stat = "identity", size = 2, color = "black", width = 0.6) +
  geom_text(aes(label = CH), size = 9, vjust = 1.5, fontface = "bold") +
  theme_minimal() +
  scale_fill_identity(labels = dt$labels, guide = guide_legend(reverse = TRUE)) + 
  scale_x_discrete(limits = dt$labels) +
  ylab("# Genes with CH") +
  xlab("Dataset") +
  theme(text = element_text(size = 36), plot.title = element_text(hjust = 0.5), 
        legend.position = "none",
        legend.box.background = element_rect(colour = "black"), legend.spacing.y = unit(0, "mm"), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 2.5) ) #aspect.ratio = 1 #axis.text.x = element_text(angle = 20, vjust=0.6)

g1
dev.off()


#############
# Plot # CH events per category
#############

png("supp_ch_events_cat.png",width = 1000, height = 900)

c1CH = sum(data1$numCH)
c2CH = sum(data2$numCH)
c3CH = sum(data3$numCH)
dt2 = data.table(labels = labels, CH = c(c1CH,c2CH,c3CH), cols = colors)

g2 = ggplot(dt2, aes(x = labels, y = CH, fill = cols)) + 
  geom_bar(stat = "identity", size = 2, color = "black", width = 0.6) +
  geom_text(data = dt2[c(1,2,3,4),], aes(label = CH), size = 9, vjust = -0.3, fontface = "bold") +
  geom_text(data = dt2[c(5),], aes(label = CH), size = 9, vjust = 1.3, fontface = "bold") +
  theme_minimal() +
  scale_fill_identity(labels = dt$labels, guide = guide_legend(reverse = TRUE)) + 
  scale_x_discrete(limits = dt$labels) +
  ylab("# CH events") +
  xlab("Dataset") +
  theme(text = element_text(size = 36), plot.title = element_text(hjust = 0.5), 
        legend.position = "none",
        legend.box.background = element_rect(colour = "black"), legend.spacing.y = unit(0, "mm"), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 2.5),  ) #aspect.ratio = 1 # axis.text.x = element_text(angle = 20, vjust=0.6)

g2
dev.off()


## Figure legend
png("legend_3_cats.png",width = 600, height = 600)
legendDT = data.table(cols = colors, data = labels, x = c(1,2,3), y = c(1,2,3))
g3 = ggplot(legendDT, aes(x = x, y = y, fill = data)) +
  geom_point(size = 10, shape = 21, color = "black") +
  scale_fill_manual(values = legendDT$cols, labels = labels) +
  theme_minimal() +
  theme(text = element_text(size= 30), legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"), legend.spacing.y = unit(0, "mm"))

g3
dev.off()
