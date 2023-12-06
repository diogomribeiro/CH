# 29-Sep-2022 Diogo Ribeiro @ UNIL
# Perform direct comparison OR

library(data.table)
require(ggplot2)
options("scipen"=100, "digits"=2)

inputFile = "data/lof_3_cats.out"

data <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
colnames(data) = c("ExternalList","Size1","ItemGroup","Size2","Overlap","OddsRatio","Pvalue")

MAX_X_AXIS = 6

data[ExternalList == "hart2014_essential"]$ExternalList = "Essential in culture (Hart 2014)"
data[ExternalList == "hart2017_essential"]$ExternalList = "Essential CRISPR (Hart 2017)"
data[ExternalList == "gnomad_essential"]$ExternalList = "Essential GnomAD"
data[ExternalList == "gnomad_nonessential"]$ExternalList = "Non-essential GnomAD"
data[ExternalList == "gnomad_homozygous_tolerant"]$ExternalList = "Homozygous LoF tolerant"
data[ExternalList == "mgi_essential"]$ExternalList = "Essential in mice (MGI)"
data[ExternalList == "hart2014_nonessential"]$ExternalList = "Non-essential CRISPR (Hart 2014)"
data[ExternalList == "ADaM_essential"]$ExternalList = "ADaM essential"

o = c("Homozygous LoF tolerant","Non-essential GnomAD","Non-essential CRISPR (Hart 2014)","Essential in culture (Hart 2014)","Essential CRISPR (Hart 2017)", "Essential in mice (MGI)","Essential GnomAD","ADaM essential")
data$ExternalList = factor(as.character(data$ExternalList), levels = o)
data = data[order(-data$ExternalList)]


direct_fisher <- function(data, data1, data2) {
  dataset = data.table()
  for (extList in unique(data$ExternalList)){
    d1 = data[ExternalList == extList][ItemGroup == data1]
    d2 = data[ExternalList == extList][ItemGroup == data2]
    TP = d1$Overlap
    FP = d1$Size2 - d1$Overlap
    FN = d2$Overlap
    TN = d2$Size2 - d2$Overlap
    m = matrix(c(TP,FP,FN,TN),nrow=2)
    f = fisher.test(m, conf.level = 0.95)
    dataset = rbind(dataset,  data.table(ExternalList = extList, odds = f$estimate, pval = f$p.value, confmin = f$conf.int[1], confmax = f$conf.int[2]) )
  }
  
  dataset$odds = log2(dataset$odds)
  dataset$confmin = log2(dataset$confmin)
  dataset$confmax = log2(dataset$confmax)
  
  dataset[odds < -MAX_X_AXIS]$odds = -MAX_X_AXIS
  dataset[odds > MAX_X_AXIS]$odds = MAX_X_AXIS
  dataset[confmin < -MAX_X_AXIS]$confmin = -MAX_X_AXIS
  dataset[confmax > MAX_X_AXIS]$confmax = MAX_X_AXIS
  dataset$idx = seq(1,nrow(dataset)*3,3)
  return(dataset)
}


data1 = "lof_ch_genes"
data2 = "lof_possible_no_ch_genes"
dataset1 = direct_fisher(data,data1,data2)

data1 = "lof_high_ch_genes"
data2 = "lof_high_possible_no_ch_genes"
dataset2 = direct_fisher(data,data1,data2)

data1 = "rand_ch_genes"
data2 = "rand_possible_no_ch_genes"
dataset3 = direct_fisher(data,data1,data2)

labels = c("Full data","High confidence","Random phasing")
colors = c("#c6dbef","#2171b5","#fc9272")


#############
# Essential genes
#############

png("main_b_essential.png",width = 900, height = 875)

f1 = dataset1[c(1,2,3,4,5),]
f2 = dataset2[c(1,2,3,4,5),]
f3 = dataset3[c(1,2,3,4,5),]

g1 = ggplot() +
  geom_rect(data = f1[idx %% 2 == 0], aes(xmin = -Inf, xmax = Inf, ymin = idx-1.5, ymax = idx+1.5), fill = "#d9d9d9", size = 1) +
  geom_vline(xintercept = 0, linetype = 2, size = 1.5) +
  geom_segment(data = f1, aes(x = confmin, xend = confmax, y = idx+0.8, yend = idx+0.8 ), color = "black", size = 9 ) + #, alpha = -log10(pval)
  geom_segment(data = f1, aes(x = confmin, xend = confmax, y = idx+0.8, yend = idx+0.8 ), color = colors[1], size = 7 ) + #, alpha = -log10(pval)
  geom_point(data = f1, aes(x = odds, y = idx+0.8 ), fill = colors[1], size = 13.5, shape = 21) +  #, alpha = -log10(pval)
  geom_segment(data = f2, aes(x = confmin, xend = confmax, y = idx+0, yend = idx+0 ), color = "black", size = 9 ) + #alpha = -log10(pval) 
  geom_segment(data = f2, aes(x = confmin, xend = confmax, y = idx+0, yend = idx+0 ), color = colors[2], size = 7 ) + #alpha = -log10(pval) 
  geom_point(data = f2, aes(x = odds, y = idx+0), fill = colors[2], size = 13.5, shape = 21) + # alpha = -log10(pval) 
  geom_segment(data = f3, aes(x = confmin, xend = confmax, y = idx-0.8, yend = idx-0.8 ), color = "black", size = 9 ) + #, alpha = -log10(pval) 
  geom_segment(data = f3, aes(x = confmin, xend = confmax, y = idx-0.8, yend = idx-0.8 ), color = colors[3], size = 7 ) + #, alpha = -log10(pval) 
  geom_point(data = f3, aes(x = odds, y = idx-0.8), fill = colors[3], size = 13.5, shape = 21) + #, alpha = -log10(pval) 
  xlab("log2(Odds ratio)") +
  ylab("") +
  scale_x_continuous(limits = c(-MAX_X_AXIS,MAX_X_AXIS)) +
  scale_y_continuous(breaks = seq(1,nrow(f1)*3,3), labels = f1$ExternalList) + #dataset$ExternalList
  theme_minimal() + 
  theme(text = element_text(size=30), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2.5))

g1
dev.off()

#############
# Non-Essential genes
#############

png("main_b_nonessential.png",width = 900, height = 550)

d1 = dataset1[c(6,7,8),]
d2 = dataset2[c(6,7,8),]
d3 = dataset3[c(6,7,8),]

g2 = ggplot() +
  geom_rect(data = d1[idx == 19], aes(xmin = -Inf, xmax = Inf, ymin = idx-1.5, ymax = idx+1.5), fill = "#d9d9d9", size = 1) +
  geom_vline(xintercept = 0, linetype = 2, size = 1.5) +
  geom_segment(data = d1, aes(x = confmin, xend = confmax, y = idx+0.8, yend = idx+0.8 ), color = "black", size = 9 ) + #, alpha = -log10(pval)
  geom_segment(data = d1, aes(x = confmin, xend = confmax, y = idx+0.8, yend = idx+0.8 ), color = colors[1], size = 7 ) + #, alpha = -log10(pval)
  geom_point(data = d1, aes(x = odds, y = idx+0.8 ), fill = colors[1], size = 13.5, shape = 21) +  #, alpha = -log10(pval)
  geom_segment(data = d2, aes(x = confmin, xend = confmax, y = idx+0, yend = idx+0 ), color = "black", size = 9 ) + #alpha = -log10(pval) 
  geom_segment(data = d2, aes(x = confmin, xend = confmax, y = idx+0, yend = idx+0 ), color = colors[2], size = 7 ) + #alpha = -log10(pval) 
  geom_point(data = d2, aes(x = odds, y = idx+0), fill = colors[2], size = 13.5, shape = 21) + # alpha = -log10(pval) 
  geom_segment(data = d3, aes(x = confmin, xend = confmax, y = idx-0.8, yend = idx-0.8 ), color = "black", size = 9 ) + #, alpha = -log10(pval) 
  geom_segment(data = d3, aes(x = confmin, xend = confmax, y = idx-0.8, yend = idx-0.8 ), color = colors[3], size = 7 ) + #, alpha = -log10(pval) 
  geom_point(data = d3, aes(x = odds, y = idx-0.8), fill = colors[3], size = 13.5, shape = 21) + #, alpha = -log10(pval) 
  xlab("log2(Odds ratio)") +
  ylab("") +
  guides(fill = guide_legend(reverse = TRUE)) + 
  scale_x_continuous(limits = c(-2,2)) +
  scale_y_continuous(breaks = c(16,19,22), labels = d1$ExternalList) + #dataset$ExternalList
  theme_minimal() + 
  theme(text = element_text(size=30), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2.5))

g2
dev.off()
