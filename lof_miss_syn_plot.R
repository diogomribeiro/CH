#03-Oct-2022 Diogo Ribeiro @ UNIL
# CH results summary

library(data.table)
library(ggplot2)

dataLoF = fread("data/lof.nosingle.results")
dataSyn = fread("data/synonymous.nosingle.results")
dataMis = fread("data/missense.nosingle.results")

## Basic plot
nTotalGenes = nrow(dataSyn)
nLoFCH = nrow(dataLoF[numCH > 0])
nLoFPossible = nrow(dataLoF[numCH == 0][num2mut > 0])
nLoFOther = nTotalGenes - nLoFPossible - nLoFCH
nSynCH = nrow(dataSyn[numCH > 0])
nSynPossible = nrow(dataSyn[numCH == 0][num2mut > 0])
nSynOther = nTotalGenes - nSynPossible - nSynCH
nMisCH = nrow(dataMis[numCH > 0])
nMisPossible = nrow(dataMis[numCH == 0][num2mut > 0])
nMisOther = nTotalGenes - nMisPossible - nMisCH

nLoFCH + nLoFPossible + nLoFOther == nTotalGenes
nSynCH + nSynPossible + nSynOther == nTotalGenes
nMisCH + nMisPossible + nMisOther == nTotalGenes

sum(dataLoF$numCH)
sum(dataSyn$numCH)
sum(dataMis$numCH)

png("sup_basic_mis_syn.png",width = 1000, height = 900)

basic = data.table(labels = c("LoF","Syn","Mis"), CH = c(nLoFCH,nSynCH,nMisCH))
g1 = ggplot(basic, aes(x = labels, y = CH, fill = labels)) +
  geom_bar(stat = "identity", size = 1, color = "black", width = 0.7) +
  geom_text(aes(label = CH), size = 10, vjust = -1, fontface = "bold") +
  theme_minimal() +
  ylab("# Genes with CH") +
  xlab("Variant type") +
  scale_x_discrete(limits = c("LoF","Mis","Syn"),labels = c("LoF","Missense","Synonymous")) +
  scale_fill_manual(values = c("#737373","#bdbdbd","#f0f0f0")) +
  ylim(c(0,max(basic$CH) + 1000)) +
  theme(text = element_text(size = 38), plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        # panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 2.5),  ) #aspect.ratio = 1

g1
dev.off()


#####################
# Proportion CH / assessed
#####################

png("sup_prop_ch_mis_syn.png",width = 1000, height = 900)

dt = data.table(group = c("LoF","LoF","LoF","Syn","Syn","Syn","Miss","Miss","Miss"), 
                Gene = c("3: ch","2: possible","1: other","3: ch","2: possible","1: other","3: ch","2: possible","1: other"), 
                values = c(nLoFCH,nLoFPossible,nLoFOther,nSynCH,nSynPossible,nSynOther,nMisCH,nMisPossible,nMisOther))

g2 = ggplot(dt,aes(x = group, fill = Gene, y = values)) + 
  geom_bar(stat = "identity", size = 2, color = "black", width = 0.6) +
  theme_minimal() +
  scale_fill_manual(values = c("#ccebc5","#fbb4ae","#c6dbef"), labels = c("Not assessed","Not Compound","Compound")) +
  scale_x_discrete(limits = c("Syn","Miss","LoF"), labels = c("Synonymous","Missense","LoF")) + 
  ylab("# Genes") +
  xlab("Variant type") +
  theme(text = element_text(size = 30), plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 2.5), aspect.ratio = 1  )

g2
dev.off()

#####################
# Observed vs Expected
#####################
png("main_d.png",width = 1000, height = 800)

dataLoF$LoF = dataLoF$numCH/dataLoF$expectedCH
dataSyn$Syn = dataSyn$numCH/dataSyn$expectedCH
dataMis$Mis = dataMis$numCH/dataMis$expectedCH

m1 = merge(dataLoF[,.(gene,LoF)],dataSyn[,.(gene,Syn)], by = "gene")
mergedData = merge(m1, dataMis[,.(gene,Mis)], by = "gene")
meltedData = melt(mergedData)
summary(mergedData$LoF)
summary(mergedData$Syn)
summary(mergedData$Mis)

median(meltedData[variable == "LoF"][!is.na(value)]$value)

g3 = ggplot(meltedData,aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot(width = 0.3, size = 1) +
  theme_minimal() +
  xlab("Variant type") +
  ylab("Observed CH / Expected CH") +
  scale_x_discrete( limits = c("LoF","Mis","Syn"), labels = c("LoF","Missense","Synonymous")) +
  scale_fill_manual(values = c("#737373","#f0f0f0","#bdbdbd")) +
  theme(text = element_text(size = 36), plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none",
        panel.background = element_rect(colour = "black", fill = "white", size = 2.5)  ) # aspect.ratio = 1

g3
dev.off()


t = wilcox.test(meltedData[variable == "LoF"][!is.na(value)]$value,meltedData[variable == "Mis"][!is.na(value)]$value)
t$p.value

mean(meltedData[variable == "LoF"][!is.na(value)]$value)
mean(meltedData[variable == "Mis"][!is.na(value)]$value)
mean(meltedData[variable == "Syn"][!is.na(value)]$value)

t = wilcox.test(meltedData[variable == "Syn"][!is.na(value)]$value,meltedData[variable == "Mis"][!is.na(value)]$value)
t$p.value