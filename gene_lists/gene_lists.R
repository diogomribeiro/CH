
library(data.table)
#data = fread("/home/dribeiro/git/compound_het/results/lof.nosingle.results")
#data = fread("/home/dribeiro/git/compound_het/results/missense.nosingle.results")
data = fread("/home/dribeiro/git/compound_het/results/synonymous.nosingle.results")

mgi_essential_data = fread("/home/dribeiro/git/compound_het/gene_lists/mgi_essential", header = F)
gnomad_essential_data = fread("/home/dribeiro/git/compound_het/gene_lists/gnomad_essential", header = F)
ADaM_essential_data = fread("/home/dribeiro/git/compound_het/gene_lists/ADaM_essential", header = F)
hart2014_essential_data = fread("/home/dribeiro/git/compound_het/gene_lists/hart2014_essential", header = F)
hart2017_essential_data = fread("/home/dribeiro/git/compound_het/gene_lists/hart2017_essential", header = F)

hart2014_nonessential_data = fread("/home/dribeiro/git/compound_het/gene_lists/hart2014_nonessential", header = F)
gnomad_tolerant_data = fread("/home/dribeiro/git/compound_het/gene_lists/gnomad_homozygous_tolerant", header = F)
gnomad_nonessential_data = fread("/home/dribeiro/git/compound_het/gene_lists/hart2014_nonessential", header = F)

ch_genes = data

ch_genes$mgi_essential = 0
ch_genes[geneName %in% mgi_essential_data$V1]$mgi_essential = 1

ch_genes$gnomad_essential = 0
ch_genes[geneName %in% gnomad_essential_data$V1]$gnomad_essential = 1

ch_genes$adam_essential = 0
ch_genes[geneName %in% ADaM_essential_data$V1]$adam_essential = 1

ch_genes$hart2014_essential = 0
ch_genes[geneName %in% hart2014_essential_data$V1]$hart2014_essential = 1

ch_genes$hart2017_essential = 0
ch_genes[geneName %in% hart2017_essential_data$V1]$hart2017_essential = 1

ch_genes$hart2014_nonessential = 0
ch_genes[geneName %in% hart2014_nonessential_data$V1]$hart2014_nonessential = 1

ch_genes$gnomad_tolerant = 0
ch_genes[geneName %in% gnomad_tolerant_data$V1]$gnomad_tolerant = 1

ch_genes$gnomad_nonessential = 0
ch_genes[geneName %in% gnomad_nonessential_data$V1]$gnomad_nonessential = 1

ch_genes[numCH > 0][mgi_essential == 1 | gnomad_essential == 1 | adam_essential == 1 | hart2014_essential == 1 | hart2017_essential == 1]


write.table(ch_genes[mgi_essential == 1 | gnomad_essential == 1 | adam_essential == 1 | hart2014_essential == 1 | hart2017_essential == 1]$geneName,"~/git/compound_het/gene_lists/any_essential.txt",quote=F,row.names=F,col.names=F)

ch_genes

#write.table(ch_genes,"~/git/compound_het/gene_lists/lof_gene_lists.tsv",quote=F,row.names=F,sep="\t")
#write.table(ch_genes,"~/git/compound_het/gene_lists/missense_gene_lists.tsv",quote=F,row.names=F,sep="\t")
write.table(ch_genes,"~/git/compound_het/gene_lists/synonymous_gene_lists.tsv",quote=F,row.names=F,sep="\t")

