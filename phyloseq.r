library("phyloseq")
library("scales")
library("grid")
library("ggplot2")

##---------------------------------------------------SETUP PATHS

base_dir <- "~/iihg/metagenomics/just_musselbed_2016/"

biom_file <- paste(base_dir,"open_ref_otus/otu_table_mc2_w_tax_no_pynast_failures_json.biom", sep='')
map_file <- paste(base_dir, "mapping.txt", sep="")
tree_file <- paste(base_dir, "open_ref_otus/rep_set.tre", sep="")

# RUN THIS AT CMD LINE BEFORE IMPORTING BIOM FROM QIIME 1.9
#biom convert -i otu_table_mc2_w_tax_no_pynast_failures.biom -o otu_table_mc2_w_tax_no_pynast_failures_json.biom --table-type="OUT table" --to-json

##--------------------------------------------------IMPORT AND MERGE

biom_otu_tax <- import_biom(biom_file)
map <- import_qiime_sample_data(map_file)
tree <- read_tree(tree_file)

river_bed <- merge_phyloseq(biom_otu_tax, map, tree)
colnames(tax_table(river_bed)) <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")

##-------------------------------------------------DATA CHECKING

any(taxa_sums(river_bed) == 0)   #empty OTUs?
any(sample_sums(river_bed) == 0)  #empty Samples?

readsumsdf = data.frame(
  nreads = sort(taxa_sums(river_bed), TRUE),
  sorted = 1:ntaxa(river_bed),
  type = "OTUs"
)
readsumsdf = rbind(readsumsdf, data.frame(
  nreads = sort(sample_sums(river_bed), TRUE),
  sorted = 1:nsamples(river_bed),
  type="Samples")
)

p = ggplot(readsumsdf, aes(x=sorted, y=nreads)) + geom_bar(stat="identity") + facet_wrap(~type, scales='free')
p + scale_y_log10()

##------------------------------------------------ALPHA DIVERSITY BOXPLOTS (done before filtering out rare OTUs)

plot_richness(river_bed, x="mussel_bed", color="mussel_bed", measure=c("Observed","Chao1","Shannon")) + geom_boxplot()
rb3 = subset_samples(river_bed, depth=="THREE")
rb5 = subset_samples(river_bed, depth=="FIVE")
##------------------------------------------------ALPHA DIV STATISTICAL TESTS

rich = estimate_richness(river_bed)

obs_mussel = rich$Observed[river_bed@sam_data$mussel_bed == 'Mussel']
obs_nomussel = rich$Observed[river_bed@sam_data$mussel_bed == 'No Mussel']

t_obs = t.test(obs_mussel,obs_nomussel)

chao1_mussel = rich$Chao1[river_bed@sam_data$mussel_bed == 'Mussel']
chao1_nomussel = rich$Chao1[river_bed@sam_data$mussel_bed == 'No Mussel']

t_chao1 = t.test(chao1_mussel,chao1_nomussel)

shannon_mussel = rich$Shannon[river_bed@sam_data$mussel_bed == 'Mussel']
shannon_nomussel = rich$Shannon[river_bed@sam_data$mussel_bed == 'No Mussel']

t_shannon = t.test(shannon_mussel,shannon_nomussel)

##---------------------------------------------ORDINATION 

#----all samples depths
library(vegan)
rb100 <- prune_taxa(names(sort(taxa_sums(river_bed),TRUE)[1:100]), river_bed)
rb100_nmds_bray = ordinate(physeq = rb100, method = 'NMDS', distance = 'bray')
plot_ordination(rb100, rb100_nmds_bray, type='sample',color='mussel_bed')

stressplot(rb100_nmds_bray)  ##does the stress converge and fit the model well?
ordihull(rb100_nmds_bray, groups=rb100@sam_data$mussel_bed, col='grey90')  ## hulls by condition

#----3CM depth only
library(vegan)
rb3_100 <- prune_taxa(names(sort(taxa_sums(rb3),TRUE)[1:100]), rb3)
rb3_100_nmds_bray = ordinate(physeq = rb3_100, method = 'NMDS', distance = 'bray')
plot_ordination(rb3_100, rb3_100_nmds_bray, type='sample',color='mussel_bed')

stressplot(rb100_nmds_bray)  ##does the stress converge and fit the model well?

#-----5CM depth only
library(vegan)
rb5_100 <- prune_taxa(names(sort(taxa_sums(rb5),TRUE)[1:100]), rb5)
rb5_100_nmds_bray = ordinate(physeq = rb5_100, method = 'NMDS', distance = 'bray')
plot_ordination(rb5_100, rb5_100_nmds_bray, type='sample',color='mussel_bed')

stressplot(rb5_100_nmds_bray)  ##does the stress converge and fit the model well?

##--------------------------------------------BAR PLOTS

rb_50 <- prune_taxa(names(sort(taxa_sums(river_bed),TRUE)[1:50]), river_bed)
rbm_50 <- merge_samples(rb_50, "mussel_bed")
sample_data(rbm_50)$mussel_bed <- levels(sample_data(rb_50)$mussel_bed)

rbm_50 <- transform_sample_counts(rbm_50, function(x) 100*x / sum(x))
title <- "Differences in Most Abundant Phyla in Mussel-populated Sediment versus Control Sediment"
p = plot_bar(rbm_50, x = "mussel_bed", fill = "Family", facet_grid = ~Phylum, title=title) + ylab("% Abundance (Sequences)") + ylim(0,25) + xlab("Sample Type")
p + geom_bar(aes(), stat="identity", position="stack")

bb <- psmelt(rbm_50)  ## capture the dataframe from plotting for export

##---------------------------------------------TREE PLOTS

plot_tree(rb_50, color="mussel_bed", label.tips="Order", size="abundance", ladderize=TRUE)

##---------------------------------------------HEATMAPS

#TOP 50 MOST ABUNDANT TAXA
rb_pruned <- prune_taxa(names(sort(taxa_sums(river_bed),TRUE)[1:50]), river_bed)
plot_heatmap(rb_pruned, "NMDS", "bray", "mussel_bed", "Phylum")

##---------------------------------------------DESEQ2 ANALYSIS

library("DESeq2")

mussel_dds = phyloseq_to_deseq2(river_bed, ~ mussel_bed)
mussel_dds = DESeq(mussel_dds, test="Wald", fitType="parametric")

res = results(mussel_dds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(river_bed)[rownames(sigtab), ], "matrix"))
sigtab = sigtab[order(sigtab$padj), ]

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggplot(sigtab, aes(x=Phylum, y=log2FoldChange, color=padj)) + geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

