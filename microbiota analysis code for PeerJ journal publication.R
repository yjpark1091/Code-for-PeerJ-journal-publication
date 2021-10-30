library(phyloseq)
library(qiime2R)
library(microbiome)
library(microbiomeSeq)
library(metagMisc)
library(ggplot2)
library(ggpubr)
library(DESeq2)
library(vegan)
library(knitr)
library(dplyr)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(plyr)
library(reshape2)
library(viridis)
library(igraph)
library(RColorBrewer)


#1. data import
physeq <- qza_to_phyloseq(
  features="rarefied_filtered_table-dada2.qza",
  metadata="sample-metadata_1st_revision.tsv",
  taxonomy="taxonomy.qza",
  tree="tree.qza")

#genes level sort & normalization
physeq_genus <- tax_glom(physeq, "Genus")
physeq_genus <- phyloseq_rm_na_tax(physeq_genus)

transform <- microbiome::transform
physeq_genus <- transform(physeq_genus, "compositional")

#2. Alpha diversity
ps1 <- prune_taxa(taxa_sums(physeq_genus) >= 0, physeq_genus)

tab <- microbiome::alpha(ps1, index = c("shannon", "chao1", "observed", "evenness_pielou"))
kable(head(tab))

ps1.meta <- meta(ps1)
kable(head(ps1.meta))

ps1.meta$Shannon <- tab$diversity_shannon 
ps1.meta$Observed <- tab$observed
ps1.meta$Chao1 <- tab$chao1
ps1.meta$Evenness <- tab$evenness_pielou

mode(ps1.meta$Type)
type <- unique(ps1.meta$Type)
type.pairs <- combn(seq_along(type), 2, simplify = FALSE, FUN = function(i)type[i])
print(type.pairs)

p_shannon <- ggboxplot(ps1.meta, x = "Type", y = "Shannon",
                add = "point", fill = "Type", palette = c())
p_shannon <- p_shannon + stat_compare_means(comparisons = type.pairs)
print(p_shannon)

p_chao <- ggboxplot(ps1.meta, x = "Type", y = "Chao1",
                add = "point", fill = "Type", palette = c())
p_chao <- p_chao + stat_compare_means(comparisons = type.pairs)
print(p_chao)

p_even <- ggboxplot(ps1.meta, x = "Type", y = "Evenness",
                    add = "point", fill = "Type", palette = c())
p_even <- p_even + stat_compare_means(comparisons = type.pairs)
print(p_even)


# 3. beta diversity - PERMANOVA, investigate a top taxa which influence a seperating the group
otu <- abundances(physeq)
meta <- meta(physeq)
permanova_bray <- adonis(t(otu) ~ Type,
                    data = meta, permutations=999, method = "bray")
print(as.data.frame(permanova_bray$aov.tab)["Type", "Pr(>F)"])

dist.uf <- phyloseq::distance(physeq, method = "unifrac")
permanova_uf <- adonis2(dist.uf ~ Type,
                       data = meta, permutations=999)
print(as.data.frame(permanova_uf)["Type", "Pr(>F)"])

dist <- vegdist(t(otu))
a <- anova(betadisper(dist, meta$Type))
print(as.data.frame(anova(betadisper(dist, meta$Type))["Group", "Pr(>F)"]))


ord_PCoA_bray <- ordinate(physeq, "PCoA", "bray")
ord_PCoA_unifrac <- ordinate(physeq, "PCoA", "unifrac")

p <- plot_ordination(physeq, ord_PCoA_bray, color = "Type") +
  geom_point(size = 5)
print(as.data.frame(permanova_bray$aov.tab)["Type", "Pr(>F)"])
p <- p + annotate("text", x = 0.3, y = -0.4, label = "PERMANOVA p-value = 0.001")
print(p + stat_ellipse())

p <- plot_ordination(physeq, ord_PCoA_unifrac, color = "Type") +
  geom_point(size = 5)
print(as.data.frame(permanova_uf)["Type", "Pr(>F)"])
p <- p + annotate("text", x = 0.3, y = -0.4, label = "PERMANOVA p-value = 0.001")
print(p + stat_ellipse())

#4.RandomForest classification
physeq_genus <- taxa_level(physeq, "Genus")
NB_sig <- differential_abundance(physeq_genus, grouping_column = "Type",output_norm="relative", pvalue.threshold=0.05, lfc.threshold=2, filename = "NB_significant")

print(p$MDA)


#5. Co-occurrence network
co_occr <- co_occurence_network(physeq_genus, grouping_column = "Type",
                                rhos = c(-0.75, -0.5, 0.5, 0.75),method = "cor",
                                qval_threshold=0.01, select.condition = "Normal",
                                scale.vertex.size=3,
                                scale.edge.width=15,
                                plotNetwork=T, plotBetweennessEeigenvalue=F)

write.csv(co_occr$net$edgelists,"Normal_co_occr.csv")


#6. Differential abundance, abundance, significance, fold change, prevalence, and AUC score
library(SIAMCAT)
physeq_genus <- taxa_level(physeq, which_level = "Genus")
transform <- microbiome::transform
physeq_genus <- transform(physeq_genus, "compositional")

label <- create.label(meta=sample_data(physeq_genus),
                      label = "Type",
                      case = "Downer", control = "Normal",
                      p.lab = "Downer", n.lab = "Normal")


siamcat <- siamcat(phyloseq = physeq_genus, label = label)

show(siamcat)

siamcat <- filter.features(siamcat,
                           filter.method = 'abundance',
                           cutoff = 0.001)

siamcat <- check.associations(
  siamcat,
  sort.by = 'fc',
  alpha = 0.01,
  mult.corr = "fdr",
  detect.lim = 10 ^-6,
  plot.type = "quantile.box",
  panels = c("fc", "prevalence", "auroc"))
