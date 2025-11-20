library(rWSBIM2122)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggthemes)
library(org.Mm.eg.db)
library(GO.db)
library(clusterProfiler)
library(ExploreModelMatrix)
library(DOSE)
library(enrichplot)
install.packages("ggVennDiagram")
library(ggVennDiagram)
### ouverture  des samples
samples <- list.files("data",
                      pattern = "*tsv.gz",
                      full.names = TRUE)



counts <- read_tsv(file = samples[1]) %>%
  select(Geneid, ends_with('.bam'))

for (sample in samples[-1]) {
  tmp <- read_tsv(sample) %>%
    select(Geneid, ends_with('.bam'))
  counts <- full_join(counts, tmp, by = "Geneid")
}

names(counts) <- sub(pattern = ".bam$", '', names(counts)) #sub = substract
names(counts) <- sub(pattern = "../processed_data/bam/", '', names(counts))


### conversion en DEseq
mat <- as.matrix(counts[, 2:13])
rownames(mat) <- counts$Geneid
head(mat)

metadata <- read_tsv("data/metadata.tsv")
metadata
mat[ , metadata$sample]
dds <- DESeqDataSetFromMatrix(countData = mat[ ,metadata$sample],
                              colData = metadata,
                              design = ~1)
dds

## Plot of the count distributions of each sample
as_tibble(assay(dds), rownames = "gene") %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 2:13) %>% 
  ggplot(aes(x = log2(counts+1), fill = sample)) +
  geom_histogram(bins = 20) +
  labs(title= "Counts distribution")+
  theme(plot.title = element_text(size = 16, face = "bold.italic", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey"))+
  facet_wrap(~ sample)


dim(dds)  # 54532 genes

plot_data <- as_tibble(assay(dds), rownames = "gene") %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 2:13)

summary(plot_data$counts)

# The count matrix contains many rows with only zeros.So we will remove
## these rows of the dds object (no statistical test can be done on these
## genes that have no variation in their counts). By removing counts containing
## only zero we will reduce the memory size of dds object which increases the speed
# of analysis.
## We remove zeros

dds <- dds[rowSums(assay(dds)) > 0, ]


dim(dds)

## we reduced the dds object to 25470 genes

plot_data <- as_tibble(assay(dds), rownames = "gene") %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 2:13)

summary(plot_data$counts)

any(plot_data$counts == 0)

as_tibble(assay(dds), rownames = "gene") %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 2:13) %>% 
  filter(counts != 0) %>% 
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  labs(title= "Counts distribution after filtering")+
  theme(plot.title = element_text(size = 16, face = "bold.italic", 
                                  hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', 
                                        colour = "grey"))+
  facet_wrap(~ sample)


### DESEQ + PCA
dds <- DESeq(dds)

rld <- rlogTransformation(dds)

pca_data <- plotPCA(rld, intgroup = c("culture", "genotype"), returnData = TRUE)

percentVar <- round(100 * attr(pca_data, "percentVar"))

custom_colours <- c("hypoxia" = "red3", "normoxia" = "steelblue3")

ggplot(pca_data, 
       aes(x = PC1, y = PC2, colour = culture, shape = genotype)) + 
  geom_point(size = 3,alpha = 0.7)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  labs(title= "PCA analysis")+
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold.italic", hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = 'solid'))+
  scale_color_manual(values = custom_colours)+
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey"))

## PCA analysis helps us to define the most appropriate design to use 
## for the analysis


# Dispersion plot
plotDispEsts(dds)

resultsNames(dds)
## each point represents one gene. The higher the expression of a gene is, the lower
## its dispersion is

dds <- DESeqDataSetFromMatrix(countData = mat[ ,metadata$sample],
                              colData = metadata,
                              design = ~ culture + genotype)

dds <- DESeqDataSetFromMatrix(countData = mat[ ,metadata$sample],
                              colData = metadata,
                              design =  ~ culture * genotype)

# Set the normoxia level as the reference level
dds$culture <- relevel(dds$culture, ref = "normoxia")
# Set the wt level as the reference level
dds$genotype <- relevel(dds$genotype, ref = "wt")

dds <- DESeq(dds)

resultsNames(dds)

### Design avec interaction
## Which genes are not affected in the same way by hypoxia in both genotypes ?
## "culturehypoxia.genotypeKO" : extra effect of KO in hypoxia compared to normoxia
res <- results(dds, name = "culturehypoxia.genotypeKO")
res_tbl <- as_tibble(res, rownames = "Geneid") %>%
  arrange(padj)
head(res_tbl)

# ENSMUSG00000005125 has a negative  log2FC meaning thaht the effect of KO is 
## much lower in hypoxia than in normoxia?? 
as_tibble(counts(dds[res_tbl$Geneid[1:4], ], normalize = TRUE),
          rownames = 'Geneid') %>%
  pivot_longer(names_to = "sample", values_to = "counts", -Geneid) %>%
  left_join(as_tibble(colData(dds))) %>%
  mutate(name = paste0(substr(genotype, 1, 5), '_', culture, '_', 1:3)) %>%
  ggplot(aes(x = name, y = counts, fill = culture)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ Geneid, scales = "free") +
  theme(axis.text.x = element_text(size = 8, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))


### histogramme
hist(res_tbl$pvalue)

##Best genes
best_genes <- res_tbl %>%
  arrange(padj)  %>%
  head(6)
best_genes

#Expression des 6 best genes 
res_tbl %>%
  arrange(padj) %>%
  head(6) %>%
  pull(Geneid) %>%
  map(~ {
    as_tibble(counts(dds[.x, ], normalize = TRUE), rownames = 'Geneid') %>%
      pivot_longer(names_to = "sample", values_to = "counts", -Geneid) %>%  
      left_join(as_tibble(colData(dds))) %>% 
      ggplot(aes(x = sample, y = counts, fill = culture)) +
      geom_bar(stat = 'identity', color = "gray30") +
      facet_wrap(~ genotype, scales = "free", ncol = 2) +
      labs(title = paste("Expression du gène :", .x)) +
      theme(axis.text.x = element_text(size = 7, angle = 90),
            axis.title.x = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 7))
  })

# Filtrer les gènes significatifs une seule fois
filtered_genes <- res_tbl %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

# Compter les gènes upregulated, downregulated et au total
upregulated <- sum(filtered_genes$log2FoldChange > 1)
downregulated <- sum(filtered_genes$log2FoldChange < -1)
significant_genes <- nrow(filtered_genes)

# Afficher les résultats
cat("Gènes upregulés :", upregulated, "\n")
cat("Gènes downregulés :", downregulated, "\n")
cat("Gènes significatifs au total :", significant_genes, "\n")

### Volcano Plot avec mise en évidence des "best_genes"
res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_colour_manual(values = c("gray", "red"),
                      labels = c("Not significant", "Significant")) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "blue") +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             color = "blue")+
  labs(x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value",
       color = "Significance")+
  annotate("text", x = 3, y = 200, label = paste("Upregulated:", upregulated),
           hjust = 0.5, color = "green", size = 4) +
  annotate("text", x = -5, y = 200, label = paste("Downregulated:", downregulated),
           hjust = 0.5, color = "red", size = 4)+
  theme_minimal()+
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  ggtitle("Design1:Effect of Hypoxia In Both Genotypes") +
  theme(legend.position="bottom")+
  geom_text(data=best_genes,aes(label=Geneid),
            vjust=1.5,size=3,color="blue") +
  geom_point(data=best_genes,
             size=1,
             shape=21,
             fill="blue",
             colour="blue")

#Visualisation of design 
sampleData <- data.frame(genotype = rep(c("WT", "KO"), each = 4),
                         culture = rep(c("Hypoxia", "Normoxia"), 4)) 

sampleData$culture <- factor(sampleData$culture, levels = c("Hypoxia", "Normoxia"))

vd <- VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ genotype + culture + genotype:culture,
  textSizeFitted = 3
)
cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

app <- ExploreModelMatrix(
  sampleData = sampleData,
  designFormula = ~ genotype + culture + genotype:culture
)
#> The `name` provided ('') does not correspond to a known icon
#> The `name` provided ('hand-o-right') does not correspond to a known icon
#> The `name` provided ('question-circle fa-1g') does not correspond to a known icon
if (interactive()) {
  shiny::runApp(app)
}

##Comment les genes varient en fonction des conditions d'Hypoxique ou normoxique ?
## Paired design 
dds2 <- DESeqDataSetFromMatrix(countData = mat[ ,metadata$sample],
                               colData = metadata,
                               design =  ~ culture + genotype)


# Set the mock level as the reference level
dds$genotype <- relevel(dds2$genotype, ref = "wt")
# Set the Epithelial cells as the reference level
dds$culture <- relevel(dds2$culture, ref = "normoxia")

dds2 <- DESeq(dds2)

resultsNames(dds2)

res2 <- results(dds2, name = "culture_normoxia_vs_hypoxia")
res_tbl2 <- as_tibble(res2, rownames = "Geneid") %>%
  arrange(padj)
head(res_tbl2)

# normalised counts of the genes with the lowest p-adjusted values are plotted bellow
as_tibble(counts(dds2[res_tbl2$Geneid[1:2], ], normalize = TRUE),
          rownames = 'Geneid') %>%
  pivot_longer(names_to = "sample", values_to = "counts", -Geneid) %>%
  left_join(as_tibble(colData(dds2))) %>%
  mutate(name = paste0(substr(genotype, 1, 5), '_', culture, '_', 1:3)) %>%
  ggplot(aes(x = name, y = log(counts), fill = culture)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ Geneid, scales = "free") +
  theme(axis.text.x = element_text(size = 8, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))

### histogramme
hist(res_tbl2$pvalue)

hist(res_tbl2$padj)

##Best genes
best_genes <- res_tbl2 %>%
  arrange(padj)  %>%
  head(6)


as_tibble(counts(dds[best_genes$Geneid[1], ], normalize = TRUE),
          rownames = 'Geneid') %>%
  pivot_longer(names_to = "sample", values_to = "counts", -Geneid) %>%  
  left_join(as_tibble(colData(dds))) %>% 
  ggplot(aes(x = sample, y = counts, fill = culture)) +
  geom_bar(stat = 'identity', color = "gray30") +
  facet_wrap( ~ genotype, scales = "free", ncol = 2) +
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))

#Visualisation of the 6 best genes
res_tbl2 %>%
  arrange(padj) %>%
  head(6) %>%
  pull(Geneid) %>%
  map(~ {
    as_tibble(counts(dds[.x, ], normalize = TRUE), rownames = 'Geneid') %>%
      pivot_longer(names_to = "sample", values_to = "counts", -Geneid) %>%  
      left_join(as_tibble(colData(dds))) %>% 
      ggplot(aes(x = sample, y = counts, fill = culture)) +
      geom_bar(stat = 'identity', color = "gray30") +
      facet_wrap(~ genotype, scales = "free", ncol = 2) +
      labs(title = paste("Expression du gène :", .x)) +
      theme(axis.text.x = element_text(size = 7, angle = 90),
            axis.title.x = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 7))
  })
##Volcano Plot avec "best_genes" mis en évidence
res_tbl2 %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_colour_manual(values = c("gray", "red"),
                      labels = c("Not significant", "Significant")) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "blue") +
  geom_vline(xintercept = c(-1, 1),
    linetype = "dashed",
    color = "blue")+
  labs(x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value",
    color = "Significance")+
  annotate("text", x = 3, y = 100, label = paste("Upregulated:", upregulated),
           hjust = 0.5, color = "green", size = 4) +
  annotate("text", x = -5, y = 100, label = paste("Downregulated:", downregulated),
           hjust = 0.5, color = "red", size = 4)+
  theme_minimal()+
  theme(legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)) +

  ggtitle("Design1:EffectofHypoxiavsNormoxiaConditions") +
  theme(legend.position="bottom")+
  geom_text(data=best_genes,aes(label=Geneid),
            vjust=1.5,size=3,color="blue") +
  geom_point(data=best_genes,
             size=1,
             shape=21,
             fill="blue",
             colour="blue")


#Adjust pvalue 
res_tblcorrected <- res_tbl %>%
  filter(pvalue > 0.85 & pvalue < 0.9) %>% 
  head(10)

hist(res_tblcorrected$pvalue)

#Volcano plot adjusted 
res_tblcorrected %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_colour_manual(values = c("gray", "red")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  theme(legend.position = "bottom")

# Visualisation design 
sampleData <- data.frame(genotype = rep(c("WT", "KO"), each = 4),
                          culture = rep(c("Normoxia", "Hypoxia"), 4))

vd <- VisualizeDesign(sampleData = metadata, 
                      designFormula = ~ genotype + culture, 
                      textSizeFitted = 4)
cowplot::plot_grid(plotlist = vd$plotlist)

app <- ExploreModelMatrix(sampleData = metadata,
                          designFormula = ~ genotype + culture)
#> The `name` provided ('') does not correspond to a known icon
#> The `name` provided ('hand-o-right') does not correspond to a known icon
#> The `name` provided ('question-circle fa-1g') does not correspond to a known icon
if (interactive()) {
  shiny::runApp(app)
}

## Analyse d'enrichissement de res_tbl de la Question 1 : Ouverture de la librairie de g?nes de souris
library(org.Mm.eg.db)
library("GO.db")

# Cr?ation du nouveau tableau avec les g?nes et le ENTREZID incorpor?
res_tbl_e <- res_tbl %>% 
  mutate(gene = mapIds(org.Mm.eg.db, Geneid, "SYMBOL", "ENSEMBL")) %>% 
  mutate(ENTREZID = mapIds(org.Mm.eg.db, Geneid, "ENTREZID", "ENSEMBL")) %>% 
  dplyr::select(Geneid, gene, ENTREZID, everything())

res_tbl_e

# Ici on retire les duplicats avec le m?me Geneid et on retire les p-value ajust?es = NA
res_tbl_E <- res_tbl_e %>%
  filter(!is.na(Geneid),
         !is.na(padj),
         !duplicated(Geneid)) %>%
  mutate(Geneid = as.character(Geneid))

res_tbl_E

# GO:0071456 = R?ponse cellulaire ? l'hypoxie 
# GO:0045765 = R?gulation de l'angiogen?se 
# GO:0001666 = R?ponse ? l'hypoxie

#we have 80 mice genes matching this Go term
GO_0071456 <- AnnotationDbi::select(org.Mm.eg.db,
                                    keys = "GO:0071456",
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "GO") %>%
  as_tibble %>%
  filter(!duplicated(ENTREZID))
GO_0071456 

#We have 36 mice genes matching this GO term
GO_0045765 <- AnnotationDbi::select(org.Mm.eg.db,
                                    keys = "GO:0045765",
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "GO") %>%
  as_tibble %>%
  filter(!duplicated(ENTREZID))
GO_0045765

#We have 123 mice genes matching this GO_term
GO_0001666 <- AnnotationDbi::select(org.Mm.eg.db,
                                    keys = "GO:0001666",
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "GO") %>%
  as_tibble %>%
  filter(!duplicated(ENTREZID))
GO_0001666

###Second analysis of enrichment using clusterprofil

library(clusterProfiler)

#GESA with GO set
ordered_genes <- abs(res_tbl_E$stat)
names(ordered_genes) <- res_tbl_E$ENTREZID
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

go_gsea <- gseGO(gene = ordered_genes,
                 OrgDb = org.Mm.eg.db,
                 scoreType = "pos") %>%
  as_tibble
go_gsea

#ORA using go set 

de_genes <- res_tbl_E$ENTREZID[res_tbl_E$padj < 0.05]
go_ora <- enrichGO(gene = de_genes,
                   universe = res_tbl_E$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   readable = TRUE) %>%
  as_tibble
go_ora

# GSEA analysis above using all GO namespaces
go_cc <- full_join(go_ora %>%
                     filter(p.adjust < 0.05) %>%
                     dplyr::select(ID, p.adjust),
                   go_gsea %>%
                     filter(p.adjust < 0.05) %>%
                     dplyr::select(ID, p.adjust),
                   by = c("ID" = "ID"))
go_cc
na.omit(go_cc)

#Comparing ORA to GESA 
filter(go_cc, is.na(p.adjust.y)) # Terme manquant dans GSEA
filter(go_cc, is.na(p.adjust.x)) # Terme manquant dans ORA

##Visualisation of enrichment analysis-highlight the genes in a particular set of interest on the volcano plot
#For GO_0071456
sel <- res_tbl_E$ENTREZID %in% GO_0071456$ENTREZID
plot(res_tbl_E$log2FoldChange, -log10(res_tbl_E$padj))
points(res_tbl_E$log2FoldChange[sel],
       -log10(res_tbl_E$padj)[sel],
       col = "red")
grid()

# For GO_0045765 
sel <- res_tbl_E$ENTREZID %in% GO_0045765 $ENTREZID
plot(res_tbl_E$log2FoldChange, -log10(res_tbl_E$padj))
points(res_tbl_E$log2FoldChange[sel],
       -log10(res_tbl_E$padj)[sel],
       col = "red")
grid()

#For GO_0001666
sel <- res_tbl_E$ENTREZID %in% GO_0001666 $ENTREZID
plot(res_tbl_E$log2FoldChange, -log10(res_tbl_E$padj))
points(res_tbl_E$log2FoldChange[sel],
       -log10(res_tbl_E$padj)[sel],
       col = "red")
grid()

#Visualization of the enrichment analysis/Diagramme de venn 
# Extraire les gènes significatifs pour ORA et GSEA
ora_genes <- go_ora$ID[go_ora$p.adjust < 0.05]
gsea_genes <- go_gsea$ID[go_gsea$p.adjust < 0.05]

# Créer les données pour le diagramme de Venn
venn_data <- list(
  ORA = ora_genes,  # Gènes significatifs de ORA
  GSEA = gsea_genes  # Gènes significatifs de GSEA
)

# Générer le diagramme de Venn avec des couleurs personnalisées
ggVennDiagram(venn_data, 
              label_alpha = 0.4,    # Transparence des étiquettes
              label_col = "black",  # Couleur des étiquettes
              category.names = c("ORA", "GSEA"),  # Noms des catégories
              fill_color = c("skyblue", "orange"))  # Couleurs personnalisées pour ORA et GSEA


###QUESTION2: Enrichment analysis

res_tbl_e2 <- res_tbl2 %>% 
  mutate(gene = mapIds(org.Mm.eg.db, Geneid, "SYMBOL", "ENSEMBL")) %>% 
  mutate(ENTREZID = mapIds(org.Mm.eg.db, Geneid, "ENTREZID", "ENSEMBL")) %>% 
  dplyr::select(Geneid, gene, ENTREZID, everything())

res_tbl_e2

res_tbl_E2 <- res_tbl_e2 %>%
  filter(!is.na(Geneid),
         !is.na(padj),
         !duplicated(Geneid)) %>%
  mutate(Geneid = as.character(Geneid))

res_tbl_E2

#GESA with GO set
ordered_genes2 <- abs(res_tbl_E2$stat)
names(ordered_genes2) <- res_tbl_E2$ENTREZID
ordered_genes2 <- sort(ordered_genes, decreasing = TRUE)

go_gsea2 <- gseGO(gene = ordered_genes2,
                 OrgDb = org.Mm.eg.db,
                 scoreType = "pos") %>%
  as_tibble
go_gsea2

#ORA using go set 

de_genes2 <- res_tbl_E2$ENTREZID[res_tbl_E2$padj < 0.05]
go_ora2 <- enrichGO(gene = de_genes2,
                   universe = res_tbl_E2$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   ont = "ALL",
                   readable = TRUE) %>%
  as_tibble
go_ora2

# GSEA analysis above using all GO namespaces
go_cc2 <- full_join(go_ora2 %>%
                     filter(p.adjust < 0.01) %>%
                     dplyr::select(ID, p.adjust),
                   go_gsea2 %>%
                     filter(p.adjust < 0.01) %>%
                     dplyr::select(ID, p.adjust),
                   by = c("ID" = "ID"))
go_cc2
na.omit(go_cc2)

#Comparing ORA to GESA 
filter(go_cc2, is.na(p.adjust.y)) 
filter(go_cc2, is.na(p.adjust.x))


##Visualisation of enrichment analysis-highlight the genes in a particular set of interest on the volcano plot
#For GO_0071456
sel2 <- res_tbl_E2$ENTREZID %in% GO_0071456$ENTREZID
plot(res_tbl_E2$log2FoldChange, -log10(res_tbl_E2$padj))
points(res_tbl_E2$log2FoldChange[sel2],
       -log10(res_tbl_E2$padj)[sel2],
       col = "red")
grid()

# For GO_0045765 
sel2 <- res_tbl_E2$ENTREZID %in% GO_0045765 $ENTREZID
plot(res_tbl_E2$log2FoldChange, -log10(res_tbl_E2$padj))
points(res_tbl_E2$log2FoldChange[sel2],
       -log10(res_tbl_E2$padj)[sel2],
       col = "red")
grid()

#For GO_0001666
sel2 <- res_tbl_E2$ENTREZID %in% GO_0001666 $ENTREZID
plot(res_tbl_E2$log2FoldChange, -log10(res_tbl_E2$padj))
points(res_tbl_E$log2FoldChange[sel2],
       -log10(res_tbl_E2$padj)[sel2],
       col = "red")
grid()



#Visualization of the enrichment analysis/Diagramme de venn 
# Extraire les gènes significatifs pour ORA et GSEA
ora_genes2 <- go_ora2$ID[go_ora2$p.adjust < 0.05]
gsea_genes2 <- go_gsea2$ID[go_gsea2$p.adjust < 0.05]

# Créer les données pour le diagramme de Venn
venn_data <- list(
  ORA = ora_genes2,  # Gènes significatifs de ORA
  GSEA = gsea_genes2  # Gènes significatifs de GSEA
)

# Générer le diagramme de Venn avec des couleurs personnalisées
ggVennDiagram(venn_data, 
              label_alpha = 0.4,    # Transparence des étiquettes
              label_col = "black",  # Couleur des étiquettes
              category.names = c("ORA", "GSEA"),  # Noms des catégories
              fill_color = c("skyblue", "orange"))  # Couleurs personnalisées pour ORA et GSEA


