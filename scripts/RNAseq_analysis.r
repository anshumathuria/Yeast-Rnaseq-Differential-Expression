###############################################################
# RNA-seq Differential Expression Analysis
#
# Organism : Saccharomyces cerevisiae
# Author   : Anshu Mathuria
#
# Workflow
# 1. Load Required Libraries
# 2. Read Counting using featureCounts
# 3. RPM & RPKM Normalization
# 4. Differential Expression Analysis (edgeR)
# 5. Functional Enrichment Preparation
# 6. Heatmap Visualization
# 7. Volcano Plot Visualization
###############################################################

# Step 1: Load Required Libraries
library(Rsubread); library(edgeR); library(limma); library(pheatmap); library(ggplot2)

# Step 2: Read Counting using featureCounts
gtf_file <- "GCF_000146045.2_R64_genomic.gtf"

conA_FeatCount_rep1 <- featureCounts(files="conA_rep1.sam", annot.ext=gtf_file, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id")
conA_FeatCount_rep2 <- featureCounts(files="conA_rep2.sam", annot.ext=gtf_file, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id")
conB_FeatCount_rep1 <- featureCounts(files="conB_rep1.sam", annot.ext=gtf_file, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id")
conB_FeatCount_rep2 <- featureCounts(files="conB_rep2.sam", annot.ext=gtf_file, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id")

# Combine Read Counts
Total_FeatCount <- data.frame(conA_FeatCount_rep1$counts, conA_FeatCount_rep2$counts, conB_FeatCount_rep1$counts, conB_FeatCount_rep2$counts)
write.table(Total_FeatCount, file="Total_FeatCount_conA_conB_rep1_rep2.txt", sep="\t", quote=FALSE)

# Step 3: RPM Normalization
sf_conA_rep1 <- sum(Total_FeatCount$conA_rep1.sam) / 1000000
sf_conA_rep2 <- sum(Total_FeatCount$conA_rep2.sam) / 1000000
sf_conB_rep1 <- sum(Total_FeatCount$conB_rep1.sam) / 1000000
sf_conB_rep2 <- sum(Total_FeatCount$conB_rep2.sam) / 1000000

Total_FeatCount$conA_rep1_RPM <- Total_FeatCount$conA_rep1.sam / sf_conA_rep1
Total_FeatCount$conA_rep2_RPM <- Total_FeatCount$conA_rep2.sam / sf_conA_rep2
Total_FeatCount$conB_rep1_RPM <- Total_FeatCount$conB_rep1.sam / sf_conB_rep1
Total_FeatCount$conB_rep2_RPM <- Total_FeatCount$conB_rep2.sam / sf_conB_rep2

# -------------------------------------------------------------
# RPKM Normalization & Gene Length Extraction
# -------------------------------------------------------------
# EXTRACT GENE LENGTHS directly from the featureCounts output instead of loading a CSV
gene_data <- data.frame(Gene = conA_FeatCount_rep1$annotation$GeneID, Length_bp = conA_FeatCount_rep1$annotation$Length)

gene_data$gene_length_kb <- gene_data$Length_bp / 1000
Total_FeatCount$Gene <- rownames(Total_FeatCount)
merged_data <- merge(Total_FeatCount, gene_data, by="Gene")

merged_data$conA_rep1_RPKM <- merged_data$conA_rep1_RPM / merged_data$gene_length_kb
merged_data$conA_rep2_RPKM <- merged_data$conA_rep2_RPM / merged_data$gene_length_kb
merged_data$conB_rep1_RPKM <- merged_data$conB_rep1_RPM / merged_data$gene_length_kb
merged_data$conB_rep2_RPKM <- merged_data$conB_rep2_RPM / merged_data$gene_length_kb

write.csv(merged_data, "Normalized_RPKM.csv", row.names=FALSE)

# Step 4: Differential Expression Analysis using edgeR
x <- read.delim("Total_FeatCount_conA_conB_rep1_rep2.txt", header=TRUE)
group <- c(rep("ConA",2), rep("ConB",2))
y2 <- DGEList(counts=x, group=group)
keep <- filterByExpr(y2)
y2 <- y2[keep, , keep.lib.sizes=FALSE]
y2 <- calcNormFactors(y2)
design2 <- model.matrix(~group)
y2 <- estimateDisp(y2, design2)
DEet <- exactTest(y2, pair=c("ConA","ConB"))
DEet <- topTags(DEet, n=Inf, adjust.method="BH")
et_table <- DEet$table

# Extract Significant DEGs
up_et <- subset(et_table, logFC >= 1 & FDR <= 0.05)
down_et <- subset(et_table, logFC <= -1 & FDR <= 0.05)

write.csv(et_table, "All_DEGs.csv")
write.csv(up_et, "Upregulated_Genes.csv")
write.csv(down_et, "Downregulated_Genes.csv")

# Step 5: Functional Enrichment Preparation
up_gene_list <- rownames(up_et)
down_gene_list <- rownames(down_et)
write.table(up_gene_list, "up_regulated_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(down_gene_list, "down_regulated_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Step 6: Heatmap Visualization
all_DEG_list <- c(up_gene_list, down_gene_list)
logCPM_All <- cpm(y2, log=TRUE, normalized.lib.sizes=TRUE)[all_DEG_list, ]
annotation_col <- data.frame(Condition = factor(group))
rownames(annotation_col) <- colnames(logCPM_All)

png("Heatmap.png", width=1000, height=900)
pheatmap(logCPM_All, scale="row", annotation_col=annotation_col, show_rownames=FALSE, cluster_rows=TRUE, cluster_cols=TRUE, color=colorRampPalette(c("blue","white","red"))(100), main="Heatmap of Differentially Expressed Genes")
dev.off()

# Step 7: Volcano Plot
et_table$diffexpressed <- "Not Significant"
et_table$diffexpressed[et_table$logFC > 1 & et_table$FDR < 0.05] <- "Up"
et_table$diffexpressed[et_table$logFC < -1 & et_table$FDR < 0.05] <- "Down"

volcano_plot <- ggplot(et_table, aes(x=logFC, y=-log10(FDR), color=diffexpressed)) + geom_point(size=2) + theme_minimal() + scale_color_manual(values=c("blue", "grey", "red")) + labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10(FDR)", color="Expression")
ggsave("Volcano_Plot.png", volcano_plot, width=8, height=6, dpi=300)

###############################################################
# Analysis Completed Successfully
###############################################################
