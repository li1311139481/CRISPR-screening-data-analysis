#!/usr/bin/env Rscript

# Summarize MAGeCK outputs and generate QC/visualization reports.
# Usage:
# Rscript scripts/mageck_flute.R <workdir> <sample_label_csv> <comparison_name>

.libPaths(R.home("library"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript mageck_flute.R <workdir> <sample_label_csv> <comparison_name>")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(data.table)
  library(MAGeCKFlute)
})

workdir <- as.character(args[1])
sample_level <- unlist(strsplit(args[2], ","))
header_name <- args[3]

setwd(workdir)
workdir <- getwd()
message("[INFO] Working directory: ", workdir)
message("[INFO] Sample order: ", paste(sample_level, collapse = ","))
message("[INFO] Comparison name: ", header_name)

# Convert gene symbols and write annotated full gene summary.
rra.gene_summary <- read.table(sprintf("%s/%s_test.gene_summary.txt", workdir, header_name), header = TRUE)

id <- rra.gene_summary$id %>%
  bitr(fromType = "SYMBOL", toType = c("ENTREZID", "GENENAME"), OrgDb = "org.Mm.eg.db", drop = FALSE) %>%
  distinct(SYMBOL, .keep_all = TRUE)

gene_summary_all <- merge(data.table(id), data.table(rra.gene_summary), by.x = "SYMBOL", by.y = "id") %>% .[order(pos.rank)]
colnames(gene_summary_all) <- c(
  colnames(gene_summary_all)[1:3],
  paste(colnames(gene_summary_all)[4:length(colnames(gene_summary_all))], header_name, sep = "_")
)
write.csv(gene_summary_all, file = sprintf("%s/%s_gene_summary.csv", workdir, header_name), row.names = FALSE)

# Create a unified positive/negative gene summary for plotting and hit prioritization.
rra.pos.gene <- rra.gene_summary[which(rra.gene_summary$pos.lfc >= 0), c("id", "num", "pos.score", "pos.p.value", "pos.fdr", "neg.rank", "pos.rank", "pos.lfc")]
colnames(rra.pos.gene) <- c("id", "num", "score", "p.value", "fdr", "negrank", "posrank", "lfc")
rra.neg.gene <- rra.gene_summary[which(rra.gene_summary$pos.lfc < 0), c("id", "num", "neg.score", "neg.p.value", "neg.fdr", "neg.rank", "pos.rank", "neg.lfc")]
colnames(rra.neg.gene) <- c("id", "num", "score", "p.value", "fdr", "negrank", "posrank", "lfc")
gene_bind <- rbind(rra.pos.gene, rra.neg.gene)

id2 <- gene_bind$id %>%
  bitr(fromType = "SYMBOL", toType = c("ENTREZID", "GENENAME"), OrgDb = "org.Mm.eg.db", drop = FALSE) %>%
  distinct(SYMBOL, .keep_all = TRUE)

gene_summary <- cbind(id2, gene_bind[, c("num", "lfc", "score", "fdr", "negrank", "posrank")])
colnames(gene_summary) <- c("SYMBOL", "ENTREZID", "GENENAME", "NumOfValidGuide", "Log2FC", "RRAscore", "fdr", "neg.rank", "pos.rank")
write.csv(gene_summary, file = sprintf("%s/%s_mix_gene_summary.csv", workdir, header_name), row.names = FALSE)

# QC plots from count summary.
countsummary <- read.delim(sprintf("%s/%s_count.countsummary.txt", workdir, header_name), check.names = FALSE)
countsummary$Label <- factor(countsummary$Label, ordered = TRUE, levels = sample_level)

p1 <- MapRatesView(countsummary) + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave("MapRatesView.pdf", plot = p1, width = 5, height = 5)

p2 <- IdentBarView(countsummary, x = "Label", y = "GiniIndex", ylab = "Gini index", main = "Evenness of sgRNA reads")
ggsave("GiniIndexView.pdf", plot = p2, width = 5, height = 5)

countsummary$Missed <- round(countsummary$Zerocounts / countsummary$TotalsgRNAs, 2)
p3 <- ggplot(data = countsummary, mapping = aes(x = Label, y = Missed)) +
  geom_bar(stat = "identity", fill = "#004DA7", color = "black") +
  geom_text(aes(label = scales::percent(Missed)), vjust = -0.5, size = 5, color = "black") +
  labs(title = header_name, y = "Missed ratio") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15)
  )
ggsave("MissedsgRNAView.pdf", plot = p3, width = 5, height = 5)

# Volcano plots and summary PDF.
rra.sgrna_summary <- read.table(sprintf("%s/%s_test.sgrna_summary.txt", workdir, header_name), header = TRUE)

p4 <- VolcanoView(rra.gene_summary, x = "pos.lfc", y = "pos.score", Label = "id", ylab = "-log10(RRA.score)", main = "Pos_filter")
ggsave("VolcanoView.pos.score.pdf", plot = p4, width = 5, height = 5)

p5 <- VolcanoView(rra.gene_summary, x = "neg.lfc", y = "neg.score", Label = "id", ylab = "-log10(RRA.score)", main = "Neg_filter")
ggsave("VolcanoView.neg.score.pdf", plot = p5, width = 5, height = 5)

rra.pos.gene2 <- rra.gene_summary[which(rra.gene_summary$pos.lfc > 0), c("id", "pos.score", "pos.p.value", "pos.fdr", "pos.rank", "pos.lfc")]
colnames(rra.pos.gene2) <- c("id", "score", "p.value", "fdr", "rank", "lfc")
rra.neg.gene2 <- rra.gene_summary[which(rra.gene_summary$pos.lfc < 0), c("id", "neg.score", "neg.p.value", "neg.fdr", "neg.rank", "neg.lfc")]
colnames(rra.neg.gene2) <- c("id", "score", "p.value", "fdr", "rank", "lfc")
gene_bind2 <- rbind(rra.pos.gene2, rra.neg.gene2)
pos.top <- as.character(rra.pos.gene2[which(rra.pos.gene2$rank <= 10), 1])
neg.top <- as.character(rra.neg.gene2[which(rra.neg.gene2$rank <= 10), 1])
p6 <- VolcanoView(gene_bind2, x = "lfc", y = "score", Label = "id", topnames = c(pos.top, neg.top), ylab = "-Log10(RRA.score)") +
  ggtitle(header_name) +
  theme(plot.title = element_text(size = 8))
ggsave("VolcanoView.pdf", plot = p6, width = 5, height = 5)

pdf(sprintf("%s_MAGeCKFlute_report.pdf", header_name), width = 8, height = 6)
plot(p1)
plot(p2)
plot(p3)
plot(p4)
plot(p5)
plot(p6)
dev.off()

message("[INFO] all done")
