# Load required libraries
library(edgeR)
library(openxlsx)

# test the DEG of oil with ctrl. 4 tests: with uv --- high vs. ctrl; low vs. ctrl; without uv --- high vs. ctrl; low vs. ctrl

# define the function of DEA
DEA <- function(control, treat, count.control, count.treat) {
  num.control <- ncol(count.control)
  num.treat <- ncol(count.treat)
  
  count.table <- cbind(count.control, count.treat)
  
  gene.list <- rownames(count.table)
  group <- relevel(factor(c(rep(control, num.control), rep(treat, num.treat))), ref = control)
  design <- model.matrix(~group)
  
  # Put the data into a DGEList object
  y <- DGEList(counts = count.table, genes = gene.list)
  
  # Filtering
  countsPerMillion <- cpm(y)
  countCheck <- countsPerMillion > 1
  keep <- which(rowSums(countCheck) > 1)
  y <- y[keep, ]
  
  # Normalization
  y <- calcNormFactors(y, method="TMM")
  
  y$samples$group <- group
  
  rownames(design) <- colnames(y)
  
  # Estimating the dispersion
  
  # estimate the NB dispersion for the dataset
  y <- estimateDisp(y, design, robust = TRUE)
  
  # Differential expression
  
  # determine differentially expressed genes
  # fit genewise glms
  fit <- glmFit(y, design)
  
  # conduct likelihood ratio tests for tumour vs normal tissue differences and show the top genes
  lrt <- glmLRT(fit)
  
  # the DEA result for all the genes
  # dea <- lrt$table
  toptag <- topTags(lrt, n = nrow(y$genes), p.value = 1)
  dea <- toptag$table  # just to add one more column of FDR
  dea <- dea[order(dea$FDR, -abs(dea$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC
  
  return(dea)
}

# Specify the input and output paths
path.input = 'data/oil_uv'
path.output = 'output/oil_uv'

# Load data and merge all tables together

## metafile
meta.data <- read.xlsx(file.path(path.input, 'metadata.xlsx'))
samples <- meta.data$sample
group.oil <- meta.data$oil
group.uv <- meta.data$uv

## read all the files
merge.file<-NULL
file.names <- dir(path.input, pattern = '.tsv')
for(i in 1:length(file.names)){
  file <- read.table(file.path(path.input, file.names[i]),header=TRUE, row.names = 1, sep="\t", stringsAsFactors=FALSE)
  if (i == 1) {
    merge.file <- file
  } else {
    merge.file <- cbind(merge.file, file) 
  }
}

# reorder the merged table based on the metadata
merge.table <- merge.file[, samples]
gene.list <- rownames(merge.table)

# ================== effect of low oil without uv==================
uv <- "no"
control <- "ctrl"
treat <- "low"
count.control <- merge.table[, which(colnames(merge.table) == meta.data$sample[meta.data$uv == uv & meta.data$oil == control])]
count.treat <- merge.table[, which(colnames(merge.table) == meta.data$sample[meta.data$uv == uv & meta.data$oil == treat])]

dea <- DEA(control, treat, count.control, count.treat)

write.table(dea, paste(path.output, '/dea_low_without_uv.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')

# ================== effect of high oil without uv==================
uv <- "no"
control <- "ctrl"
treat <- "high"
count.control <- merge.table[, which(colnames(merge.table) == meta.data$sample[meta.data$uv == uv & meta.data$oil == control])]
count.treat <- merge.table[, which(colnames(merge.table) == meta.data$sample[meta.data$uv == uv & meta.data$oil == treat])]

dea <- DEA(control, treat, count.control, count.treat)

write.table(dea, paste(path.output, '/dea_high_without_uv.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')

# ================== effect of low oil with uv==================
uv <- "yes"
control <- "ctrl"
treat <- "low"
count.control <- merge.table[, which(colnames(merge.table) == meta.data$sample[meta.data$uv == uv & meta.data$oil == control])]
count.treat <- merge.table[, which(colnames(merge.table) == meta.data$sample[meta.data$uv == uv & meta.data$oil == treat])]

dea <- DEA(control, treat, count.control, count.treat)

write.table(dea, paste(path.output, '/dea_low_with_uv.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')

# ================== effect of high oil with uv==================
uv <- "yes"
control <- "ctrl"
treat <- "high"
count.control <- merge.table[, which(colnames(merge.table) == meta.data$sample[meta.data$uv == uv & meta.data$oil == control])]
count.treat <- merge.table[, which(colnames(merge.table) == meta.data$sample[meta.data$uv == uv & meta.data$oil == treat])]

dea <- DEA(control, treat, count.control, count.treat)

write.table(dea, paste(path.output, '/dea_high_with_uv.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')

