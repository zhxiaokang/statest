# Load required libraries
library(edgeR)
library(openxlsx)

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

# Put the data into a DGEList object
y <- DGEList(counts = merge.table, genes = gene.list)
y$samples$oil <- group.oil
y$samples$uv <- group.uv
y$samples$group <- paste(group.oil, group.uv, sep = '_')

# Filtering
countsPerMillion <- cpm(y)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) > 1)
y <- y[keep, ]

# Normalization
y <- calcNormFactors(y, method="TMM")

# check the meta info:
print(y$samples)

# ================== effect of uv on oil ==================
# design matrix: define how to run the test
design <- model.matrix(~uv + oil + uv:oil, data = meta.data)

# fit the model
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design)

# run the test
qlf <- glmQLFTest(fit, coef=5:6)

# the DEA result for all the genes
toptag <- topTags(qlf, n = nrow(y$genes), p.value = 1)
dea <- toptag$table  # just to add one more column of FDR

# save the DEA result to files
write.table(dea, paste(path.output, '/dea_uv_effect_on_oil.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')



# ================== effect of oil without uv==================
# subtract the data
merge.table.uvno <- merge.table[, 1:9]
meta.data.unvo <- meta.data[1:9, ]

# Put the data into a DGEList object
y <- DGEList(counts = merge.table.uvno, genes = gene.list)
y$samples$oil <- group.oil[1:9]

# Filtering
countsPerMillion <- cpm(y)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) > 1)
y <- y[keep, ]

# Normalization
y <- calcNormFactors(y, method="TMM")

# check the meta info:
print(y$samples)

# design matrix: define how to run the test
design <- model.matrix(~oil, data = meta.data.unvo)

# fit the model
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design)

# run the test
qlf <- glmQLFTest(fit, coef=2:3)

# the DEA result for all the genes
toptag <- topTags(qlf, n = nrow(y$genes), p.value = 1)
dea <- toptag$table  # just to add one more column of FDR

# save the DEA result to files
write.table(dea, paste(path.output, '/dea_oil_effect_without_uv.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')

# ================== effect of oil with uv==================
# subtract the data
merge.table.uvno <- merge.table[, 10:18]
meta.data.unvo <- meta.data[10:18, ]

# Put the data into a DGEList object
y <- DGEList(counts = merge.table.uvno, genes = gene.list)
y$samples$oil <- group.oil[10:18]

# Filtering
countsPerMillion <- cpm(y)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) > 1)
y <- y[keep, ]

# Normalization
y <- calcNormFactors(y, method="TMM")

# check the meta info:
print(y$samples)

# design matrix: define how to run the test
design <- model.matrix(~oil, data = meta.data.unvo)

# fit the model
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design)

# run the test
qlf <- glmQLFTest(fit, coef=2:3)

# the DEA result for all the genes
toptag <- topTags(qlf, n = nrow(y$genes), p.value = 1)
dea <- toptag$table  # just to add one more column of FDR

# save the DEA result to files
write.table(dea, paste(path.output, '/dea_oil_effect_with_uv.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')

# ================== interaction effect of uv and high oil ==================
# re-define the variables to avoid unexpected modifiction to them
merge.table <- merge.file[, samples]
gene.list <- rownames(merge.table)
meta.data <- read.xlsx(file.path(path.input, 'metadata.xlsx'))

# subtract the data
merge.table <- merge.table[, subset(meta.data, oil == "ctrl" | oil == "high")$sample]
meta.data <- subset(meta.data, oil == "ctrl" | oil == "high")

group.oil <- meta.data$oil
group.uv <- meta.data$uv

group.oil <- relevel(factor(group.oil), ref = "ctrl")
group.uv <- relevel(factor(group.uv), ref = "no")

# Put the data into a DGEList object
y <- DGEList(counts = merge.table, genes = gene.list)
y$samples$oil <- group.oil
y$samples$uv <- group.uv
y$samples$group <- paste(group.oil, group.uv, sep = '_')

# Filtering
countsPerMillion <- cpm(y)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) > 1)
y <- y[keep, ]

# Normalization
y <- calcNormFactors(y, method="TMM")

# check the meta info:
print(y$samples)

# design matrix: define how to run the test
design <- model.matrix(~oil * uv, data = meta.data)

# fit the model
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design)

# run the test
qlf <- glmQLFTest(fit, coef=4)

# the DEA result for all the genes
toptag <- topTags(qlf, n = nrow(y$genes), p.value = 1)
dea <- toptag$table  # just to add one more column of FDR

# save the DEA result to files
write.table(dea, paste(path.output, '/dea_interact_high_uv.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')




