# UV yes V.S. UV no (oil as ctrl)

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

# ================== effect of uv without oil==================
# groups to compare: ctrl_no V.S. ctrl_yes
# so relevant samples are: sample 1,2,3 V.S. 4,5,6
# subtract the data
merge.table.ctrl <- merge.table[, c(1:3, 10:12)]
meta.data.ctrl <- meta.data[c(1:3, 10:12), ]


# ================================== glmQLFTest ====================

# Put the data into a DGEList object
y <- DGEList(counts = merge.table.ctrl, genes = gene.list)
y$samples$uv <- group.uv[c(1:3, 10:12)]

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
design <- model.matrix(~uv, data = meta.data.ctrl)

# fit the model
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design)

# run the test
qlf <- glmQLFTest(fit)

# the DEA result for all the genes
toptag <- topTags(qlf, n = nrow(y$genes), p.value = 1)
dea <- toptag$table  # just to add one more column of FDR

# save the DEA result to files
write.table(dea, paste(path.output, '/dea_uv_effect_ctrl_oil.tsv', sep = ''), row.names = F, quote = FALSE, sep = '\t')




