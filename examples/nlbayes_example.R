library(rjson)
library(org.Hs.eg.db)
library(nlbayes)

# download files: 
#   - This differential expression table was generated using the GEO2R tool, contrasting the E2F3 treated
#     samples against the control samples. At GEO2R, we need to select the `Gene.ID` column that contains
#     Entrez (NCBI) gene ids.
#     url: https://umbibio.math.umb.edu/nlbayes/assets/data/experiments/GSE3151.E2F3.top.table.tsv
#   - Since the experiment above was performed on mammary epithelial cell cultures, we may choose a breast
#     tissue specific network. 
#     url: https://umbibio.math.umb.edu/nlbayes/assets/data/networks/gtex_chip/homo_sapiens/tissue_specific/breast.rels.json

# load the network. It must contain modes of regulation as integers: -1, 0, 1
network <- fromJSON(file='breast.rels.json')

# load a table with differential expression
deg.table <- read.table('GSE3151.E2F3.top.table.tsv', sep = '\t', header = TRUE, quote = NULL)

# filter DE genes by setting cutoff values on p-value and/or fold-change
top.deg.table <- deg.table[(abs(deg.table$logFC) >= 1.) & (deg.table$adj.P.Val <= 0.001), ]

# remove duplicated entries
top.deg.table <- top.deg.table[!duplicated(top.deg.table$Gene.ID), ]

# construct the evidence vector
evidence <- sign(top.deg.table$logFC)
names(evidence) <- top.deg.table$Gene.ID

# build the inference model object container
inference.model <- ORNOR.inference(network, evidence, n.graphs = 5, verbosity=2)

# sample the posterior distributions
inference.model <- sample.posterior(inference.model, N = 2000, gr.level = 1.15, burnin = TRUE)

# we would like to have the results annotated with TF symbols for easy analysis
annotation.src <- select(org.Hs.eg.db, keys(org.Hs.eg.db, 'ENTREZID'), 'SYMBOL')
annotation <- as.list(annotation.src$SYMBOL)
names(annotation) <- annotation.src$ENTREZID

# we run this method to generate the output inference results table
inference.model <- postprocess.result(inference.model, annotation)

# finally, we write the results to a tsv file
tf.inference <- inference.model$result.info$tf.inference
write.table(tf.inference, 'inference_result.tsv', sep = '\t', row.names = FALSE)


# optionally, we may sample the posterior distribution right away with the ORNOR.inference method
# inference.model <- ORNOR.inference(network, evidence, n.graphs = 5, verbosity=2, do.sample = TRUE, N = 2000, gr.level = 1.15, burnin = TRUE)
# inference.model <- postprocess.result(inference.model, annotation)

