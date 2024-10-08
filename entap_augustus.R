# matching AUGUSTUS and ENTAP output into a unique gtf

# libraries
library(data.table)
library(dplyr)

# reading the files
tsv <- fread(file = "entap_results.tsv", header = FALSE)
tsv <- tsv[-1,]
tsv <- tsv[,c(1,13)]
gtf <- fread(file = "cbez_genome_annot.gtf")

# updated gtf
new_gtf <- left_join(gtf, tsv, by = c("V9" = "V1"))

# write gtf
fwrite(x = new_gtf, quote = FALSE, sep = '\t', row.names = FALSE,
       col.names = FALSE, file = "cbez_final_annot.gtf")