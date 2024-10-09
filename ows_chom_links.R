# BUSCO table parsing for macrosynteny

# libraries
library(data.table)

# reading raw BUSCO tables
ows <- fread(file = "links_full_table.tsv")
chom <- fread(file = "chom_full_table.tsv")

df_list <- list(ows, chom)

# filtering the tsv's 
df_list[[2]] <- df_list[[2]][- grep("Duplicated|Fragmented", df_list[[2]][[2]]),]

for (i in 1:length(df_list)) {
  df_list[[i]] <- df_list[[i]][,-c(2, 6, 7, 8, 9, 10)]
}

# getting the common BUSCO genes between species
com <- df_list[[1]][df_list[[1]][[1]] %in% df_list[[2]][[1]], 1]

# comparasions with Chom
comp_ows_chom <- cbind(df_list[[2]][V1 %in% com$V1]
                       [order(V1), c("V3", "V4", "V5")],
                       df_list[[1]][V1 %in% com$V1]
                       [order(V1), c("V3", "V4", "V5")])

# writing tables
fwrite(x = comp_ows_chom, quote = FALSE, sep = '\t', row.names = FALSE,
       col.names = FALSE, file = "ows_chom.txt")
