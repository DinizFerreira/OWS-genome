# BUSCO table parsing for macrosynteny

# libraries
library(data.table)
library(dplyr)

# reading raw BUSCO tables
ows <- fread(file = "full_table.tsv")
chom <- fread(file = "chom_full_table.tsv")

df_list <- list(ows, chom)

# filtering the tsv's 
df_list[[2]] <- df_list[[2]][- grep("Duplicated|Fragmented", df_list[[2]][[2]]),]

for (i in 1:length(df_list)) {
  df_list[[i]] <- df_list[[i]][,-c(2, 6, 7, 8)]
}

# getting the common BUSCO genes between species
com <- df_list[[1]][df_list[[1]][[1]] %in% df_list[[2]][[1]], 1]

# comparasions with Chom
comp_ows_chom <- merge(df_list[[2]][V1 %in% com$V1, .(V1, chom_V3 = V3)], 
                       df_list[[1]][V1 %in% com$V1, .(V1, ows_V3 = V3)], 
                       by = "V1", sort = TRUE)

colnames(comp_ows_chom) <- c('BUSCO id','Chromosome name','OWS scaffold')

comp_ows_chom <- comp_ows_chom %>%
  mutate(`Chromosome name`=ifelse(`Chromosome name`=="Chom_Scaffold_315", "X", `Chromosome name`),
         `Chromosome name`=ifelse(`Chromosome name`=="Chom_Scaffold_159", "2", `Chromosome name`),
         `Chromosome name`=ifelse(`Chromosome name`=="Chom_Scaffold_498", "3", `Chromosome name`),
         `Chromosome name`=ifelse(`Chromosome name`=="Chom_Scaffold_470", "4", `Chromosome name`),
         `Chromosome name`=ifelse(`Chromosome name`=="Chom_Scaffold_416", "5", `Chromosome name`),
         `Chromosome name`=ifelse(`Chromosome name`=="Chom_Scaffold_497", "6", `Chromosome name`))

comp_ows_chom <- comp_ows_chom %>%
  mutate(`Muller element` = case_when(
    `Chromosome name` == "2" ~ "B",
    `Chromosome name` == "X" ~ "F",
    `Chromosome name` == "3" ~ "A",
    `Chromosome name` == "4" ~ "E",
    `Chromosome name` == "5" ~ "D",
    `Chromosome name` == "6" ~ "C"))

# writing tables
fwrite(x = comp_ows_chom, quote = FALSE, sep = '\t', row.names = FALSE,
       col.names = FALSE, file = "ows_chom.txt")
