# Libraries
library("ggplot2")
library("ggpubfigs")
library("data.table")
library("gridExtra")
library("ggtext")

# data

# lengths of contigs
# awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' cbez_original_filtered.fsa > scaff_len.txt
len <- fread("scaff_len.txt")
sort_len <- sort(len$V1, decreasing = TRUE)
cum_len <- cumsum(sort_len/1e6)
num_contigs <- seq_along(cum_len)
df2 <- data.table(num_contigs = num_contigs, cum_len = cum_len)


# BUSCO values for Cbez
busco <- rep(c("Pre-polished Genome Assembly", "Post-polished Genome Assembly", "Final Genome Assembly","Proteome", "Transcriptome"), 
             each = 4)
gene_cat <- rep(c("Complete and Single-Copy", 
                         "Complete and Duplicated", 
                         "Fragmented",
                         "Missing"), 3)
per <- c(93.8, 0.9, 2, 3.3,
         98.4, 0.7, 0.4, 0.5,
         98.4, 0.5, 0.4, 0.7,
         90.7, 6.5, 1.3, 1.5,
         43.2, 38.8, 4.1, 13.9)

df <- data.table(busco, gene_cat, per)

# BUSCO values for other Calliphoridae
busco2 <- rep(c("Chrysomya bezziana", "Calliphora vomitoria", "Calliphora grahami",
               "Bellardia pandia", "Cochliomyia hominivorax", "Protocalliphora azurea",
               "Lucilia cuprina", "Lucilia sericata", "Chrysomya rufifacies", "Phormia regina",
               "Calliphora vicina"), 
             each = 4)
gene_cat2 <- rep(c("Complete and Single-Copy", 
                   "Complete and Duplicated", 
                   "Fragmented",
                   "Missing"), 3)
per2 <- c(98.4, 0.5, 0.4, 0.7,
          98.2, 0.5, 0.3, 1,
          98.2, 0.4, 0.7, 0.7,
          97.4, 1.1, 0.3, 1.2,
          96, 2.1, 0.5, 1.4,
          97.5, 0.2, 0.2, 1.6,
          95.6, 0.6, 1.6, 2.2,
          92, 3.1, 1.5, 3.4,
          92.9, 0.5, 3.6, 3,
          88.2, 1.6, 5.6, 4.6,
          60.1, 0.2, 17.3, 22.4)
df3 <- data.table(busco2, gene_cat2, per2)
df3[, busco2 := factor(busco2, levels = unique(busco2))] # to plot in the same order of busco2 and not alphabetically

# plots

# cumulative length plot from Cbez
p1 <- ggplot(df2, aes(x = num_contigs, y = cum_len)) +
  geom_line(linewidth = 1.5) +
  scale_y_continuous(limits = c(0, max(df2$cum_len))) +
  scale_x_continuous(limits = c(0, max(df2$num_contigs)), breaks = c(0, 200, 400, 600, 800, 991)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  labs(x = "Contig Index",
       y = "Cumulative Length (Mbp)", 
       title = "Cumulative Length") +
  theme_classic()

# BUSCO - genome, proteome (BRAKER), trascriptome (trinity) from Cbez
p2 <- ggplot(data = df, aes(x = busco, y = per, fill = gene_cat)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) + # color blind friendly
  labs(y = "Percentage of BUSCO genes from diptera_odb10 (n = 2385)") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12)) +
  guides(fill = guide_legend(title = "Gene Category")) +
  theme_classic()

# Nx plot from Cbez
len <- fread("scaff_len.txt")
len_mbp <- len$V1/1e6
sort_len <- sort(len_mbp, decreasing = TRUE)

gen_len <- sum(sort_len)  # total genome length
cum_len <- cumsum(sort_len)

x_values <- 1:100  # N1 to N100
cont_len_nx <- numeric(100)  # calculate the Nx values

# Nx values
for (i in x_values) {
  val <- (gen_len*i)/100  # cumulative length for Nx
  idx <- which(cum_len >= val)[1]
  cont_len_nx[i] <- ifelse(!is.na(idx), sort_len[idx], NA)
}

# Create data table
nx_values <- data.table(x = x_values, contig_length = cont_len_nx)

N50 <- sum(sort_len)/2
N50 <- sort_len[min(which(cum_len >= N50))]

# Nx Plot
p3 <- ggplot(nx_values, aes(x = x, y = contig_length)) +
  geom_line(size = 1.5) +
  scale_x_continuous(limits = c(1, 100), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  geom_hline(yintercept = N50, linetype = "dashed", size = 0.5, color = "gray") +
  labs(title = "Nx Plot",
       x = "Nx", 
       y = "Contig Length (Mbp)") +
  theme_minimal()

# BUSCO for blowfly species
p4 <- ggplot(data = df3, aes(x = busco2, y = per2, fill = gene_cat2)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  scale_fill_manual(values = friendly_pal("ito_seven")) + # Color-blind friendly
  scale_x_discrete(labels = function(x) {
    ifelse(x == "Chrysomya bezziana", expression(bold("Chrysomya bezziana")), x)
  }) +
  labs(y = "Percentage of BUSCO genes from diptera_odb10 (n = 2385)") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12)) +
  guides(fill = guide_legend(title = "Gene Category")) +
  theme_classic()


grid.arrange(p2, p3, p1, p4, ncol = 2)
