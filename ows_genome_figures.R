# Libraries
library("ggplot2")
library("ggpubfigs")
library("data.table")
library("gridExtra")

# data
# NCBI table
df <- fread("ncbi_diptera.tsv")
df <- na.omit(df)
df <- df[grepl("Nanopore", df$`Assembly Sequencing Tech`),]
df <- df[!duplicated(df$`Organism Name`)]

nr <- list(`Organism Name` = "Chrysomya bezziana",
           `Assembly Stats Number of Scaffolds` = 992,
           `Assembly Stats Scaffold N50` = 1851041)

df <- rbind(df, nr, fill = TRUE)

# length for density plot
# awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' cbez_N_genome_final.fa > scaff_len.txt
len <- fread("scaff_len.txt")

# BUSCO
busco <- rep(c("Genome", "Proteome", "Transcriptome"), each = 4)
`Gene Category` <- rep(c("Complete and Single-Copy", 
                         "Complete and Duplicated", 
                         "Fragmented",
                         "Missing"), 3)
per <- c(98.5, 0.5, 0.4, 0.6,
         90.7, 6.5, 1.3, 1.5,
         43.2, 38.8, 4.1, 13.9)

df2 <- data.table(busco, cat, per)

# plots
# density of scaffold lengths (fig. 1 a)
m <- mean(len$V1)

p1 <- ggplot(len, aes(x = V1)) +
  geom_density() +
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  labs(x = "Sequence length (bp)",
       y = "Density") +
  geom_vline(aes(xintercept = m),
             color = "gray", linetype = "dashed") +
  #annotate("text", hjust = 0, col = "black",
   #        label=(paste0("Mean length = ", round(m)))) +
theme_classic()

# N50 comparasion (fig. 1 b)
p2 <- ggplot(df, aes(x = `Organism Name`, y = `Assembly Stats Scaffold N50`, 
               color = `Organism Name`)) + 
  geom_jitter(width = 1, size = 3, 
              aes(color = ifelse(`Organism Name` == "Chrysomya bezziana", 
                                                      "Chrysomya bezziana", 
                                                      "Other Diptera"))) + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  scale_x_discrete(expand = c(0.5, 0.5)) + 
  scale_color_manual(values = c("Other Diptera" = "grey", 
                                "Chrysomya bezziana" = "#87aff5")) +
  guides(fill = "none", color = "none")

# Number of scaffolds comparasion (fig. 1 c)
p3 <- ggplot(df, aes(x = `Organism Name`, y = `Assembly Stats Number of Scaffolds`, 
               color = `Organism Name`)) + 
  geom_jitter(width = 1, size = 3, 
              aes(color = ifelse(`Organism Name` == "Chrysomya bezziana", 
                                                      "Chrysomya bezziana", 
                                                      "Other Diptera"))) + 
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  scale_x_discrete(expand = c(0.5, 0.5)) + 
  scale_color_manual(values = c("Other Diptera" = "grey", 
                                "Chrysomya bezziana" = "#87aff5")) +
  guides(fill = "none", color = "none")

# BUSCO - genome, proteome (BRAKER), trascriptome (trinity) (fig. 1 d)
p4 <- ggplot(data = df2, aes(x = busco, y = per, fill = `Gene Category`)) +
  geom_bar(stat = "identity", alpha = 1/3) +
  scale_fill_manual(values = friendly_pal("wong_eight")) + # color blind friendly
  labs(y = "Percentage of BUSCO genes from diptera_odb10") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12)) +
  theme_classic()

grid.arrange(p1, p2, p3, p4, ncol = 2)
