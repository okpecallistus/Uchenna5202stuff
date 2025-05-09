#STEP 1: LOAD REQUIRED LIBRARIES
install.packages(c("phyloseq", "vegan", "ggplot2", "tidyverse"))
install.packages("phyloseq", dependencies = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
a
library(phyloseq)


library(phyloseq)
library(vegan)
library(ggplot2)
library(tidyverse)

#Step 2 Read Data
hh_phylum<-read.csv("https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/HH_combined_bracken_phylum_fraction.csv")
View(hh_phylum)
hh_species<-read.csv("https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/HH_combined_bracken_species_fraction.csv")
View(hh_species)
pC_phylim<-read.csv("https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/PC_combined_bracken_phylum_fraction.csv")
View(pC_phylim)
pC_species<-read.csv("https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/PC_combined_bracken_species_fraction.csv")
View(pC_species)

#STEP 3 COnvert to Phyloseq Object

otu_matrix <- as.matrix(hh_phylum)
OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# Optional: if you have metadata, you can load it too
# sample_data <- read.csv("your_metadata.csv", row.names = 1)
# SAMPLE <- sample_data(sample_data)

physeq_obj <- phyloseq(OTU)  # Add SAMPLE if available: phyloseq(OTU, SAMPLE)

# Step 1: Load the CSV again with row names set correctly
hh_phylum <- read.csv("https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/HH_combined_bracken_phylum_fraction.csv", 
                      row.names = 1, check.names = FALSE)

# Step 2: Confirm the data is numeric
str(hh_phylum)  # This should show: 'num' or 'num [1:XX, 1:YY]' for all columns

# Step 3: Convert to matrix
otu_matrix <- as.matrix(hh_phylum)

# Step 4: Create the phyloseq OTU table
OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)

#Step 4: Alpha Diversity

# Calculate common alpha diversity metrics
alpha_div <- estimate_richness(physeq_obj, measures = c("Shannon", "Simpson"))
head(alpha_div)

# Plot Shannon Diversity
library(ggplot2)
alpha_div$SampleID <- rownames(alpha_div)
ggplot(alpha_div, aes(x = SampleID, y = Shannon)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Shannon Diversity Index per Sample",
       x = "Sample ID", y = "Shannon Index")

#Plot Simpson Index
ggplot(alpha_div, aes(x = SampleID, y = Simpson)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Simpson Diversity Index per Sample",
       x = "Sample ID", y = "Simpson Index")

colnames(hh_phylum)

#Create Metadata
#install.packages(c("vegan", "ggplot2", "readr"))  # only once
library(vegan)
library(ggplot2)
library(readr)
head(hh_species[1:5])

#Transpose Data (This is because most functions(like veagn::diversity()))expect samples as rows, and Features (species) as columns.
hh_species_t<-t(hh_species)
dim(hh_species_t) #check shape
head(hh_species_t[,1:5])

#Fix to ensure that everything is numeric
# Convert all columns to numeric safely
hh_species_t_numeric <- hh_species_t #make a copy
view(hh_species_t_numeric)
#Loop over columns to ensure numeric conversion
for(col in colnames(hh_species_t_numeric)) {hh_species_t_numeric[,col]<-as.numeric(hh_species_t_numeric[,col])}
#Check column and row names
head(colnames(hh_species_t_numeric)) #species name
head(rownames(hh_species_t_numeric)) #sample IDs
# Preserve sample names as row names
rownames(hh_species_t_numeric) <- rownames(hh_species_t)

str(hh_species_t_numeric)


 #Step 4: Calculate Shannon and Simpson Diversity using vegan::diversity()
shannon_index<-diversity(hh_species_t, index = "shannon")
simpson_index<-diversity(hh_species_t, index = "simpson")
 
#Starting From the Beggenning
hh_species<-read.csv("https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/HH_combined_bracken_species_fraction.csv", row.names = 1, check.names = FALSE)
view(hh_species)

#Transpose
hh_species_t<-t(hh_species)
view(hh_species_t)
#Now:
#rownames(hh_species_t) = sample IDs
#colnames(hh_species_t) = species names ✅
#Convert to Numeric Safely without losing names
#Creat a new data frame, retaining column and row names
hh_species_t_numeric<-as.data.frame(hh_species_t)
#use lapply to convert each column to numeric while preserving names
hh_species_t_numeric<-hh_species_t_numeric%>%
  mutate(across(everything(), ~as.numeric(.)))
#check that species names are still in column names
head(colnames(hh_species_t_numeric))
#Check Data Structure
str(hh_species_t_numeric)

#Lets compute Diversity indices
library(vegan)
shannon_index_species<-diversity(hh_species_t_numeric,index = "shannon")
simpson_index_species<-diversity(hh_species_t_numeric, index = "simpson")

#Prepare Data for Plotting
diversity_df <- data.frame(
  SampleID = rownames(hh_species_t_numeric),
  Shannon = shannon_index_species,
  Simpson = simpson_index_species
)

#Convert this long Plot format for easier plotting:
library(tidyr)

diversity_long_species <- pivot_longer(
  diversity_df,
  cols = c("Shannon", "Simpson"),
  names_to = "Index",
  values_to = "Value"
)
view(diversity_long_species)

#Plot a Boxplot with Jitter
library(ggplot2)

p <- ggplot(diversity_long, aes(x = Index, y = Value, fill = Index)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.6, color = "black", size = 1) +
  theme_minimal() +
  labs(
    title = "Alpha Diversity in Healthy Individuals",
    x = "Diversity Index",
    y = "Value"
  ) +
  scale_fill_manual(values = c("Shannon" = "#66c2a5", "Simpson" = "#fc8d62"))
print(p)

#Box Plot 2
library(ggplot2)

p <- ggplot(diversity_long, aes(x = Index, y = Value, fill = Index)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.6, color = "black", size = 1) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  labs(
    title = "Alpha Diversity in Healthy Individuals",
    x = "Diversity Index",
    y = "Diversity Score"
  ) +
  scale_fill_manual(values = c("Shannon" = "#66c2a5", "Simpson" = "#fc8d62"))

# Show the plot
print(p)
ggsave("diversity_boxplot.png", plot = p, width = 8, height = 6, dpi = 300)

#Alpha Diversity Calculations Computation for pC

pC_species<-read.csv("https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/PC_combined_bracken_species_fraction.csv", row.names = 1, check.names = FALSE)
view(pC_species)
#Transpose
pC_species_t<-t(pC_species)
view(pC_species_t)
#Rownames(pC_species_t) =sample ID
#Colnames(pC_species_t) = species name
#Convert to numeric safely without losing names
#Create a new data frame, retaining column and row names
pC_species_t_numeric<-as.data.frame(pC_species_t)
view(pC_species_t_numeric)
#Use lapply to convert each column to numeric while preserving names
pC_species_t_numeric<-pC_species_t_numeric%>%
  mutate(across(everything(),~as.numeric(.)))
#Check that species names are still in column names
head(colnames(pC_species_t_numeric))
#Check Data Structure
str(pC_species_t_numeric)
#Lets compute DIversity Indices
library(vegan)
pC_shannon_index_sp<-diversity(pC_species_t_numeric, index = "shannon")
pC_simpson_index_sp<-diversity(pC_species_t_numeric, index = "simpson")
#Prepare Data for Plotting
pC_diversity_df<-data.frame(SampleID = rownames(pC_species_t_numeric),Shannon = pC_shannon_index_sp, Simpson = pC_simpson_index_sp)
view(pC_diversity_df)
#Convert this long plot format for easier plotting:
library(tidyr)
pC_diversity_long_sp<-pivot_longer(pC_diversity_df, cols = c("Shannon", "Simpson"), names_to = "Index", values_to ="Value")
view(pC_diversity_long_sp)
#Box Plot
library(ggplot2)
pC_sp <- ggplot(pC_diversity_long_sp, aes(x = Index, y = Value, fill = Index)) +
  geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, color = "black", size = 1) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(
    title = "Alpha Diversity (Shannon and Simpson) in Cancer patients Microbiome",
    x = "Diversity Index",
    y = "Diversity Score"
  ) +
  scale_fill_manual(values = c("Shannon" = "#E69F00", "Simpson" = "#56B4E9"))  # color-blind friendly
print(pC_sp)
ggsave("alpha_diversity_boxplot.png", plot = pC_sp, width = 8, height = 6, dpi = 300)



#Comparative Diversity Plots for HH & pC
install.packages(c("readxl", "vegan", "ggplot2", "dplyr", "ggpubr"))
library(readxl)
library(vegan)
library(ggplot2)
library(dplyr)
library(ggpubr)
view(diversity_df)
view(pC_diversity_df)
#Group Data
pC_diversity_df$Group<-"Cancer"
diversity_df$Group<-"Healthy"
#Combine Both Datasets
combined_div<-rbind(diversity_df, pC_diversity_df)
view(combined_div)
#Plot: Simpson BOXplot with Statistical Test
library(ggplot2)
library(ggpubr)

p_simpson <- ggplot(combined_div, aes(x = Group, y = Simpson, fill = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1, color = "black") +
  
  # Statistical test with bar and asterisks
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("Healthy", "Cancer")), 
                     tip.length = 0.02, size = 5) +
  
  # Apply a clean theme and add gridlines + border
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.position = "none"
  ) +
  
  labs(
    title = "Simpson Diversity: Healthy vs Cancer",
    x = "Group",
    y = "Simpson Index"
  ) +
  
  # Color-blind friendly & elegant colors (Okabe-Ito palette)
  scale_fill_manual(values = c("Healthy" = "#009E73", "Cancer" = "#D55E00"))
print(p_simpson)
ggsave("combined_Simpson_Diversity_Boxplot.png", plot = p_simpson, width = 8, height = 6, dpi = 300)

# Wilcoxon rank-sum test (non-parametric)
wilcox.test(Simpson ~ Group, data = combined_div)

#Shannon Box Plot
p_shannon <- ggplot(combined_div, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1, color = "black") +
  
  # Statistical test with bar and asterisks
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("Healthy", "Cancer")), 
                     tip.length = 0.02, size = 5) +
  
  # Apply a clean theme and add gridlines + border
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.position = "none"
  ) +
  
  labs(
    title = "Shannon Diversity: Healthy vs Cancer",
    x = "Group",
    y = "Shannon Index"
  ) +
  
  # Color-blind friendly & elegant colors (Okabe-Ito palette)
  scale_fill_manual(values = c("Healthy" = "#56B4E9", "Cancer" = "#E69F00"))
print(p_shannon)
ggsave("combined_Shannon_Diversity_Boxplot.png", plot = p_shannon, width = 8, height = 6, dpi = 300)

# Wilcoxon rank-sum test (non-parametric)
wilcox.test(Shannon ~ Group, data = combined_div)

#Permanova calculations
# Create a metadata data frame
metadata_df <- data.frame(Group = group_labels)

# Run PERMANOVA
adonis_result <- adonis2(bray_dist ~ Group, data = metadata_df, permutations = 999)

# View result
print(adonis_result)

#PART B Taxanomic COmparison
#Q: Which Bacteria are more Abundant in Healthy VS Cancer?
#Option A: Wilcoxon test (for small datasets)
#Combine data
install.packages(c("readxl", "dplyr", "ggplot2", "VennDiagram", "indicspecies"))
install.packages(c("VennDiagram", "indicspecies"))
library(readxl)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(indicspecies)

#STEP 2: Read Excel Data (Relative Abundance Tables)
library(readr)

# Re-import correctly
# Load necessary library
library(readr)

# Load required library
library(ggplot2)

cat("=== Step 1: Reading and cleaning Healthy data ===\n")
h_raw <- read.csv(
  "https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/HH_combined_bracken_species_fraction.csv",
  row.names = 1, check.names = FALSE
)
h_raw_clean <- h_raw[!rownames(h_raw) %in% c("Homo sapiens"), ]
h_t_clean <- as.data.frame(t(h_raw_clean))
h_t_clean[] <- lapply(h_t_clean, function(x) as.numeric(as.character(x)))
h_t_clean$Group <- "Healthy"
cat("✅ Healthy dataset cleaned and transposed.\n")

cat("=== Step 2: Reading and cleaning Cancer data ===\n")
p_raw <- read.csv(
  "https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/PC_combined_bracken_species_fraction.csv",
  row.names = 1, check.names = FALSE
)
p_raw_clean <- p_raw[!rownames(p_raw) %in% c("Homo sapiens"), ]
p_t_clean <- as.data.frame(t(p_raw_clean))
p_t_clean[] <- lapply(p_t_clean, function(x) as.numeric(as.character(x)))
p_t_clean$Group <- "Cancer"
cat("✅ Cancer dataset cleaned and transposed.\n")

cat("=== Step 3: Harmonizing and merging datasets ===\n")
all_species <- union(colnames(h_t_clean), colnames(p_t_clean))
missing_in_hh <- setdiff(all_species, colnames(h_t_clean))
missing_in_pc <- setdiff(all_species, colnames(p_t_clean))
h_t_clean[missing_in_hh] <- 0
p_t_clean[missing_in_pc] <- 0
h_t_clean <- h_t_clean[, all_species]
p_t_clean <- p_t_clean[, all_species]
combined_data <- rbind(h_t_clean, p_t_clean)
cat("✅ Combined dataset created. Structure:\n")
print(str(combined_data))

cat("=== Step 4: Preparing for Wilcoxon testing ===\n")
group_labels <- combined_data$Group
abundance_data <- combined_data[, !colnames(combined_data) %in% "Group"]
cat("✅ Confirm Homo sapiens excluded: ", "Homo sapiens" %in% colnames(abundance_data), "\n")

cat("=== Step 5: Running Wilcoxon tests ===\n")
pvals <- apply(abundance_data, 2, function(x) {
  wilcox.test(x ~ group_labels)$p.value
})
padj <- p.adjust(pvals, method = "fdr")
results_df <- data.frame(
  Species = names(pvals),
  p_value = pvals,
  p_adj = padj
)
cat("✅ Wilcoxon test completed. Top species:\n")
print(head(results_df[order(results_df$p_adj), ]))

cat("=== Step 6: Plotting top 5 species ===\n")

top_species <- results_df$Species[order(results_df$p_adj)][1:5]

#Create a function to plot and save each boxplot

# Load ggplot2 if not already loaded
library(ggplot2)

# Loop through top species and save each plot
for (sp in top_species) {
  p <- ggplot(combined_data, aes(x = Group, y = .data[[sp]], fill = Group)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste("Abundance of", sp), y = "Relative Abundance") +
    scale_fill_manual(values = c("Healthy" = "#56B4E9", "Cancer" = "#E69F00"))
  
  # Generate a safe filename
  fname <- paste0("boxplot_", gsub("[^a-zA-Z0-9]", "_", sp), ".png")
  
  # Save the plot
  ggsave(filename = fname, plot = p, width = 6, height = 4, dpi = 300)
  
  # Optional: print confirmation
  cat("✅ Saved:", fname, "\n")
}

#Grouped Boxplot by Species and Group
# Load necessary libraries
#library(ggplot2)
#library(tidyr)
#library(dplyr)

# Select top species
top_species <- results_df$Species[order(results_df$p_adj)][1:5]

# Subset combined data to top species + Group
abundance_subset <- combined_data[, c(top_species, "Group")]

# Reshape to long format for ggplot
abundance_long <- pivot_longer(
  abundance_subset,
  cols = all_of(top_species),
  names_to = "Species",
  values_to = "Abundance"
)

# Plot as grouped boxplot
ggplot(abundance_long, aes(x = Species, y = Abundance, fill = Group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  theme_minimal() +
  labs(title = "Relative Abundance of Top 5 Species", y = "Relative Abundance") +
  scale_fill_manual(values = c("Healthy" = "#56B4E9", "Cancer" = "#E69F00")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Enhancement of my Graphs
#install.packages("ggpubr")  # Only once
library(ggpubr)


# Load necessary libraries
#library(ggplot2)
#library(tidyr)
#library(dplyr)
#library(ggpubr)

# Select top species
top_species <- results_df$Species[order(results_df$p_adj)][1:5]

# Subset and reshape
abundance_subset <- combined_data[, c(top_species, "Group")]

abundance_long <- pivot_longer(
  abundance_subset,
  cols = all_of(top_species),
  names_to = "Species",
  values_to = "Abundance"
)

# Plot with statistical comparison
ggplot(abundance_long, aes(x = Species, y = Abundance, fill = Group)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) +
  stat_compare_means(
    aes(group = Group),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  scale_y_continuous(trans = "log10") +  # log scale to handle low values
  theme_minimal() +
  labs(
    title = "Relative Abundance of Top 5 Species",
    y = "Relative Abundance (log10 scale)",
    x = "Species"
  ) +
  scale_fill_manual(values = c("Healthy" = "#56B4E9", "Cancer" = "#E69F00")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Load necessary libraries
#library(ggplot2)
#library(tidyr)
#library(dplyr)
#library(ggpubr)

# Step 1: Select top 5 species
top_species <- results_df$Species[order(results_df$p_adj)][1:5]

# Step 2: Subset and reshape data
abundance_subset <- combined_data[, c(top_species, "Group")]

abundance_long <- pivot_longer(
  abundance_subset,
  cols = all_of(top_species),
  names_to = "Species",
  values_to = "Abundance"
)

# Step 3: Plot with enclosed axes and stat bar
p <- ggplot(abundance_long, aes(x = Species, y = Abundance, fill = Group)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) +
  stat_compare_means(
    aes(group = Group),
    method = "wilcox.test",
    label = "p.signif",
    label.y.npc = "top",
    bracket.size = 0.6,
    tip.length = 0.02,
    size = 4,
    hide.ns = TRUE
  ) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c("Healthy" = "#56B4E9", "Cancer" = "#E69F00")) +
  labs(
    title = "Differential Abundance of Top 5 Species",
    x = "Species",
    y = "Relative Abundance (log10 scale)"
  ) +
  theme_classic(base_size = 13) +  # Use classic theme for boxed plot
  theme(
    axis.line = element_line(color = "black", size = 0.8),
    panel.border = element_rect(fill = NA, color = "black", size = 1),  # Full box around plot
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",       # Legend on the right
    legend.title = element_blank()
  )

# Step 4: Print to viewer
print(p)

# Step 5: Save to PNG
ggsave("top5_species_boxplot.png", plot = p, width = 10, height = 6, dpi = 300)
cat("✅ Plot saved as 'top5_species_boxplot.png'\n")



#Venn DIagram (Shared Species)
library(VennDiagram)

# Identify species with non-zero counts
hh_species <- colnames(h_t_clean)[colSums(h_t_clean[, -ncol(h_t_clean)]) > 0]
pc_species <- colnames(p_t_clean)[colSums(p_t_clean[, -ncol(p_t_clean)]) > 0]

venn.plot <- venn.diagram(
  x = list(Healthy = hh_species, Cancer = pc_species),
  category.names = c("Healthy", "Cancer"),
  filename = "venn_species_overlap.png",
  output = TRUE
)

# Identify only numeric columns
numeric_cols <- sapply(h_t_clean, is.numeric)

# Get species with non-zero abundance (only numeric columns)
hh_species <- colnames(h_t_clean)[numeric_cols & colSums(h_t_clean[, numeric_cols], na.rm = TRUE) > 0]


# Step 1: Identify numeric columns only
numeric_cols <- sapply(h_t_clean, is.numeric)

# Step 2: Compute column sums and select non-zero columns from the numeric ones
non_zero_numeric <- colSums(h_t_clean[, numeric_cols], na.rm = TRUE) > 0

# Step 3: Get the names of those species (i.e., columns)
hh_species <- colnames(h_t_clean)[numeric_cols][non_zero_numeric]


#Venn DIagram (Shared Species)
library(VennDiagram)

# Identify species with non-zero counts
hh_species <- colnames(h_t_clean)[colSums(h_t_clean[, -ncol(h_t_clean)]) > 0]
pc_species <- colnames(p_t_clean)[colSums(p_t_clean[, -ncol(p_t_clean)]) > 0]

venn.plot <- venn.diagram(
  x = list(Healthy = hh_species, Cancer = pc_species),
  category.names = c("Healthy", "Cancer"),
  filename = "venn_species_overlap.png",
  output = TRUE
)

library(VennDiagram)

# Step 1: Identify numeric columns only
hh_numeric_cols <- sapply(h_t_clean, is.numeric)
pc_numeric_cols <- sapply(p_t_clean, is.numeric)

# Step 2: Get species (columns) with non-zero abundance
hh_species <- colnames(h_t_clean)[hh_numeric_cols][colSums(h_t_clean[, hh_numeric_cols]) > 0]
pc_species <- colnames(p_t_clean)[pc_numeric_cols][colSums(p_t_clean[, pc_numeric_cols]) > 0]

# Step 3: Generate and save Venn diagram
venn.plot <- venn.diagram(
  x = list(Healthy = hh_species, Cancer = pc_species),
  category.names = c("Healthy", "Cancer"),
  filename = "venn_species_overlap.png",  # Saves to current working directory
  output = TRUE,
  imagetype = "png",
  height = 2000,
  width = 2000,
  resolution = 300,
  col = "black",
  fill = c("#56B4E9", "#E69F00"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.pos = 0
)
print(venn.plot)


#Heatmao Code:
library(pheatmap)

# Choose top 20 species based on lowest adjusted p-values
top_heatmap_species <- results_df$Species[order(results_df$p_adj)][1:20]

# Prepare abundance matrix for heatmap (samples as rows, species as columns)
heatmap_data <- combined_data[, top_heatmap_species]
rownames(heatmap_data) <- paste0(combined_data$Group, "_", seq_len(nrow(combined_data)))

# Group annotation (used for coloring rows by Healthy or Cancer)
annotation <- data.frame(Group = combined_data$Group)
rownames(annotation) <- rownames(heatmap_data)

# Plot heatmap with species names shown (they are columns)
pheatmap(
  t(heatmap_data),                # Transpose so species are on y-axis
  annotation_col = annotation,   # Now columns are samples
  scale = "row",                 # Normalize species abundance per row
  show_rownames = TRUE,          # Show species names
  show_colnames = FALSE,         # Hide sample names if too many
  fontsize_row = 8,              # Adjust if species names are too long
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Heatmap of Top 20 Differential Species"
)
#ALternative HEatmap savings
pheatmap(
  t(heatmap_data),
  annotation_col = annotation,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Heatmap of Top 20 Differential Species",
  filename = "top20_species_heatmap.png",  # 👈 Saves as PNG
  width = 8, height = 6                    # Inches; adjust as needed
)

#Color blind friendly heatmap




library(pheatmap)
library(viridis)

pheatmap(
  t(heatmap_data),
  annotation_col = annotation,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color = viridis(100, option = "D"),  # Use viridis palette
  main = "Heatmap of Top 20 Differential Species",
  filename = "top20_species_heatmap_colorblind_friendly.png",
  width = 8, height = 6
)


  #Full R Code Block for Indicator bacteria for both species
# Load required libraries
library(indicspecies)
library(vegan)
library(ggplot2)

# Prepare abundance matrix (numeric columns only)
abundance_matrix <- combined_data[, sapply(combined_data, is.numeric)]
group_factor <- combined_data$Group

# Run indicator species analysis
set.seed(123)  # Reproducibility
ind_result <- multipatt(abundance_matrix, group_factor, control = how(nperm = 999))

# Summarize result (optional to display)
summary(ind_result, indvalcomp = TRUE, alpha = 0.05)

# Extract results and get top 20 by p-value
ind_df <- as.data.frame(ind_result$sign)
ind_df$Species <- rownames(ind_df)
top_indicators <- ind_df[order(ind_df$p.value), ][1:20, ]

# Determine correct column names (these start with 's.' and indicate group membership)
s_columns <- grep("^s\\.", colnames(ind_df), value = TRUE)

# If there are only two groups, map them safely:
if (length(s_columns) == 2) {
  group_labels <- c("Healthy", "Cancer")  # Adjust if your group names differ
  
  # Assign group labels
  top_indicators$Group <- apply(top_indicators[, s_columns], 1, function(row) {
    assigned <- which(row == 1)
    if (length(assigned) == 1) {
      return(group_labels[assigned])
    } else {
      return("Both")
    }
  })
} else {
  stop("More than two groups detected — please adjust group labeling.")
}


# Add group labels
top_indicators$Group <- ifelse(top_indicators$s.1 == 1, "Healthy", 
                               ifelse(top_indicators$s.2 == 1, "Cancer", "Both"))

# Plot top 20 indicator species
plot <- ggplot(top_indicators, aes(x = reorder(Species, -stat), y = stat, fill = Group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top 20 Indicator Species by Group",
    y = "Indicator Value",
    x = "Species"
  ) +
  scale_fill_manual(values = c("Healthy" = "#56B4E9", "Cancer" = "#E69F00", "Both" = "gray")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

# Display plot
print(plot)

# Save plot to PNG
ggsave("top_20_indicator_species.png", plot = plot, width = 10, height = 6, dpi = 300)
cat("✅ Plot saved as top_20_indicator_species.png\n")





# Load libraries
library(ggplot2)
library(indicspecies)
library(dplyr)

# Step 1: Prepare data (assuming `combined_data` exists with Group column and numeric species columns)
abundance_matrix <- combined_data[, -ncol(combined_data)]
group_factor <- combined_data$Group
# Exclude Group column and ensure all remaining are numeric
abundance_matrix <- combined_data[, sapply(combined_data, is.numeric)]


# Step 2: Run indicator species analysis
group_factor <- combined_data$Group
ind_result <- multipatt(abundance_matrix, group_factor, control = how(nperm = 999))

# Step 3: Extract results
ind_df <- as.data.frame(ind_result$sign)
ind_df$Species <- rownames(ind_df)

# Step 4: Determine group assignments (safely detect columns that start with "s.")
s_columns <- grep("^s\\.", colnames(ind_df), value = TRUE)
group_labels <- levels(group_factor)

# Assign group labels based on s.1 and s.2 (or more generally)
ind_df$Group <- apply(ind_df[, s_columns], 1, function(row) {
  assigned <- which(row == 1)
  if (length(assigned) == 1) {
    return(group_labels[assigned])
  } else {
    return("Both")
  }
})

# Step 5: Filter significant results and get top 30 for each group
top_healthy <- ind_df %>%
  filter(Group == "Healthy", p.value <= 0.05) %>%
  arrange(desc(stat)) %>%
  slice_head(n = 30)

top_cancer <- ind_df %>%
  filter(Group == "Cancer", p.value <= 0.05) %>%
  arrange(desc(stat)) %>%
  slice_head(n = 30)

# Step 6: Plot Healthy
p_healthy <- ggplot(top_healthy, aes(x = reorder(Species, stat), y = stat)) +
  geom_bar(stat = "identity", fill = "#56B4E9") +
  coord_flip() +
  labs(title = "Top 30 Indicator Species in Healthy Group", x = "Species", y = "Indicator Value") +
  theme_minimal(base_size = 12)

# Step 7: Plot Cancer
p_cancer <- ggplot(top_cancer, aes(x = reorder(Species, stat), y = stat)) +
  geom_bar(stat = "identity", fill = "#E69F00") +
  coord_flip() +
  labs(title = "Top 30 Indicator Species in Cancer Group", x = "Species", y = "Indicator Value") +
  theme_minimal(base_size = 12)

# Step 8: Display and Save
print(p_healthy)
ggsave("top_30_indicator_healthy.png", p_healthy, width = 10, height = 8, dpi = 300)

print(p_cancer)
ggsave("top_30_indicator_cancer.png", p_cancer, width = 10, height = 8, dpi = 300)

summary(ind_result)
colnames(ind_result$sign)


library(dplyr)

# Convert result to data frame
ind_df <- as.data.frame(ind_result$sign)

# Filter top 30 by stat value (or you could use p.value if available)
top_healthy <- ind_df %>%
  filter(s.Healthy == 1) %>%
  slice_max(stat, n = 30) %>%
  mutate(Group = "Healthy")

top_cancer <- ind_df %>%
  filter(s.Cancer == 1) %>%
  slice_max(stat, n = 30) %>%
  mutate(Group = "Cancer")

#Subset the Abundance data from top indicator lists
top_healthy_species <- rownames(top_healthy)
top_cancer_species <- rownames(top_cancer)

# Combine for plotting
top_combined <- bind_rows(top_healthy, top_cancer)
top_combined$Species <- rownames(top_combined)

# Get only the species columns
healthy_abundances <- healthy_data[, top_healthy_species]
cancer_abundances <- cancer_data[, top_cancer_species]




#Top 30 Gut microbiome 
# Load libraries
library(dplyr)
library(ggplot2)

# Step 1: Prepare the abundance matrix and group info
group_factor <- combined_data$Group

# ✅ Remove non-numeric columns for abundance matrix
abundance_matrix <- combined_data[, sapply(combined_data, is.numeric)]

# Step 2: Run indicator species analysis
library(indicspecies)
library(vegan)
abundance_matrix <- combined_data[, -ncol(combined_data)]  # assuming last column is Group
ind_result <- multipatt(abundance_matrix, group_factor, control = how(nperm = 999))
group_factor <- combined_data$Group

# Step 2: Run indicator species analysis
library(indicspecies)
library(vegan)
ind_result <- multipatt(abundance_matrix, group_factor, control = how(nperm = 999))

# Step 3: Extract significant species
ind_df <- as.data.frame(ind_result$sign)
ind_df$Species <- rownames(ind_df)

# Step 4: Identify top 30 species per group based on indicator statistic
top_healthy <- ind_df %>%
  filter(s.Healthy == 1) %>%
  slice_max(stat, n = 30)

top_cancer <- ind_df %>%
  filter(s.Cancer == 1) %>%
  slice_max(stat, n = 30)

# Step 5: Get species names
top_healthy_species <- top_healthy$Species
top_cancer_species <- top_cancer$Species

# Step 6: Subset original data by group
healthy_data <- combined_data %>% filter(Group == "Healthy")
cancer_data  <- combined_data %>% filter(Group == "Cancer")

# Step 7: Extract abundance values for top species
healthy_abundances <- healthy_data[, top_healthy_species]
cancer_abundances  <- cancer_data[, top_cancer_species]

# Step 8: Compute mean abundance per species
healthy_means <- colMeans(healthy_abundances, na.rm = TRUE)
cancer_means  <- colMeans(cancer_abundances, na.rm = TRUE)

df_healthy <- data.frame(Species = names(healthy_means), Abundance = healthy_means, Group = "Healthy")
df_cancer  <- data.frame(Species = names(cancer_means), Abundance = cancer_means, Group = "Cancer")

# Step 9: Plot bar graphs
p_healthy <- ggplot(df_healthy, aes(x = reorder(Species, Abundance), y = Abundance)) +
  geom_bar(stat = "identity", fill = "#1b9e77") +
  coord_flip() +
  labs(title = "Top 30 Abundant Species in Healthy Group", x = "Species", y = "Mean Abundance") +
  theme_minimal(base_size = 12)

p_cancer <- ggplot(df_cancer, aes(x = reorder(Species, Abundance), y = Abundance)) +
  geom_bar(stat = "identity", fill = "#d95f02") +
  coord_flip() +
  labs(title = "Top 30 Abundant Species in Cancer Group", x = "Species", y = "Mean Abundance") +
  theme_minimal(base_size = 12)

# Step 10: Display plots
print(p_healthy)
print(p_cancer)

# Optional: Save to file
ggsave("top30_abundant_species_healthy.png", p_healthy, width = 10, height = 8, dpi = 300)
ggsave("top30_abundant_species_cancer.png", p_cancer, width = 10, height = 8, dpi = 300)

#Comparative analysis for Phylum
# Load necessary libraries
library(tidyverse)
library(ggpubr)
library(viridis)
#Read data
# Read Healthy dataset (HH)
hh_phylum <- read.csv("https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/HH_combined_bracken_phylum_fraction.csv",
                      row.names = 1, check.names = FALSE)

# Read Cancer dataset (PC)
pC_phylum <- read.csv("https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/PC_combined_bracken_phylum_fraction.csv",
                      row.names = 1, check.names = FALSE)

# Step 1: Transpose data (assumes rows are phyla, columns are samples)
hh_t <- as.data.frame(t(hh_phylum))
pC_t <- as.data.frame(t(pC_phylum))

# Step 2: Add group labels
hh_t$Group <- "Healthy"
pC_t$Group <- "Cancer"

# Step 3: Combine both datasets
combined <- bind_rows(hh_t, pC_t)
combined$Sample <- paste0("Sample_", seq_len(nrow(combined)))

# Step 4: Convert to long format for ggplot
long_data <- combined %>%
  pivot_longer(cols = -c(Group, Sample), names_to = "Phylum", values_to = "Abundance")

# Step 5a: Stacked Bar Plot (Per Sample)
ggplot(long_data, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group, scales = "free_x") +
  theme_minimal(base_size = 12) +
  labs(title = "Gut Microbiome Composition at Phylum Level",
       y = "Relative Abundance", x = "Sample") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_viridis_d(option = "D")

# Step 5b: Boxplot per Phylum with statistical comparison
ggboxplot(long_data, x = "Group", y = "Abundance", fill = "Group",
          palette = c("Healthy" = "#56B4E9", "Cancer" = "#E69F00"),
          facet.by = "Phylum", scales = "free_y") +
  stat_compare_means(method = "wilcox.test") +
  labs(title = "Phylum-Level Abundance Comparison Between Groups",
       y = "Relative Abundance")

#Check Character values
str(long_data$Abundance)




#Code this Again

# Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(pheatmap)

# Step 1: Read data (phylum in rows, samples in columns)
hh_url <- "https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/HH_combined_bracken_phylum_fraction.csv"
pc_url <- "https://raw.githubusercontent.com/okpecallistus/Uchenna5202stuff/refs/heads/main/PC_combined_bracken_phylum_fraction.csv"
View(hh_url)
hh_phylum <- read.csv(hh_url, row.names = 1, check.names = FALSE)
pC_phylum <- read.csv(pc_url, row.names = 1, check.names = FALSE)
View(hh_phylum)
View(pC_phylum)
# Step 2: Transpose and reshape to long format
hh_long <- hh_phylum %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.), Group = "Healthy") %>%
  pivot_longer(-c(Sample, Group), names_to = "Phylum", values_to = "Abundance")

pc_long <- pC_phylum %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.), Group = "Cancer") %>%
  pivot_longer(-c(Sample, Group), names_to = "Phylum", values_to = "Abundance")

# Step 3: Combine datasets
long_data <- bind_rows(hh_long, pc_long)

# -------------------------------
# ✅ Option 1: Stacked Bar Plot (Normalized)
long_data_norm <- long_data %>%
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

bar_plot<-ggplot(long_data_norm, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group, scales = "free_x") +
  theme_minimal(base_size = 12) +
  labs(title = "Gut Microbiome Composition at Phylum Level",
       y = "Relative Abundance", x = "Sample") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_viridis_d(option = "D")
print(bar_plot)
ggsave("phylum_barplot.png", bar_plot, width = 10, height = 6, dpi = 300)

# -------------------------------
# ✅ Option 2: Heatmap of Phylum Abundance

# Reshape to phylum x sample matrix
# Reshape to phylum x sample matrix: first summarize to avoid duplicates
heatmap_matrix <- long_data %>%
  group_by(Phylum, Sample) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("Phylum") %>%
  as.matrix()

# Optional: scale by phylum (row-wise)
scaled_matrix <- t(scale(t(heatmap_matrix)))

# Heatmap
png("phylum_heatmap.png", width = 1000, height = 800, res = 150) #Save heatmap
pheatmap(scaled_matrix,
         color = viridis(100),
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize = 10, main = "Phylum Abundance Heatmap")
png("phylum_heatmap.png", width = 1000, height = 800, res = 150) #Save heatmap
dev.off()


library(grid)
library(pheatmap)

p <- pheatmap(scaled_matrix,
              color = viridis(100),
              cluster_rows = TRUE, cluster_cols = TRUE,
              fontsize = 10, main = "Phylum Abundance Heatmap")

grid::grid.newpage()
grid::grid.draw(p$gtable)
#HEat MAp code

library(pheatmap)
library(grid)

# Generate heatmap and store it as a grob
p <- pheatmap(scaled_matrix,
              color = viridis(100),
              cluster_rows = TRUE, cluster_cols = TRUE,
              fontsize = 10, main = "Phylum Abundance Heatmap")

# Display in RStudio
grid::grid.newpage()
grid::grid.draw(p$gtable)

# Save to file
png("phylum_heatmap.png", width = 1000, height = 800, res = 150)
grid::grid.draw(p$gtable)
dev.off()


# 1. Show in RStudio plot viewer
pheatmap(scaled_matrix,
         color = viridis(100),
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize = 10, main = "Phylum Abundance Heatmap")

# 2. Save to file
png("phylum_heatmap.png", width = 1000, height = 800, res = 150)
pheatmap(scaled_matrix,
         color = viridis(100),
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize = 10, main = "Phylum Abundance Heatmap")
dev.off()
# Save the heatmap correctly
png("phylum_heatmap.png", width = 1000, height = 800, res = 150)

# Redraw the heatmap inside the png device
pheatmap(scaled_matrix,
         color = viridis(100),
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize = 10, main = "Phylum Abundance Heatmap")

dev.off()
getwd()
png("C:/Users/YourName/Documents/phylum_heatmap.png", width = 1000, height = 800, res = 150)



# -------------------------------
# ✅ Option 3: Group Mean Bar Plot

group_means <- long_data %>%
  group_by(Group, Phylum) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

ggplot(group_means, aes(x = reorder(Phylum, -MeanAbundance), y = MeanAbundance, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_viridis_d(option = "C") +
  labs(title = "Mean Relative Abundance of Phyla by Group",
       y = "Mean Relative Abundance", x = "Phylum") +
  theme_minimal(base_size = 12)

print(group_means)
ggsave("phylum_group_means.png", group_means, width = 8, height = 6, dpi = 300)

library(grid)

# Draw on screen
p <- pheatmap(scaled_matrix,
              color = viridis(100),
              cluster_rows = TRUE, cluster_cols = TRUE,
              fontsize = 10, main = "Phylum Abundance Heatmap")

# Save to file
png("phylum_heatmap.png", width = 1000, height = 800, res = 150)
grid.draw(p$gtable)
dev.off()







