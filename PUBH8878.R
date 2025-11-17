#### packages
library(vcfR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)



#function to keep biallelic sites
keep_biallelic_snps <- function(vcf) {
  bi <- is.biallelic(vcf) #vector marking biallelic sites
  vcf <- vcf[bi, ]         #subset
  is_snp <- nchar(vcf@fix[, "REF"]) == 1 & nchar(vcf@fix[, "ALT"]) == 1
  vcf[is_snp, ]
}

#make sure every variant has ID 
make_ids <- function(vcf) {
  ids <- vcf@fix[, "ID"]
  #if ID is missing, make one using position and alleles
  noid <- which(is.na(ids) | ids == ".")
  if (length(noid)) {
    ids[noid] <- paste0(vcf@fix[noid, "CHROM"], ":", vcf@fix[noid, "POS"],
                        "_", vcf@fix[noid, "REF"], ">", vcf@fix[noid, "ALT"])
  }
  vcf@fix[, "ID"] <- ids
  return(vcf)
}

#extract allele frequency values from the INFO field
extract_af <- function(vcf) {
  info <- vcf@fix[, "INFO"]
  data.frame(
    ID  = vcf@fix[, "ID"], #variant ID
    AF  = sub(".*AF=([^;]+).*", "\\1", info), #global AF
    AFR = sub(".*AFR_AF=([^;]+).*", "\\1", info),            
    AMR = sub(".*AMR_AF=([^;]+).*", "\\1", info),            
    EAS = sub(".*EAS_AF=([^;]+).*", "\\1", info),           
    EUR = sub(".*EUR_AF=([^;]+).*", "\\1", info),            
    SAS = sub(".*SAS_AF=([^;]+).*", "\\1", info),           
    stringsAsFactors = FALSE
  )
}

#gene regions from bash script
brca2 <- read.vcfR("brca2_region.vcf.gz")
chek2 <- read.vcfR("chek2_region.vcf.gz")
tp53  <- read.vcfR("tp53_region.vcf.gz")

#filter to biallelic SNPs and ensure IDs exist
brca2 <- make_ids(keep_biallelic_snps(brca2))
chek2 <- make_ids(keep_biallelic_snps(chek2))
tp53  <- make_ids(keep_biallelic_snps(tp53))

#extract AF information
brca2_af <- extract_af(brca2)
chek2_af <- extract_af(chek2)
tp53_af  <- extract_af(tp53)

#ddd gene labels
brca2_af$Gene <- "BRCA2"
chek2_af$Gene <- "CHEK2"
tp53_af$Gene  <- "TP53"

#combine all genes into one data frame
all_af <- rbind(brca2_af, chek2_af, tp53_af)

#reshape for plots
#wide to long format
af_long <- all_af %>%
  rename(Global_AF = AF) %>%  # rename for clarity
  pivot_longer(
    cols = c(AFR, AMR, EAS, EUR, SAS),
    names_to = "Population",
    values_to = "Pop_AF"
  ) %>%
  #convert to numeric and drop missing values
  mutate(
    Global_AF = as.numeric(ifelse(Global_AF %in% c(".", "", "NA"), NA, Global_AF)),
    Pop_AF = as.numeric(ifelse(Pop_AF %in% c(".", "", "NA"), NA, Pop_AF))
  ) %>%
  filter(!is.na(Pop_AF), Pop_AF >= 0, Pop_AF <= 1)

#summary table
summary_table <- af_long %>%
  group_by(Gene, Population) %>%
  summarise(
    Mean_AF = mean(Pop_AF), #average AF
    SD_AF = sd(Pop_AF), #variability in AF
    n = n(), #number of variants considered
    .groups = "drop"
  )
print(summary_table)

View(summary_table)
#plots
ggplot(af_long, aes(x = Population, y = Pop_AF, fill = Population)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8, width = 0.5) +
  geom_point(
    data = summary_table,
    aes(x = Population, y = Mean_AF),
    color = "black", shape = 21, size = 2.5, fill = "yellow"
  ) +
  geom_text(
    data = summary_table,
    aes(x = Population, y = Mean_AF, label = sprintf("%.4f", Mean_AF)),
    vjust = -0.3, size = 3, color = "black", fontface = "bold"
  ) +
  coord_cartesian(ylim = c(0, 0.015)) +  #zoom in to focus on variants
  facet_wrap(~ Gene, ncol = 1, scales = "fixed") +
  labs(
    title = "Allele Frequency Distributions by Population with Mean Values",
    y = "Allele Frequency (AF)",
    x = "Superpopulation"
  )


###################### F_st Analysis

#constant to prevent division by zero
eps <- 1e-9

#function to calculate mean pairwise fst between all population pairs for one gene
make_fst_table <- function(df_gene) {
  pops <- c("AFR", "AMR", "EAS", "EUR", "SAS") #define populations
  combs <- combn(pops, 2, simplify = FALSE) #create all pairwise combinations
  results <- purrr::map_dfr(combs, function(pair) { #loop
    p1 <- as.numeric(df_gene[[pair[1]]])  #AF in population 1
    p2 <- as.numeric(df_gene[[pair[2]]])  #AF frequencies in population 2
    p1 <- pmax(pmin(p1, 1 - eps), eps) #values between 0 and 1
    p2 <- pmax(pmin(p2, 1 - eps), eps)
    #Nei fst
    fst_values <- ((p1 - p2)^2) / (p1 * (1 - p1) + p2 * (1 - p2) + eps)
    #mean fst for this population pair
    tibble(
      Pop1 = pair[1],
      Pop2 = pair[2],
      FST = mean(fst_values, na.rm = TRUE)
    )
  })
  
  return(results)
}
fst_brca2 <- make_fst_table(brca2_af) %>% mutate(Gene = "BRCA2")
fst_chek2 <- make_fst_table(chek2_af) %>% mutate(Gene = "CHEK2")
fst_tp53  <- make_fst_table(tp53_af)  %>% mutate(Gene = "TP53")

#combine into df
fst_results <- bind_rows(fst_brca2, fst_chek2, fst_tp53)
print(fst_results)

#heat map
ggplot(fst_results, aes(x = Pop1, y = Pop2, fill = FST)) +
  geom_tile(color = "white") +  # draw squares for each pair
  scale_fill_gradient(low = "pink", high = "darkred") +  #color by FST
  facet_wrap(~ Gene, ncol = 1) +  #one panel p gene
  labs(
    title = "Pairwise FST Heatmap Across Populations",
    x = "Population",
    y = "Population",
    fill = "FST"
  ) 

############ PCA
#changing to long format
to_long <- function(df, gene_name) {
  df %>%
    dplyr::rename(Global_AF = AF) %>%  #rename first
    dplyr::select(ID, AFR, AMR, EAS, EUR, SAS, Global_AF) %>% 
    pivot_longer( #converts to long
      cols = c(AFR, AMR, EAS, EUR, SAS), #new column for pops
      names_to = "Population", #new column for af
      values_to = "AF"
    ) %>%
    mutate(Gene = gene_name) #gene labels
}


al <- bind_rows( #add afs
  to_long(brca2_af, "BRCA2"),
  to_long(chek2_af, "CHEK2"),
  to_long(tp53_af,  "TP53")
) %>%
  filter(!is.na(AF), AF >= 0, AF <= 1)

#pca prep
af_mat <- al %>%
  mutate(
    Var = paste(Gene, ID, sep = "::"),
    AF = as.numeric(AF)   #ensure AF is numeric
  ) %>%
  dplyr::select(Population, Var, AF) %>% #only pop, ID, AF cols
  tidyr::pivot_wider(
    names_from = Var,
    values_from = AF,
    values_fill = list(AF = 0) 
  )

#matrix
pop_labels <- af_mat$Population
X <- as.matrix(af_mat[,-1, drop=FALSE])

X <- af_mat %>%
  tibble::column_to_rownames("Population") %>% #pop as row names
  dplyr::select(where(~ var(.x, na.rm = TRUE) > 0)) #drop all 0 cols

#pca
pca <- prcomp(X, scale. = TRUE)
summary(pca)

pca <- prcomp(X, scale. = TRUE)
pdat <- data.frame(pca$x[,1:2], Population = pop_labels)

#plots
ggplot(pdat, aes(PC1, PC2, color = Population, label = Population)) +
  geom_point(size = 3) +
  geom_text_repel(show.legend = FALSE) +
  labs(
    title = "PCA of Population Mean Allele Frequencies (BRCA2, CHEK2, TP53)",
    #variance explained
    x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)")
  )
