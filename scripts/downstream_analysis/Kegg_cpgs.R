# "KEGG annotation"

  
#After building the age prediction model and identifying the top 20 CpG sites that contribute
#most strongly to the prediction, we performed functional annotation using the KEGG pathway database. 
#CpGs were classified as hypermethylated or hypomethylated based on their coefficients, 
#and enrichment analysis was conducted to explore the biological pathways they are involved in.



#Load require libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(kableExtra)
library(knitr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(readr)
library(patchwork)
library(kableExtra)

#Read csv file with top20cpgs 
top_cpgs <- read.csv("Top20_CpGs_ElasticNet.csv", 
                     stringsAsFactors = FALSE)


## Separate Cpgs: hyper- and hypomethylated 


hyper_cpgs <- top_cpgs$CpG[top_cpgs$Methylation_Change == "Hypermethylated"]
hypo_cpgs  <- top_cpgs$CpG[top_cpgs$Methylation_Change == "Hypomethylated"]


## CpG Annotation

# Load annotation
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k_df <- as.data.frame(ann450k)

# Merge top CpGs with gene names
annotated_cpgs <- top_cpgs %>%
  left_join(ann450k_df[, c("Name", "UCSC_RefGene_Name")],
            by = c("CpG" = "Name"))

# Split multiple genes per CpG
df_annotation <- annotated_cpgs %>%
  tidyr::separate_rows(UCSC_RefGene_Name, sep = ";") %>%
  dplyr::rename(Gene = UCSC_RefGene_Name) %>%
  dplyr::filter(Gene != "") %>%
  dplyr::select(Methylation_Change, Gene)




#Gene Conversion for KEGG Analysis

# Hypermethylated 
hyper_genes <- unique(df_annotation$Gene[df_annotation$Methylation_Change == "Hypermethylated"])
hyper_ids <- bitr(hyper_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Hypomethylated 
hypo_genes <- unique(df_annotation$Gene[df_annotation$Methylation_Change == "Hypomethylated"])
hypo_ids <- bitr(hypo_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


#KEGG Enrichment Analysis by Methylation Status


# KEGG for hypermethylated
hyper_kegg <- enrichKEGG(gene = hyper_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)

# KEGG for hypomethylated
hypo_kegg <- enrichKEGG(gene = hypo_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)



#Data visualization

#Cnetplots



cnetplot(hyper_kegg, showCategory = min(5, nrow(hyper_kegg)))
cnetplot(hypo_kegg, showCategory = min(5, nrow(hypo_kegg)))




#Dotplots 

prepare_top_kegg <- function(enrich_obj, top_n = 10) {
  enrich_obj@result %>%
    mutate(
      ID = as.character(ID),
      Description = as.character(Description),
      Count = as.numeric(Count),
      p.adjust = as.numeric(p.adjust)
    ) %>%
    arrange(p.adjust) %>%
    slice(1:top_n)
}

# -----------------------------
# Prepare top pathways
# -----------------------------
top_hyper <- prepare_top_kegg(hyper_kegg, top_n = 10)
top_hypo  <- prepare_top_kegg(hypo_kegg, top_n = 10)

# -----------------------------
# Create dotplot for hypermethylated CpGs
# -----------------------------
p_hyper <- ggplot(top_hyper, aes(x = reorder(Description, Count), y = Count)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(4, 10)) +  # Larger points
  labs(
    title = "KEGG Enriched Pathways (Hypermethylated CpGs)",
    x = "", 
    y = "Gene Count", 
    color = "-log10(adj p)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10)
  )

# -----------------------------
# Create dotplot for hypomethylated CpGs
# -----------------------------
p_hypo <- ggplot(top_hypo, aes(x = reorder(Description, Count), y = Count)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(4, 10)) +
  labs(
    title = "KEGG Enriched Pathways (Hypomethylated CpGs)",
    x = "", 
    y = "Gene Count", 
    color = "-log10(adj p)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10)
  )



ggsave("KEGG_Hypermethylated.png", plot = p_hyper, width = 10, height = 8, dpi = 300)
ggsave("KEGG_Hypomethylated.png", plot = p_hypo, width = 10, height = 8, dpi = 300)

ggplot(top_hyper, aes(x = reorder(Description, Count), y = Count)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(4, 10)) +  # Larger points
  labs(
    title = "KEGG Enriched Pathways (Hypermethylated CpGs)",
    x = "", 
    y = "Gene Count", 
    color = "-log10(adj p)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10)
  )


ggplot(top_hypo, aes(x = reorder(Description, Count), y = Count)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  coord_flip() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(4, 10)) +
  labs(
    title = "KEGG Enriched Pathways (Hypomethylated CpGs)",
    x = "", 
    y = "Gene Count", 
    color = "-log10(adj p)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10)
  )





#Barplots


barplot(hyper_kegg, showCategory = 10, title = "Hypermethylated CpG KEGG")
barplot(hypo_kegg, showCategory = 10, title = "Hypomethylated CpG KEGG")




#KEGG Pathway Enrichment by Methylation Status



# Hyper-methylated
hyper_table <- hyper_kegg@result
hyper_table$Methylation_Status <- "Hypermethylated"

# Hypo-methylated
hypo_table <- hypo_kegg@result
hypo_table$Methylation_Status <- "Hypomethylated"

# Combine tables
kegg_results <- rbind(
  hyper_table[, c("Description", "pvalue", "p.adjust", "Methylation_Status")],
  hypo_table[, c("Description", "pvalue", "p.adjust", "Methylation_Status")]
)

# Rename columns
colnames(kegg_results) <- c("Pathway_Name", "PValue", "Adjusted_PValue", "Methylation_Status")

# Add hyperlinks to KEGG pathways using the IDs (from the original object)
kegg_ids <- c(hyper_kegg@result$ID, hypo_kegg@result$ID)
kegg_results$Pathway_Name <- paste0(
  "<a href='https://www.kegg.jp/dbget-bin/www_bget?", kegg_ids, "' target='_blank'>",
  kegg_results$Pathway_Name,
  "</a>"
)


# Order by Adjusted_PValue
kegg_results <- kegg_results[order(kegg_results$Adjusted_PValue), ]

# Display table 
kegg_results %>%
  knitr::kable("html", escape = FALSE,
               caption = "<b style='text-align:center;'>KEGG Pathways Enriched by Methylation Status</b>") %>%
  kableExtra::kable_styling(full_width = FALSE,
                            bootstrap_options = c("striped", "hover"),
                            position = "center")


#Save results of KEGG analysis in CSV files

# Hypermethylated
hyper_table <- hyper_kegg@result
hyper_table$Methylation_Status <- "Hypermethylated"

# Hypomethylated
hypo_table <- hypo_kegg@result
hypo_table$Methylation_Status <- "Hypomethylated"

# Combine tables
kegg_results <- rbind(
  hyper_table[, c("Description", "ID", "pvalue", "p.adjust", "Methylation_Status")],
  hypo_table[, c("Description", "ID", "pvalue", "p.adjust", "Methylation_Status")]
)

# Rename columns
colnames(kegg_results) <- c("Pathway_Name", "KEGG_ID", "PValue", "Adjusted_PValue", "Methylation_Status")

# Create KEGG link column
kegg_results$KEGG_Link <- paste0("https://www.kegg.jp/dbget-bin/www_bget?", kegg_results$KEGG_ID)
kegg_results <- kegg_results[order(kegg_results$Adjusted_PValue), ]

# Save to CSV
write.csv(kegg_results, "KEGG_Pathways_Methylation.csv", row.names = FALSE)


#Gene names

final_cpg_genes <- annotated_cpgs %>%
  tidyr::separate_rows(UCSC_RefGene_Name, sep = ";") %>%
  dplyr::rename(Gene = UCSC_RefGene_Name) %>%
  dplyr::filter(Gene != "") %>%
  dplyr::select(CpG, Gene, Methylation_Change)


final_cpg_genes <- distinct(final_cpg_genes)

write.csv(final_cpg_genes,
          "Top20_CpGs_Annotated.csv",
          row.names = FALSE)



#Display table with genes
final_cpg_genes %>%
  head(20) %>%  
  knitr::kable(
    format = "html",
    caption = "<b style='text-align:center;'>Annotated Top 20 CpGs with Associated Genes</b>",
    align = "c",
    col.names = c("CpG Site", "Associated Gene", "Methylation Status")
  ) %>%
  kableExtra::kable_styling(
    full_width = FALSE,
    bootstrap_options = c("striped", "hover", "condensed"),
    position = "center",
    font_size = 13
  ) %>%
  kableExtra::column_spec(3, bold = TRUE, color = "white",
                          background = ifelse(final_cpg_genes$Methylation_Change == "Hypermethylated",
                                              "#E74C3C", "#2E86C1"))


