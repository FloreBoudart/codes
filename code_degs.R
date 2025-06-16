# Installer les packages nécessaires (si pas déjà faits)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "AnnotationHub", "org.At.tair.db"))
install.packages("readxl")
install.packages("ggplot2")
install.packages("openxlsx")

# Charger les bibliothèques
library(DESeq2)
library(clusterProfiler)
library(AnnotationHub)
library(readxl)
library(ggplot2)
library(openxlsx)

# Chemin où seront enregistrés les résultats générés par ce code
output_dir <- "/mnt/data-backup/boudart/results/r_analysis"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------
# 1. Lire les données
# ------------------------------

# === A. Importer le fichier CSV avec les comptages ===
counts_all <- read.csv("/mnt/data-backup/boudart/results/csv/counts_combined.csv", header = TRUE, row.names = 1)

# Garder uniquement les colonnes d'intérêt : WT1, WT2, WT3, OE7A, OE7B, OE7C
counts_wt_oe <- counts_all[, c("WT1_counts", "WT2_counts", "WT3_counts",
                               "OE7A_counts", "OE7B_counts", "OE7C_counts")]
#----
counts_wt_cris <- counts_all[, c("WT1_counts", "WT2_counts", "WT3_counts",
                                 "CRISPR6A_counts", "CRISPR6B_counts", "CRISPR6C_counts")]

# === B. Créer les métadonnées correspondantes ===
metadata_wt_oe <- data.frame(
  sample_wt_oe = colnames(counts_wt_oe),
  condition = factor(c(rep("Control", 3), rep("Treated", 3)), levels = c("Control", "Treated"))
)
rownames(metadata_wt_oe) <- metadata_wt_oe$sample_wt_oe

# Vérification des noms
colnames(counts_wt_oe) <- metadata_wt_oe$sample_wt_oe

# Convertir en matrice d’entiers
counts_wt_oe <- as.matrix(counts_wt_oe)
storage.mode(counts_wt_oe) <- "integer"

#----
metadata_wt_cris <- data.frame(
  sample_wt_cris = colnames(counts_wt_cris),
  condition = factor(c(rep("Control", 3), rep("Treated", 3)), levels = c("Control", "Treated"))
)
rownames(metadata_wt_cris) <- metadata_wt_cris$sample_wt_cris
colnames(counts_wt_cris) <- metadata_wt_cris$sample_wt_cris

# Convertir en matrice d’entiers
counts_wt_cris <- as.matrix(counts_wt_cris)
storage.mode(counts_wt_cris) <- "integer"


# ------------------------------
# 2. Analyse différentielle avec DESeq2
# ------------------------------

dds_wt_oe <- DESeqDataSetFromMatrix(countData = counts_wt_oe, colData = metadata_wt_oe, design = ~ condition)
dds_wt_oe <- DESeq(dds_wt_oe)
res_wt_oe <- results(dds_wt_oe, contrast = c("condition", "Treated", "Control"))

# Extraire les gènes significatifs
res_significant_wt_oe <- res_wt_oe[which(res_wt_oe$padj < 0.05), ]
degs_wt_oe <- rownames(res_significant_wt_oe)

#----
dds_wt_cris <- DESeqDataSetFromMatrix(countData = counts_wt_cris, colData = metadata_wt_cris, design = ~ condition)
dds_wt_cris <- DESeq(dds_wt_cris)
res_wt_cris <- results(dds_wt_cris, contrast = c("condition", "Treated", "Control"))

res_significant_wt_cris <- res_wt_cris[which(res_wt_cris$padj < 0.05), ]
degs_wt_cris <- rownames(res_significant_wt_cris)


# ------------------------------
# 3. Enrichissement GO
# ------------------------------

# Annotation pour Nicotiana tabacum
# L'identifiant AH117331 dans AnnotationHub correspond à une base de données d'annotation de type OrgDb pour Nicotiana tabacum
ah <- AnnotationHub()
orgNtab <- ah[["AH117331"]]

# Conversion des gènes SYMBOL → ENTREZID
degs_converted_wt_oe <- bitr(degs_wt_oe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgNtab)

# Enrichissement GO (Biological Process)
ego_wt_oe <- enrichGO(gene = degs_converted_wt_oe$ENTREZID,
                      OrgDb = orgNtab,
                      keyType = "ENTREZID",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

#-----
degs_converted_wt_cris <- bitr(degs_wt_cris, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgNtab)

ego_wt_cris <- enrichGO(gene = degs_converted_wt_cris$ENTREZID,
                        OrgDb = orgNtab,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)

# ------------------------------
# 4. Visualisation des résultats GO
# ------------------------------

# Plot
dotplot(ego_wt_oe, showCategory = 12, title = "GO Term Enrichment (OE vs WT)")
ggsave(file.path(output_dir, "GO_Term_Enrichment_OE_vs_WT_dotplot.pdf"), width = 10, height = 8)


barplot(ego_wt_oe, showCategory = 12, title = "GO Term Enrichment (OE vs WT)")
ggsave(file.path(output_dir, "GO_Term_Enrichment_OE_vs_WT_barplot.pdf"), width = 10, height = 8)

#---
dotplot(ego_wt_cris, showCategory = 12, title = "GO Term Enrichment (CRISPR vs WT)")
ggsave(file.path(output_dir, "GO_Term_Enrichment_CRISPR_vs_WT_dotplot.pdf"), width = 10, height = 8)

barplot(ego_wt_cris, showCategory = 12, title = "GO Term Enrichment (CRISPR vs WT)")
ggsave(file.path(output_dir, "GO_Term_Enrichment_CRISPR_vs_WT_barplot.pdf"), width = 10, height = 8)

# ------------------------------
# 5. Fusionner les fichiers des biotypes
# ------------------------------


# === WT ===
wt_file1 <- read_excel("/mnt/data-backup/boudart/results/excel/WT1_counts_counts_enriched_counts.xlsx")
wt_file2 <- read_excel("/mnt/data-backup/boudart/results/excel/WT2_counts_counts_enriched_counts.xlsx")
wt_file3 <- read_excel("/mnt/data-backup/boudart/results/excel/WT3_counts_counts_enriched_counts.xlsx")
merged_wt <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE),
                    list(wt_file1, wt_file2, wt_file3))

# === OE ===
oe_file1 <- read_excel("/mnt/data-backup/boudart/results/excel/OE7A_counts_counts_enriched_counts.xlsx")
oe_file2 <- read_excel("/mnt/data-backup/boudart/results/excel/OE7B_counts_counts_enriched_counts.xlsx")
oe_file3 <- read_excel("/mnt/data-backup/boudart/results/excel/OE7C_counts_counts_enriched_counts.xlsx")
merged_oe <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE),
                    list(oe_file1, oe_file2, oe_file3))

# === CRISPR ===
crispr1 <- read_excel("/mnt/data-backup/boudart/results/excel/CRISPR6A_counts_counts_enriched_counts.xlsx")
crispr2 <- read_excel("/mnt/data-backup/boudart/results/excel/CRISPR6B_counts_counts_enriched_counts.xlsx")
crispr3 <- read_excel("/mnt/data-backup/boudart/results/excel/CRISPR6C_counts_counts_enriched_counts.xlsx")
merged_crispr <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE),
                        list(crispr1, crispr2, crispr3))

# Fusion WT + OE et WT + Crispr
gene_info_wt_oe <- merge(merged_wt, merged_oe, by = "Geneid", all = TRUE)

gene_info_wt_cris <- merge(merged_wt, merged_crispr, by = "Geneid", all = TRUE)

#--------------------
#  6. Répartition des Biotypes parmi les Gènes Différentiels dans WT - OE 
#--------------------

# Sélectionner les colonnes Geneid + toutes les colonnes de biotype spécifiques
biotype_all_wt_oe <- gene_info_wt_oe %>%
  dplyr::select(Geneid,
                gene_biotype.x.x,
                gene_biotype.y.x,
                gene_biotype.x,
                gene_biotype.x.y,
                gene_biotype.y.y,
                gene_biotype.y)

# Obtenir le premier biotype non-NA pour chaque ligne
biotype_first_wt_oe <- biotype_all_wt_oe %>%
  rowwise() %>%
  mutate(Gene_Biotype = coalesce(
    gene_biotype.x.x,
    gene_biotype.y.x,
    gene_biotype.x,
    gene_biotype.x.y,
    gene_biotype.y.y,
    gene_biotype.y)) %>%
  ungroup() %>%
  dplyr::select(Geneid, Gene_Biotype)

# Filtrer les gènes différentiels (degs)
biotype_degs_wt_oe <- biotype_first_wt_oe %>%
  filter(Geneid %in% degs_wt_oe) %>%
  filter(!is.na(Gene_Biotype))

# Compter le nombre de DEGs par biotype + pourcentage
biotype_counts_degs_wt_oe <- biotype_degs_wt_oe %>%
  group_by(Gene_Biotype) %>%
  summarise(Nombre = n()) %>%
  ungroup() %>%
  mutate(Pourcentage = Nombre / sum(Nombre) * 100,
         Label = paste0(Gene_Biotype, " (", round(Pourcentage, 1), "%)"))

# Graphique final
ggplot(biotype_counts_degs_wt_oe, aes(x = reorder(Gene_Biotype, -Nombre), 
                                      y = Nombre, 
                                      fill = Gene_Biotype)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Répartition des Biotypes parmi les Gènes Différentiels (WT vs OE)",
       x = "Biotype de Gène",
       y = "Nombre de Gènes",
       fill = "Biotype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = scales::hue_pal()(nrow(biotype_counts_degs_wt_oe)),
                    labels = biotype_counts_degs_wt_oe$Label)


ggsave(file.path(output_dir, "Biotypes_DEGs_WT_vs_OE.pdf"), width = 10, height = 8)

#--------------------
# 6 bis. Répartition des Biotypes parmi les Gènes Différentiels dans WT - Crispr
#--------------------

# Sélectionner les colonnes Geneid + toutes les colonnes de biotype spécifiques
biotype_all_wt_cris <- gene_info_wt_cris %>%
  dplyr::select(Geneid,
                gene_biotype.x.x,
                gene_biotype.y.x,
                gene_biotype.x,
                gene_biotype.x.y,
                gene_biotype.y.y,
                gene_biotype.y)

# Obtenir le premier biotype non-NA pour chaque ligne
biotype_first_wt_cris <- biotype_all_wt_cris %>%
  rowwise() %>%
  mutate(Gene_Biotype = coalesce(
    gene_biotype.x.x,
    gene_biotype.y.x,
    gene_biotype.x,
    gene_biotype.x.y,
    gene_biotype.y.y,
    gene_biotype.y)) %>%
  ungroup() %>%
  dplyr::select(Geneid, Gene_Biotype)

biotype_degs_wt_cris <- biotype_first_wt_cris %>%
  filter(Geneid %in% degs_wt_cris) %>%
  filter(!is.na(Gene_Biotype))

biotype_counts_degs_wt_cris <- biotype_degs_wt_cris %>%
  group_by(Gene_Biotype) %>%
  summarise(Nombre = n()) %>%
  ungroup() %>%
  mutate(Pourcentage = Nombre / sum(Nombre) * 100,
         Label = paste0(Gene_Biotype, " (", round(Pourcentage, 1), "%)"))

ggplot(biotype_counts_degs_wt_cris, aes(x = reorder(Gene_Biotype, -Nombre), 
                                        y = Nombre, 
                                        fill = Gene_Biotype)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Répartition des Biotypes parmi les Gènes Différentiels (WT vs CRISPR)",
       x = "Biotype de Gène",
       y = "Nombre de Gènes",
       fill = "Biotype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = scales::hue_pal()(nrow(biotype_counts_degs_wt_cris)),
                    labels = biotype_counts_degs_wt_cris$Label)

ggsave(file.path(output_dir, "Biotypes_DEGs_WT_vs_CRISPR.pdf"), width = 10, height = 8)

# ------------------------------
# 7. Créer un tableau final avec DESeq2 et annotations
# ------------------------------

#--------------------
# WT - OE
#--------------------

res_df_wt_oe <- as.data.frame(res_wt_oe)
res_df_wt_oe$ENSEMBL <- rownames(res_df_wt_oe)

# Extraire les colonnes utiles : Geneid (=ENSEMBL), descriptions, chromosome
gene_info_subset_wt_oe <- gene_info_wt_oe[, c("Geneid", 
                                              "description.x.x", "description.y.x", 
                                              "description.x", "description.x.y", 
                                              "description.y.y", "description.y", 
                                              "Chr.x.x", "Chr.x", "Chr.x.y", "Chr.y.y", "Chr.y")]

# Créer une colonne "Function" en prenant la première description non vide
gene_info_subset_wt_oe$Function <- apply(gene_info_subset_wt_oe[, 2:7], 1, function(x) {
  val <- x[!is.na(x) & x != ""]
  if (length(val) > 0) return(val[1]) else return(NA)})

# Créer une colonne "Chromosome" en prenant le premier champ de chromosome non vide
gene_info_subset_wt_oe$Chromosome <- apply(gene_info_subset_wt_oe[, 8:12], 1, function(x) {
  val <- x[!is.na(x) & x != ""]
  if (length(val) > 0) return(val[1]) else return(NA)})

# Garder uniquement les colonnes utiles
gene_info_final_wt_oe <- gene_info_subset_wt_oe[, c("Geneid", "Function", "Chromosome")]


# Fusionner DESeq2 + annotations
final_table_wt_oe <- merge(res_df_wt_oe, gene_info_final_wt_oe, by.x = "ENSEMBL", by.y = "Geneid", all.x = TRUE)

# Réorganiser les colonnes dans l'ordre souhaité
final_table_wt_oe <- final_table_wt_oe[, c("ENSEMBL", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "Function", "Chromosome")]

# Créer un fichier Excel
write.xlsx(final_table_wt_oe, file = file.path(output_dir, "DEGs_OE_vs_WT_full_table.xlsx"), rowNames = FALSE)


#--------------------
# WT - Crispr
#--------------------

# === Créer tableau final WT vs CRISPR ===
res_df_wt_cris <- as.data.frame(res_wt_cris)
res_df_wt_cris$ENSEMBL <- rownames(res_df_wt_cris)

# Extraire les colonnes utiles pour annotation
gene_info_subset_wt_cris <- gene_info_wt_cris[, c("Geneid", 
                                                  "description.x.x", "description.y.x", 
                                                  "description.x", "description.x.y", 
                                                  "description.y.y", "description.y", 
                                                  "Chr.x.x", "Chr.x", "Chr.x.y", "Chr.y.y", "Chr.y")]

# Fonction pour récupérer première annotation valide
gene_info_subset_wt_cris$Function <- apply(gene_info_subset_wt_cris[, 2:7], 1, function(x) {
  val <- x[!is.na(x) & x != ""]
  if (length(val) > 0) return(val[1]) else return(NA)})

gene_info_subset_wt_cris$Chromosome <- apply(gene_info_subset_wt_cris[, 8:12], 1, function(x) {
  val <- x[!is.na(x) & x != ""]
  if (length(val) > 0) return(val[1]) else return(NA)})

# Nettoyage
gene_info_final_wt_cris <- gene_info_subset_wt_cris[, c("Geneid", "Function", "Chromosome")]

# Fusion avec résultats DESeq2
final_table_wt_cris <- merge(res_df_wt_cris, gene_info_final_wt_cris, by.x = "ENSEMBL", by.y = "Geneid", all.x = TRUE)

# Réorganisation des colonnes
final_table_wt_cris <- final_table_wt_cris[, c("ENSEMBL", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "Function", "Chromosome")]

# Export en Excel
write.xlsx(final_table_wt_oe, file = file.path(output_dir, "DEGs_CRISPR_vs_WT_full_table.xlsx"), rowNames = FALSE)


# ------------------------------
# 8. Visualisation de tous les biotypes + pourcentage (WT vs OE)
# ------------------------------


# Nettoyer et compter
biotype_all_genes_wt_oe <- biotype_first_wt_oe %>%
  filter(!is.na(Gene_Biotype))

biotype_counts_all_wt_oe <- biotype_all_genes_wt_oe %>%
  group_by(Gene_Biotype) %>%
  summarise(Nombre = n()) %>%
  ungroup() %>%
  mutate(Pourcentage = Nombre / sum(Nombre) * 100,
         Label = paste0(Gene_Biotype, " (", round(Pourcentage, 1), "%)"))

# Préparer facteur pour avoir l’ordre correct dans la légende
biotype_counts_all_wt_oe$Gene_Biotype <- factor(biotype_counts_all_wt_oe$Gene_Biotype,
                                                levels = biotype_counts_all_wt_oe$Gene_Biotype)
names_label_wt_oe <- setNames(biotype_counts_all_wt_oe$Label, biotype_counts_all_wt_oe$Gene_Biotype)

# Barplot : Axe Y = nombre de gènes, légende = % inclus
ggplot(biotype_counts_all_wt_oe, aes(x = reorder(Gene_Biotype, -Nombre), 
                                     y = Nombre, 
                                     fill = Gene_Biotype)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Répartition des Biotypes pour tous les Gènes (WT vs OE)",
       x = "Biotype de Gène",
       y = "Nombre de Gènes",
       fill = "Biotype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = scales::hue_pal()(nrow(biotype_counts_all_wt_oe)),
                    labels = names_label_wt_oe)

ggsave(file.path(output_dir, "Biotypes_All_Genes_WT_vs_OE.pdf"), width = 10, height = 8)


# ------------------------------
# 8 bis. Visualisation des biotypes pour tous les gènes (WT vs CRISPR)
# ------------------------------

# Nettoyer et compter
biotype_all_genes_wt_cris <- biotype_first_wt_cris %>%
  filter(!is.na(Gene_Biotype))

biotype_counts_all_wt_cris <- biotype_all_genes_wt_cris %>%
  group_by(Gene_Biotype) %>%
  summarise(Nombre = n()) %>%
  ungroup() %>%
  mutate(Pourcentage = Nombre / sum(Nombre) * 100,
         Label = paste0(Gene_Biotype, " (", round(Pourcentage, 1), "%)"))

# Préparer facteur pour avoir l’ordre correct dans la légende
biotype_counts_all_wt_cris$Gene_Biotype <- factor(biotype_counts_all_wt_cris$Gene_Biotype,
                                                  levels = biotype_counts_all_wt_cris$Gene_Biotype)
names_label_wt_cris <- setNames(biotype_counts_all_wt_cris$Label, biotype_counts_all_wt_cris$Gene_Biotype)

# Barplot : Axe Y = nombre de gènes, légende = % inclus
ggplot(biotype_counts_all_wt_cris, aes(x = reorder(Gene_Biotype, -Nombre), 
                                       y = Nombre, 
                                       fill = Gene_Biotype)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Répartition des Biotypes pour tous les Gènes (WT vs CRISPR)",
       x = "Biotype de Gène",
       y = "Nombre de Gènes",
       fill = "Biotype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = scales::hue_pal()(nrow(biotype_counts_all_wt_cris)),
                    labels = names_label_wt_cris)

ggsave(file.path(output_dir, "Biotypes_All_Genes_WT_vs_CRISPR.pdf"), width = 10, height = 8)

# ------------------------------
# 9. Visualisation des GO term en dotplot mais avec l'axe x = foldchange (WT vs OE)
# ------------------------------

# Convertir ego_wt_oe en dataframe
ego_wt_oe_df <- as.data.frame(ego_wt_oe)

# Vérifiez les noms des colonnes
colnames(ego_wt_oe_df)

# Vérifiez les valeurs manquantes dans la colonne p.adjust (p_adj)
sum(is.na(ego_wt_oe_df$p.adjust))

# Appliquez le filtre (si p.adjust est bien la colonne de p-value ajustées)
ego_wt_oe_filtre <- ego_wt_oe_df[ego_wt_oe_df$p.adjust < 0.05, ]


# Ajouter une colonne log2FoldChange à partir de FoldEnrichment
ego_wt_oe_filtre$log2FoldChange <- log2(ego_wt_oe_filtre$FoldEnrichment)


ggplot(ego_wt_oe_filtre, aes(x = log2FoldChange, y = Description, size = Count, color = p.adjust)) +
  geom_point() +  # Utilisation des points pour le dotplot
  theme_minimal() +  # Thème minimal pour la présentation
  scale_color_gradient(low = "blue", high = "red") +  # Gradient de couleurs pour p.adjust
  labs(
    title = "GO Term Enrichment WT vs OE",
    x = "log2(Fold Enrichment)",
    y = "GO Term Description",
    color = "Adjusted P-Value",
    size = "Gene Count"
  ) +
  theme(axis.text.y = element_text(size = 15)) + 
  theme(panel.border = element_rect(color = "grey", fill = NA, size = 1))


# Sauvegarder le graphique
ggsave(file.path(output_dir, "GO_Term_Enrichment_with_log2FoldChange_OE_vs_WT_dotplot.pdf"), width = 10, height = 8)


# ------------------------------
# 9 bis. Visualisation des GO term en dotplot mais avec l'axe x = foldchange  (WT vs CRISPR)
# ------------------------------

# Convertir ego_wt_oe en dataframe
ego_wt_cris_df <- as.data.frame(ego_wt_cris)

# Vérifiez les noms des colonnes
colnames(ego_wt_cris_df)

# Vérifiez les valeurs manquantes dans la colonne p.adjust (p_adj)
sum(is.na(ego_wt_cris_df$p.adjust))

# Appliquez le filtre (si p.adjust est bien la colonne de p-value ajustées)
ego_wt_cris_filtre <- ego_wt_cris_df[ego_wt_cris_df$p.adjust < 0.05, ]


# Ajouter une colonne log2FoldChange à partir de FoldEnrichment
ego_wt_cris_filtre$log2FoldChange <- log2(ego_wt_cris_filtre$FoldEnrichment)


ggplot(ego_wt_cris_filtre, aes(x = log2FoldChange, y = Description, size = Count, color = p.adjust)) +
  geom_point() +  # Utilisation des points pour le dotplot
  theme_minimal() +  # Thème minimal pour la présentation
  scale_color_gradient(low = "blue", high = "red") +  # Gradient de couleurs pour p.adjust
  labs(
    title = "GO Term Enrichment WT vs CRISPR",
    x = "log2(Fold Enrichment)",
    y = "GO Term Description",
    color = "Adjusted P-Value",
    size = "Gene Count"
  ) +
  theme(axis.text.y = element_text(size = 15)) + 
  theme(panel.border = element_rect(color = "grey", fill = NA, size = 1))

# Sauvegarder le graphique
ggsave(file.path(output_dir, "GO_Term_Enrichment_with_log2FoldChange_CRISPR_vs_WT_dotplot.pdf"), width = 10, height = 8)


