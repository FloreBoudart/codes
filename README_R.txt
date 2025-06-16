
Ce code R permet de faire un GO Term Enrichment sur nos données et génère bien d'autres résultats. 

Actuellement il est paramétré pour travailler sur 3 échantillons de WT, OE7 et CRISPR6. 
Le code est commenté ce qui aide a comprendre globalement sa structure. 

Il y a 9 points principaux qui permettent d'analyser les données et voici ce qu'ils font:  

1.	Lire les données 
		-> Récupère le fichier csv contenant les counts
		-> Créer les fichiers counts et métadata pour comparer wt-oe et wt-cripsr
2.	Analyse différentielle avec DESeq2
		-> Normalisation et extrait les gènes significatifs 
3.	Enrichissement GO
		-> Récupère la base de donnée correspondant a Nicotiana Tabacum 
		-> Créer les fichiers Enrichissement Go pour wt-oe et wt-cripsr
4.	Visualisation des résultats GO
		-> Sous forme de dotplot (y = terms go ; x =  generatio)
		-> Sous forme de barplot (y = terms go ; x = Count ) 
		-> Enregistrement des fichiers sous forme pdf
5.	Fusionner les fichiers des biotypes (venant des fichiers excel provenant du pipeline) 
		-> Fusionner les fichiers excel selon l'échantillion
		-> Fusionner les fichiers pour obtenir wt-oe et wt-crispr
6.	Répartition des Biotypes parmi les Gènes Différentiels
		-> Récuperer les biotypes, filter les gènes différentiels, les compter et 
			calculer le pourcentage 
		-> Visualisation sous forme de barplot avec le pourcentage des biotypes 
			wt-oe et wt-crispr (y = nombres de gènes ; x = terms go)
		-> Enregistre les fichiers sous forme pdf 
7.	Créer un tableau final avec DESeq2 et annotations
		-> Récuperer les degs (gènes différentiellement exprimés)les descriptions 
			des fichiers excel 
		-> Renommer les descriptions "function" 
		-> Enregistrer les tableaux excel, wt-oe et wt-crispr, qui contienent les colonnes : 
			ENSEMBL, baseMean,log2FoldChange, lfcSE, stat, pvalue, padj, Function, Chromosome
8.	Visualisation de tous les biotypes + pourcentage 
		-> Nettoyer, compter et faire les pourcentage de tous les biotypes 
		-> Barplot pour wt-oe et wt-crispr (y = Nombres de gènes; x = biotypes)
			+ pourcentage dans la légende 
		-> Enregistrement des fichiers sous forme pdf 
9.  Visualisation des GO term en dotplot mais avec l'axe x = foldchange (WT vs OE)
		-> Récuperer les fichiers Enrichissement Go
		-> Filter pour garder que les padj significatifs 
		-> Ajouter une colonne log2FoldChange à partir de FoldEnrichment
		-> Dotplot pour wt-oe et wt-crispr (y = terms go ; x =  log2FoldChange)
		-> Enregistrer les fichiers sous forme pdf 


Pour utiliser ce code :

1.	Changer le chemin aux lignes 17, 28, 144 à 160
2. 	Changer le noms des échantillons 
3.	Lancer le code 

Ce code R a été écrit dans le cadre d'un stage




-----------------------------------------------------------------------------

Si vous voulez changer les chemins manuels en automatique, il faudra changer :



1. Début de code à remplacer par : 

base_dir <- getwd()
output_dir <- file.path(base_dir, "results", "r_analysis")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


2. Ligne 28 remplacé par : 

counts_all <- read.csv(file.path(base_dir, "results", "csv", "counts_combined.csv"), header = TRUE, row.names = 1)



3. Ligne 124-125 : (+ suivre le même raisonnement pour les lignes 128 à 136)

dotplot_obj <- dotplot(ego_wt_oe, showCategory = 12, title = "GO Term Enrichment (OE vs WT)")
ggsave(filename = file.path(output_dir, "GO_Term_Enrichment_OE_vs_WT_dotplot.pdf"), plot = dotplot_obj, width = 10, height = 8)



4. Lignes 143 à 148 : (+ suivre le même raisonnement pour les lignes 151 à 162)

wt_file1 <- read_excel(file.path(base_dir, "results", "excel", "WT1_counts_counts_enriched_counts.xlsx"))
wt_file2 <- read_excel(file.path(base_dir, "results", "excel", "WT2_counts_counts_enriched_counts.xlsx"))
wt_file3 <- read_excel(file.path(base_dir, "results", "excel", "WT3_counts_counts_enriched_counts.xlsx"))
merged_wt <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE),
                    list(wt_file1, wt_file2, wt_file3))




5. Lignes 210 à 225 : (+ suivre le même raisonnement pour les lignes 265 à 279)

p <- ggplot(biotype_counts_degs_wt_oe, aes(x = reorder(Gene_Biotype, -Nombre), 
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

ggsave(filename = file.path(output_dir, "Biotypes_DEGs_WT_vs_OE.pdf"), plot = p, width = 10, height = 8)



6. Ligne 320 : (+ suivre le même raisonnement pour la ligne 357)

write.xlsx(final_table_wt_oe, file = file.path(output_dir, "DEGs_OE_vs_WT_full_table.xlsx"), rowNames = FALSE)



7. Lignes 382 à 396 : (+ suivre le même raisonnement pour les lignes 420 à 432)

p <- ggplot(biotype_counts_all_wt_oe, aes(x = reorder(Gene_Biotype, -Nombre), 
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

ggsave(filename = file.path(output_dir, "Biotypes_All_Genes_WT_vs_OE.pdf"), plot = p, width = 10, height = 8)



8. Lignes 457 à 473 : (+ suivre le même raisonnement pour les lignes 497 à 512)

p <- ggplot(ego_wt_oe_filtre, aes(x = log2FoldChange, y = Description, size = Count, color = p.adjust)) +
  geom_point() +
  theme_minimal() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "GO Term Enrichment WT vs OE",
    x = "log2(Fold Enrichment)",
    y = "GO Term Description",
    color = "Adjusted P-Value",
    size = "Gene Count"
  ) +
  theme(
    axis.text.y = element_text(size = 15),
    panel.border = element_rect(color = "grey", fill = NA, size = 1)
  )

# Sauvegarde automatique
ggsave(filename = file.path(output_dir, "GO_Term_Enrichment_with_log2FoldChange_OE_vs_WT_dotplot.pdf"),plot = p,width = 10,height = 8)
