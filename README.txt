
Pipeline Nextflow RNA-Seq


Ce pipeline sert à effectuer des analyses RNA-Seq sur le serveur LIBST de l'UCLouvain. Il génère
les résultats basiques pour pouvoir faire n'importe quelles statistiques par la suite. 

Il a été écrit en 3 mois par Boudart Flore, de début février 2025 à fin avril 2025 dans le cadre 
de stage et de TFE de dernière année en bachelier biotechnique option bioinformatique à la HEH 
de Mons. 

Il est fonctionnel sur les données qui possèdent un génome de référence et une annotation, 
quelques modifications seront surement à faire par la suite pour une utilisation plus facile et 
correcte/complète pour tous. 


Ce pipeline est écrit en Nextflow, qui est un gestionnaire de workflows. Celui-ci comporte les
étapes les plus importante pour effectuer une analyse en RNA-Seq. 

Il est composé de plusieurs programmes :

1.	FastQC -> Vérification de la qualité des données brutes
2.	Fastp -> Nettoyage des mauvaises séquences
3.	Trimmomatic -> Permet de couper les adaptateurs et modifier la tailles des séquences
4.	Hisat2 -> Fait l'indexation du génome de référence et l'alignement global des séquences 
		nettoyées sur ce dernier
5.	Samtools -> Utile pour l'analyse de données de séquençage à haut débit, il génère les 
		fichier bam, les trie, les index et en sort les statistiques d'alignement
6.	FeaturesCounts -> Compte le nombre de features (caractéristiques) qu'il y a dans les fichiers 
		sorte.bam en utilisant le fichier gtf
7.	CreateCSV (code python) -> Récupère le fichier généré par featuresCounts et le transforme en 
		tableau en ne gardant que les colonnes des identifiants de gènes et les counts. Cela sera
		utilisé pour les statistiques en R
8.	CompleteCounts (code python) -> Récupère le fichier généré par featuresCounts et y rajoute 
		les colonnes "description", "gbkey" et "gene_biotype" qui sont récupérées du fichier gtf. 
		Qui sont utile pour savoir quelle est le rôle de chaque gène
9.	ConvertToExcel(code python) -> Transforme le fichier de Complete_counts en version Excel 
		pour lire et se balader dans le fichier plus facilement
10.	MultiQC -> Regroupe tous les résultats du pipeline



Les résultats qui sont générés par chaque outils :  

1.	FastQC -> génère des fichiers html et des zip
2.	Fastp -> genère des clean.fastq.gz
3.	Trimmomatic -> genère des trim.fastq.gz
4.	Hisat2 -> Indexation au format .ht2 
		   -> Alignement au format .sam 
5.	Samtools -> Fichiers .sam sont transformé en .bam 
			 -> Fichiers .bam sont trié et ressortent en .sorted.bam
			 -> Fichiers .sorted.bam sont indexé en .bai
			 -> Les statistiques d'alignement ressortent dans les fichiers .flagstat
6.	FeaturesCounts -> Comptage au format counts.txt et les résumé au format .txt.summary
7.	CreateCSV (code python) -> Un seul fichier au format counts_combined.csv
8.	CompleteCounts (code python) -> des fichiers au format counts_counts_enriched.txt
9.	ConvertToExcel(code python) -> des fichiers au format counts_counts_enriched_counts.xlsx
10.	MultiQC -> Permet la visualisation en format html 
 


Pour l'utiliser, il faut des données brutes sorties de NGS au format fastq.gz. Ainsi que le 
fichier d'annotation et le génome de référence qui sont associés. Pour les déziper, effectues 
la commande "gunzip nom_du_fichier.gz"  

L'arborescence des fichiers doit être comme dans l'exemple ci-dessous. Tu peux vérifier avec la 
commande "tree -L 2", elle permet d'afficher l'arbre de fichier sur 2 étages/branches.

│   
└── analyse
    ├── conda
    │   └── env_conda.yaml
    ├── help
    │   ├── README.txt
    │   └── README_R.txt
    ├── genome
    │   ├── Genome_de_référence.fasta
    │   └── Annotation.gtf
    ├── script
    │   ├── create_csv.py
    │   ├── enrich_counts.py
    │   ├── convert_excel.py
    │   └── code_degs.R 
    ├── fastq
    │   ├── sample1_R1.fastq.gz
    │   ├── sample1_R2.fastq.gz
    │   ├── sample2_R1.fastq.gz
    │   ├── sample2_R2.fastq.gz
    │   ├── sample3_R1.fastq.gz
    │   └── sample3_R2.fastq.gz  
    ├── nextflow.config
    ├── run.sh
    ├── samplesheet.csv
    └── pipeline_rna.nf


Exemple du contenu de "samplesheet.csv" qui contient le nom de l'achantillon, et le chemin de 
chaque reads. 

sample,fastq_1,fastq_2
sample1,/chemin/vers/fastq/sample1_R1.fastq.gz,/chemin/vers/fastq/sample1_R2.fastq.gz
sample2,/chemin/vers/fastq/sample2_R1.fastq.gz,/chemin/vers/fastq/sample2_R2.fastq.gz
sample3,/chemin/vers/fastq/sample3_R1.fastq.gz,/chemin/vers/fastq/sample3_R2.fastq.gz



Exemple du contenu de "run.sh" permet de modifier les paramètres selon vos analyses. 

#!/usr/bin/env bash

nextflow run pipeline_rna.nf \
    --reads "samplesheet.csv" \
    --fasta "/chemin/vers/fichier.fasta" \
    --gtf "/chemin/vers/fichier.gtf" \
    --outdir "results" \
    -profile conda \
    -resume



Pour lancer une analyse, il faut : 
1.	Ouvrir l’application MobaxTerm
2.	Se conneceter sur le serveur avec son mot de passe
3.	Si cela n’est pas déjà fait, créer un environnement conda avec la commande conda create env -f env_conda.yaml
4.	Activer un environnement conda avec la commande conda activate env 
5.	Se placer au bon endroit dans le serveur, là où sont les données 
6.	Créer le fichier samplesheet.csv avec le nom des échantillons et les chemins des reads
7.	Créer le fichier run.sh en précisant les bons chemins des fichiers fasta et gtf ainsi que de préciser le nom du fichier de résultats.
8.	Effectuer la commande chmod+x run.sh pour activer ce fichier
9.	Lancer le pipeline nextflow depuis le run.sh avec la commande ./run.sh 
10.	S'il y a une erreur au point 9, faire la commande suivante : sed -i 's/\r//' run.sh
11.	Relancer le pipeline nextflow avec la commande ./run.sh 
12.	Attendre que ça tourne



Cet outil a été pensé par Ben, Hachez Charles et son équipe. Ils ont stoppé la licence provenant de 
la plateforme CLC Genomics qu'ils n'utilisaient que pour faire des analyses RNA-Seq et Chip-Seq. 
Comme alternative, ils ont utilisé la pipeline Nextflow nf-core RNA-Seq. Mais ils ont voulu 
avoir leur propre pipeline pour faire leur analyses rapidement et simplement, c'est ainsi que 
ce projet est nait.  
