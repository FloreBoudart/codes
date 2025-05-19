#!/usr/bin/env nextflow

// Chargement des paramètres principaux
params.samplesheet = "samplesheet.csv"
params.reads  = null
params.fasta  = null
params.gtf    = null
params.outdir = null
params.multiqc_input = "${params.outdir}"
params.outdir_abs = file(params.outdir).toAbsolutePath().toString()
params.threads = 4

// Générer le chemin d'indexdir automatiquement à partir du chemin du fichier FASTA
def fastaFile = file(params.fasta)
def fastaBaseDir = fastaFile.getParent()
params.indexdir = "${fastaBaseDir}/index"
	
	
workflow.onComplete {
    if (params.outdir == null || params.outdir == "") {
        // Get path from the first FASTQ in the samplesheet
        def firstLine = file(params.samplesheet).readLines()[1]
        def firstFastq = firstLine.tokenize(',')[1]
        def baseDir = file(firstFastq).getParentFile().getParent()

        // Use that to build the outdir
        def autoOutdir = "${baseDir}/results_config"
        println "Setting output directory to: ${autoOutdir}"
        params.outdir = autoOutdir
    }
}


//FastQC : s'assure de la qualité des données brutes
//récupérer les échantillons dans les deux sens de lecture, donc r1 et r2
//Sortie de fichiers HTML et fichiers bruts au format fastq.gz
process FastQC {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    path "*_fastqc.*"

    script:
    """
    fastqc ${r1} 
	fastqc ${r2}
    """
}

//Fastp : nettoyage des données
process Fastp {
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_R1_clean.fastq.gz"), path("${sample}_R2_clean.fastq.gz")

    script:
    """
    fastp -i ${r1} -I ${r2} -o ${sample}_R1_clean.fastq.gz -O ${sample}_R2_clean.fastq.gz
    """
}

//Trimmomatic : filtre, coupe les adaptateurs
process Trimmomatic {
    publishDir "${params.outdir}/trim", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_R1_trim.fastq.gz"), path("${sample}_R2_trim.fastq.gz")

    script:
    """
    trimmomatic PE -threads 4 ${r1} ${r2} \
        ${sample}_R1_trim.fastq.gz /dev/null \
        ${sample}_R2_trim.fastq.gz /dev/null \
        HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:100
    """
}

//Construction de l'index avec Hisat2
//Résultats enregistré dans genome/index
process BuildIndex {
    publishDir "${params.indexdir}", mode: 'copy'

    input:
    path fasta

    output:
    path "genome_index.*.ht2"

    script:
    """
    hisat2-build ${fasta} genome_index
    """
}

//Hisat2 : fait l'alignement des échantillons sur l'index
process AlignHisat2 {
    publishDir "${params.outdir}/alignments", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)
    each path(index_files)

    output:
    tuple val(sample), path("${sample}.sam")

    script:
    """
    hisat2 -p 4 -x genome_index -1 ${r1} -2 ${r2} -S ${sample}.sam
    """
}

//Samtools view :transforme les fichiers sam en bam
process SamtoolsView {
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("${sample}.bam")

    script:
    """
    samtools view -b ${sam} -o ${sample}.bam
    """
}

//Samtools sort: Trie les fichiers bam
process SamtoolsSort {
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.sorted.bam")

    script:
    """
    samtools sort ${bam} -o ${sample}.sorted.bam
    """
}

//Samtools index: index les fichiers bam
process SamtoolsIndex {
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    path("${bam.simpleName}.bai")

    script:
    """
    samtools index ${bam} > ${bam.simpleName}.bai
    """
}

//créer le fichier de statistique d'alignement avec samtools flagstat
process SamtoolsFlagstat {
    publishDir("${params.outdir}/bam", mode: 'copy')

    input:
    tuple val(sample), path(bam)

    output:
    path("${sample}.flagstat")

    script:
    """
    samtools flagstat ${bam} > ${sample}.flagstat
    """
}

//FeaturesCounts : compter les lectures 
process FeatureCounts {
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(sample), path(bam)
    each path(gtf)

    output:
    path("${sample}_counts.txt")

    script:
    """
    featureCounts -T 4 -p -a ${gtf} -o ${sample}_counts.txt ${bam}
    """
}

//FeaturesCounts : création du fichier de statistique d'alignement 
process FeatureCounts_summary {
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(sample), path(bam)
    each path(gtf)

    output:
    path("${sample}_counts.txt.summary")

    script:
    """
    featureCounts -T 4 -p -a ${gtf} -o ${sample}_counts.txt ${bam}
    """
}

//le code python Récupère les fichiers de comptage
//Créé un fichier au format csv avec les noms de gènes et les comptages de tous les échantillons 
process CreateCSV {
    publishDir("${params.outdir}/csv", mode: 'copy')

    input:
    path(counts_txt)

    output:
    path("counts_combined.csv")

    script:
    """
    python3 ${workflow.projectDir}/scripts/create_csv.py ${counts_txt}
    """
}

//Code python qui récupère les fichiers de comptages et rajoute les colonnes
//descriptions, gbkey, gene et gene_biotype qui sont récupéré dans le fichier gtf 
//et les enregistre au format txt 
process CompleteCounts {
    publishDir("${params.outdir}/counts", mode: 'copy')

    input:
    path(counts_txt)
    each path(gtf)

    output:
    path("${counts_txt.baseName}_counts_enriched.txt")

    script:
    """
    python3 ${workflow.projectDir}/scripts/enrich_counts.py ${gtf} ${counts_txt} ${counts_txt.baseName}_counts_enriched.txt
    """
}

//code python qui utilise les fichier txt cité avant 
//les transforme en fichier xlsx qui ressemble à Excel
process ConvertToExcel {
    publishDir("${params.outdir}/excel", mode: 'copy')

    input:
    path(counts_enriched_txt)

    output:
    path("${counts_enriched_txt.baseName}_counts.xlsx")

    script:
    """
    python3 ${workflow.projectDir}/scripts/convert_excel.py ${counts_enriched_txt} ${counts_enriched_txt.baseName}_counts.xlsx
    """
}

//MultiQC : récupère les résultats de tous les outils et les place dans un fichier html
//Fichier qu'on peut ouvrir sur un navigateur internet et visualiser les résultats en graphique ou en tableau 
process MultiQC {
    publishDir("${params.outdir}/multiqc", mode: 'copy')

	input:
    val(params.outdir_abs)

    output:
    path("multiqc_result.html")

    script:
    """
    multiqc ${params.outdir_abs}/* -o . -n multiqc_result.html
    """
}


//Workflow : la partie qui active chaque processus et de les lier pour 
//qu'ils se lancent les uns à la suite des autres 
//Seul MultiQc est différent car il doit attendre les résultats de tous avant de s'activer 
workflow {

    /// 1. Charger les échantillons depuis le fichier CSV
    Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample, file(row.fastq_1), file(row.fastq_2))
        }
        .set { samples_ch }
		
    // QC
    qc_ch = FastQC(samples_ch)

    // Nettoyage et Trim
    cleaned_ch = Fastp(samples_ch)
    trimmed_ch = Trimmomatic(cleaned_ch)
	
	// Indexation du génome
    genome_ch = Channel.fromPath(params.fasta)
    index_ch  = BuildIndex(genome_ch)
	
    // Alignement
    aligned_ch = AlignHisat2(trimmed_ch, index_ch)

    // BAM
    bam_ch        = SamtoolsView(aligned_ch)
    sorted_bam_ch = SamtoolsSort(bam_ch)
    index_bam     = SamtoolsIndex(sorted_bam_ch)
	flagstat_ch = SamtoolsFlagstat(sorted_bam_ch)

    // Quantification
    gtf_ch = Channel.fromPath(params.gtf)
    counts_ch = FeatureCounts(sorted_bam_ch, gtf_ch)
    summary_ch = FeatureCounts_summary(sorted_bam_ch, gtf_ch)
	
	// Fusion des comptages
    csv_counts_ch = CreateCSV(counts_ch.collect())

    // Enrichissement des comptages
    enriched_ch = CompleteCounts(FeatureCounts.out, gtf_ch)

    // Conversion en Excel
    excel_ch = ConvertToExcel(enriched_ch)

    // Rapport MultiQC
	multiqc_ready = qc_ch.mix(cleaned_ch)
                      .mix(trimmed_ch)
                      .mix(aligned_ch)
                      .mix(sorted_bam_ch)
                      .mix(counts_ch)
					  .mix(enriched_ch)
                      .collect()
	multiqc_ch = MultiQC(multiqc_ready)
}
