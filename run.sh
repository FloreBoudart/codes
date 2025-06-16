#!/usr/bin/env bash

nextflow run pipeline_rna.nf \
    --reads "samplesheet.csv" \
    --fasta "/chemin/vers/fichier.fasta" \
    --gtf "/chemin/vers/fichier.gtf" \
    --outdir "results" \
    -profile conda \
    -resume