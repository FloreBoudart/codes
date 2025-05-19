# -*- coding: utf-8 -*-
#!/usr/bin/env python3
import sys
import re

if len(sys.argv) < 4:
    print("Usage: enrich_counts.py <gtf_file> <counts_file> <output_file>")
    sys.exit(1)

#récuperer les fichiers 
gtf_file = sys.argv[1]
counts_file = sys.argv[2]
output_file = sys.argv[3]

gene_info = {}

# Lire le fichier GTF et extraire les informations
with open(gtf_file, 'r') as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        columns = line.strip().split("\t")
        if len(columns) > 8 and columns[2] == "gene":
            attributes = columns[8]
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            description_match = re.search(r'description "([^"]+)"', attributes)
            gbkey_match = re.search(r'gbkey "([^"]+)"', attributes)
            gene_biotype_match = re.search(r'gene_biotype "([^"]+)"', attributes)

            if gene_id_match:
                gene_id = gene_id_match.group(1)
                gene_info[gene_id] = {
                    "description": description_match.group(1) if description_match else "",
                    "gbkey": gbkey_match.group(1) if gbkey_match else "",
                    "gene_biotype": gene_biotype_match.group(1) if gene_biotype_match else ""
                }


with open(counts_file, 'r') as counts, open(output_file, 'w') as output:
    first_line = counts.readline().strip()
    
    if first_line.startswith("#"):
        output.write(f"{first_line}\n")
        second_line = counts.readline().strip()
    else:
        second_line = first_line

    new_header = f"{second_line}\tdescription\tgbkey\tgene_biotype"
    output.write(f"{new_header}\n")
    
    for line in counts:
        if line.startswith("#"):
            output.write(line)
        else:
            columns = line.strip().split("\t")
            gene_id = columns[0]

            description = gene_info.get(gene_id, {}).get("description", "")
            gbkey = gene_info.get(gene_id, {}).get("gbkey", "")
            gene_biotype = gene_info.get(gene_id, {}).get("gene_biotype", "")

            output.write(f"{line.strip()}\t{description}\t{gbkey}\t{gene_biotype}\n")

