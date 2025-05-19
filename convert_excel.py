# -*- coding: utf-8 -*-
#!/usr/bin/env python3


# convert_excel.py
import sys
import pandas as pd

# Récupérer les fichiers
input_file = sys.argv[1]
output_file = sys.argv[2]

# Lire le fichier
with open(input_file, 'r') as f:
    lines = f.readlines()

# Trouver l'index du header (première ligne sans '#')
for i, line in enumerate(lines):
    if not line.startswith('#'):
        header_index = i
        break

# Charger les données en utilisant pandas
df = pd.read_csv(input_file, sep='\t', skiprows=header_index)

# Sauvegarder le fichier Excel
df.to_excel(output_file, index=False)

print(f"Fichier converti : {output_file}")
