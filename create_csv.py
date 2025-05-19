# -*- coding: utf-8 -*-
#!/usr/bin/env python3

import sys
import pandas as pd
from pathlib import Path

input_files = [Path(f) for f in sys.argv[1:]]
if len(input_files) < 1:
    sys.exit("Aucun fichier fourni a combiner")

# Lire le premier fichier (base)
df_main = pd.read_csv(input_files[0], sep='\t', comment='#')
df_main = df_main.iloc[:, [0, -1]]  # 1ère colonne = GeneID, dernière = counts
df_main.columns = ['GeneID', input_files[0].stem]

# Ajouter les autres fichiers
for file in input_files[1:]:
    df_next = pd.read_csv(file, sep='\t', comment='#')
    df_main[file.stem] = df_next.iloc[:, -1]

# Sauvegarder
df_main.to_csv("counts_combined.csv", index=False)
