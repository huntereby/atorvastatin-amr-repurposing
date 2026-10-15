#!/usr/bin/env python3
"""Step 05 - Structural similarity across the statin class (Figure 3).

Each statin is represented as an isomeric SMILES in the form in which it is
administered (simvastatin and lovastatin as lactone prodrugs, the others as
open hydroxy acids; data/raw/statin_structures.csv). Every structure is checked
against its PubChem molecular formula before fingerprinting, and the script
stops if any structure fails.

Morgan fingerprints (radius 2, 2,048 bits; ECFP4-equivalent) are generated with
RDKit and compared with the Tanimoto coefficient for all 28 unique pairs. Pairs
are grouped as within type I, within type II, or between families and
summarised descriptively (median, IQR, range).

Outputs (results/tables/):
    statin_tanimoto_matrix.csv
    statin_pairwise_tanimoto.csv
    statin_family_tanimoto_summary.csv

Run from the repository root:
    python scripts/05_structural_similarity.py
"""

import csv
import statistics
import sys
from itertools import combinations
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors

IN_FILE = Path("data", "raw", "statin_structures.csv")
OUT_DIR = Path("results", "tables")
GROUP_ORDER = ["Within Type I", "Within Type II", "Between families"]


def quantile(values, q):
    """Linear-interpolation quantile (R type 7, the R default)."""
    xs = sorted(values)
    pos = (len(xs) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (pos - lo)


with IN_FILE.open(newline="") as fh:
    statins = list(csv.DictReader(fh))

generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fingerprints = {}
for s in statins:
    mol = Chem.MolFromSmiles(s["smiles"])
    if mol is None:
        sys.exit(f"RDKit could not parse {s['name']}")
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if formula != s["formula"]:
        sys.exit(f"Formula check failed for {s['name']}: RDKit {formula} != expected {s['formula']}")
    fingerprints[s["name"]] = generator.GetFingerprint(mol)
print(f"Validated {len(statins)} structures against their molecular formulas")

names = [s["name"] for s in statins]
family = {s["name"]: s["family"] for s in statins}


def tanimoto(a, b):
    return DataStructs.TanimotoSimilarity(fingerprints[a], fingerprints[b])


OUT_DIR.mkdir(parents=True, exist_ok=True)
with (OUT_DIR / "statin_tanimoto_matrix.csv").open("w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow([""] + names)
    for a in names:
        writer.writerow([a] + [repr(tanimoto(a, b)) for b in names])

pairs = []
for a, b in combinations(names, 2):
    group = f"Within {family[a]}" if family[a] == family[b] else "Between families"
    pairs.append({"drug1": a, "drug2": b, "family1": family[a], "family2": family[b],
                  "comparison_group": group, "tanimoto": tanimoto(a, b)})
pairs.sort(key=lambda p: -p["tanimoto"])
with (OUT_DIR / "statin_pairwise_tanimoto.csv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(pairs[0]))
    writer.writeheader()
    writer.writerows(pairs)

with (OUT_DIR / "statin_family_tanimoto_summary.csv").open("w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["comparison_group", "n_pairs", "median_tanimoto", "q1", "q3", "min", "max"])
    for group in GROUP_ORDER:
        vals = [p["tanimoto"] for p in pairs if p["comparison_group"] == group]
        row = [group, len(vals), statistics.median(vals), quantile(vals, 0.25),
               quantile(vals, 0.75), min(vals), max(vals)]
        writer.writerow(row)
        print(f"{group:17s} n={len(vals):2d}  median {row[2]:.2f}  range {row[5]:.2f}-{row[6]:.2f}")

print("Most / least similar pairs:")
for p in pairs[:2] + pairs[-2:]:
    print(f"  {p['drug1']:12s} {p['drug2']:12s} {p['tanimoto']:.2f}")
