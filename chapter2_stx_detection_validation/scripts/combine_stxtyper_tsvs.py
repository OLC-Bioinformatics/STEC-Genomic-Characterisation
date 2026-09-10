from pathlib import Path
import pandas as pd

folder = Path("stxtyper_results")   # folder containing all TSVs

dfs = []

for file in sorted(folder.glob("*.tsv")):
    if file.name.startswith("._"):
        continue

    df = pd.read_csv(file, sep="\t")

    # Rename #name to Genome
    df = df.rename(columns={"#name": "Genome"})

    dfs.append(df)

combined = pd.concat(dfs, ignore_index=True)

combined.to_csv("stxtyper_combined.tsv", sep="\t", index=False)

print(f"Combined {len(dfs)} files.")
print(f"Total rows: {len(combined)}")