import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

import argparse
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser(description="Build Bray-Curtis distances for author-supplied plate positions.")
parser.add_argument("--positions", type=Path, required=True)
parser.add_argument("--metadata", type=Path, required=True)
parser.add_argument("--features", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
POSITION_FILE, META_FILE, FEATURE_TABLE = args.positions, args.metadata, args.features
WELLS_OUT, BC_OUT = args.output/"wells.tsv", args.output/"bc_matrix.tsv"

# --- Build well table ---
pos = pd.read_csv(POSITION_FILE)
pos["sid"] = "e" + pos["Sample"].str.split("_").str[:2].str.join("_")
raw_positions = pos.copy()
if not pos.Position.str.fullmatch(r"[A-H](?:[1-9]|1[0-2])").all():
    raise ValueError("Invalid 96-well position")
if pos.duplicated(["Plate", "Position"]).any():
    raise ValueError("A physical well has multiple input rows")
pos = pos.drop_duplicates("sid", keep="first")  # libraries plated twice: keep the first listed well
pos["row"] = pos["Position"].str[0].map(lambda c: ord(c) - ord("A"))
pos["col"] = pos["Position"].str[1:].astype(int)
meta = pd.read_csv(META_FILE, sep="\t")
wells = pos.merge(meta, left_on="sid", right_on="profile", how="inner", validate="one_to_one")
pos.loc[~pos.sid.isin(wells.sid)].to_csv(args.output/"unmapped_positions.tsv", sep="\t", index=False)
raw_positions.loc[raw_positions.duplicated("sid", keep=False)].to_csv(args.output/"repeated_profiles.tsv", sep="\t", index=False)
if wells[["lat", "lon", "site", "compartment"]].isna().any().any():
    raise ValueError("Mapped wells lack required scientific metadata")
wells = wells[["sid", "Plate", "Position", "row", "col", "site", "compartment", "lat", "lon", "dna_conc", "dna_kit"]]
print(f"{len(wells)} wells mapped, {wells.sid.nunique()} unique samples, dropped {len(pos) - len(wells)} unmapped")

# --- Stream feature table, keeping only plated samples ---
sids = sorted(wells.sid.unique())
with open(FEATURE_TABLE) as f:
    f.readline()
    header = f.readline().rstrip("\n").split("\t")
usecols = [0] + [header.index(s) for s in sids]
chunks = []
dtypes = {s: np.float32 for s in sids}
for chunk in pd.read_csv(FEATURE_TABLE, sep="\t", skiprows=1, usecols=usecols, index_col=0, dtype=dtypes, chunksize=50000):
    chunk = chunk.loc[chunk.sum(axis=1) > 0]  # drop features absent from all plated samples
    chunks.append(chunk)
ft = pd.concat(chunks)[sids]
print(f"feature table subset: {ft.shape[0]} features x {ft.shape[1]} samples")

# --- Bray-Curtis on relative abundance ---
X = ft.T.to_numpy(dtype=np.float64)
depth = X.sum(axis=1)
print(f"read depth: min={depth.min():.0f} median={np.median(depth):.0f} max={depth.max():.0f}")
X = X / depth[:, None]
bc = squareform(pdist(X, metric="braycurtis"))
bc_df = pd.DataFrame(bc, index=sids, columns=sids)

# --- Save ---

wells.to_csv(WELLS_OUT, sep="\t", index=False)
bc_df.to_csv(BC_OUT, sep="\t")
print(f"saved {WELLS_OUT}, {BC_OUT}")


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()

provenance = {
    "positions_sha256": digest(POSITION_FILE),
    "metadata_sha256": digest(META_FILE),
    "features_sha256": digest(FEATURE_TABLE),
    "source_position_rows": len(raw_positions),
    "unique_source_profiles": len(pos),
    "repeated_profile_rows_removed": len(raw_positions)-len(pos),
    "mapped_profiles": len(wells),
    "unmapped_profiles": len(pos)-len(wells),
    "features_present": len(ft),
    "read_depth_min": float(depth.min()),
    "read_depth_median": float(np.median(depth)),
    "read_depth_max": float(depth.max()),
    "duplicate_profile_rule": "first listed well",
    "distance": "Bray-Curtis on relative abundance of all canonical ASVs",
}
(args.output/"input_provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
