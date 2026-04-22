#!/usr/bin/env bash
#SBATCH --job-name=t1-kofam
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/logs/t1_kofam.%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/logs/t1_kofam.%j.err
#SBATCH --partition=debug
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --exclude=node001,node004,node007

set -u -o pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics

WD=/data/emptyquarter/ecology-paper-runs/t1_mechanism/kofam
mkdir -p $WD
cd $WD

KOFAM_HMM=/data/databases/anvio_kegg_full/HMMs/Kofam.hmm
KOLIST=/data/databases/anvio_kegg_full/ko_list.txt

# Target KOs for the osmoprotectant mechanism
cat > target_kos.txt << 'TKOEOF'
K00697	otsA	trehalose-6-phosphate synthase
K01087	otsB	trehalose-6-phosphate phosphatase
K05343	treS	trehalose synthase
K06044	treY	malto-oligosyltrehalose synthase
K01236	treZ	malto-oligosyltrehalose trehalohydrolase
K13057	treT	trehalose glycosyltransferring synthase
K00108	betA	choline dehydrogenase
K00130	betB	betaine aldehyde dehydrogenase
K00499	cdh	choline dehydrogenase alt
K17755	gbsA	glycine betaine biosynthesis A
K17756	gbsB	glycine betaine biosynthesis B
K06718	ectA	ectoine biosynthesis A
K06720	ectB	ectoine biosynthesis B
K06721	ectC	ectoine biosynthesis C
K10674	ectD	ectoine hydroxylase
K01858	ino1	myo-inositol-1-phosphate synthase
K05845	opuA	glycine betaine/proline ABC
K05846	opuB	choline ABC
K05847	opuC	glycine betaine/carnitine/choline
K05848	opuD	general osmoprotectant
K03451	proP	proline/glycine betaine symporter
K02168	betT	high-affinity choline
K02000	proV	proU ATP-binding
K02001	proW	proU membrane
K02002	proX	proU substrate-binding
K02176	treA	trehalase
K02779	treB	PTS trehalose permease
K01226	treC	trehalose-6-phosphate hydrolase
K03486	treR	trehalose operon regulator
TKOEOF

# Extract trusted thresholds and the HMMs themselves
awk 'NR>1 {print $1"\t"$2}' $KOLIST > ko_thresholds.tsv

# Extract only the target KOs from the big HMM database
echo "[t1-kofam] extracting $(wc -l < target_kos.txt) target KO HMMs ..."
> target.hmm
for KO in $(awk '{print $1}' target_kos.txt); do
  # hmmfetch requires an index
  :
done
# Use hmmfetch with an index
if [ ! -f ${KOFAM_HMM}.ssi ]; then
  hmmfetch --index $KOFAM_HMM 2>&1 | tail -3
fi
awk '{print $1}' target_kos.txt > ko_list.txt
hmmfetch -f $KOFAM_HMM ko_list.txt > target.hmm 2>/dev/null
NKO=$(grep -c '^NAME' target.hmm)
echo "[t1-kofam] extracted $NKO KO HMMs (expected 29)"

# Pre-compress
hmmpress -f target.hmm 2>&1 | tail -3

# Define genomes + .faa files
declare -a GENOMES=()
# 4 original CSP1-2 MAGs (Prodigal-called from checkm2)
for mag in V27Dr2__SemiBin_73 V30PRr1__SemiBin_7 V32PRr1__SemiBin_26 V38PRr3__SemiBin_38; do
  FAA=/data/emptyquarter/ecology-paper-runs/mags/csp_checkm2/out/protein_files/${mag}.faa
  [ -s "$FAA" ] && GENOMES+=("CSP12_${mag}|$FAA")
done
# 10 EQ-native from co-assembly (Prokka .faa)
for d in /data/emptyquarter/ecology-paper-runs/csp_mag/prokka_out/*/; do
  name=$(basename $d)
  FAA=$(ls $d*.faa 2>/dev/null | head -1)
  [ -s "$FAA" ] && GENOMES+=("EQnative_${name}|$FAA")
done
# 7 NCBI RefSeq (already have _protein.faa.gz, decompress if needed)
for d in /data/emptyquarter/ecology-paper-runs/t1_mechanism/refseq/genomes/*/; do
  genus=$(basename $d | cut -d_ -f1)
  FAA=$(ls $d*_protein.faa 2>/dev/null | head -1)
  [ -z "$FAA" ] && FAA=$(ls $d*_protein.faa.gz 2>/dev/null | head -1 | xargs -I G sh -c 'gunzip -kf G; echo "${G%.gz}"')
  [ -s "$FAA" ] && GENOMES+=("RefSeq_${genus}|$FAA")
done

echo
echo "[t1-kofam] Genomes to scan: ${#GENOMES[@]}"
for G in "${GENOMES[@]}"; do echo "  $G"; done

# Run hmmsearch per genome
mkdir -p per_genome
> all_hits.tsv
echo -e "genome\tKO\tprotein_id\tscore\tevalue\tabove_threshold" > all_hits.tsv

for G in "${GENOMES[@]}"; do
  IFS='|' read -r NAME FAA <<< "$G"
  OUT=per_genome/${NAME}.tbl
  if [ ! -s "$OUT" ]; then
    hmmsearch --cpu 16 --tblout $OUT target.hmm "$FAA" > /dev/null 2>&1
  fi
  N_RAW=$(awk '/^[^#]/ {print}' $OUT | wc -l)
  # Filter by trusted threshold from ko_thresholds.tsv
  # tblout format: target_name acc query_name acc E-value score bias ...
  python3 - << PYEOF
import re
thr = {}
with open("ko_thresholds.tsv") as fh:
    for line in fh:
        f = line.rstrip().split('\t')
        if len(f) >= 2 and f[1] and f[1] != "-":
            try: thr[f[0]] = float(f[1])
            except: pass
rows = []
with open("$OUT") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.split()
        if len(f) < 6: continue
        prot, _, ko, _, evalue, score, *_ = f
        try:
            score = float(score); evalue = float(evalue)
        except: continue
        above = (ko in thr) and (score >= thr[ko])
        rows.append((ko, prot, score, evalue, above))
with open("all_hits.tsv","a") as fh:
    for ko, prot, score, evalue, above in rows:
        fh.write(f"$NAME\t{ko}\t{prot}\t{score}\t{evalue}\t{int(above)}\n")
PYEOF
  n_above=$(awk -F'\t' -v g="$NAME" '$1==g && $6==1' all_hits.tsv | wc -l)
  echo "  $NAME: $N_RAW raw hits, $n_above above trusted threshold"
done

echo
echo "[t1-kofam] Summary written to all_hits.tsv"
wc -l all_hits.tsv

# Pretty report
python3 << 'PYEOF'
import pandas as pd
df = pd.read_csv("all_hits.tsv", sep="\t")
ko_labels = {}
with open("target_kos.txt") as fh:
    for line in fh:
        f = line.rstrip().split('\t')
        if len(f) >= 3:
            ko_labels[f[0]] = f[1]

# Count passing hits per (genome, KO)
agg = df[df.above_threshold==1].groupby(["genome","KO"]).size().unstack(fill_value=0)
# Reorder cols by original target_kos order
ordered = [k for k in list(ko_labels) if k in agg.columns]
agg = agg[ordered]
# Rename to gene names
agg.columns = [ko_labels[k] for k in agg.columns]
agg.to_csv("t1_kofam_gene_counts.tsv", sep="\t")

print("\n===== T1 KOfam (trusted-threshold) gene counts =====")
# Split into biosynth and uptake
bio = ['otsA','otsB','treS','treY','treZ','treT','betA','betB','cdh','gbsA','gbsB','ectA','ectB','ectC','ectD','ino1']
upt = ['opuA','opuB','opuC','opuD','proP','betT','proV','proW','proX','treA','treB','treC','treR']
bio_cols = [c for c in bio if c in agg.columns]
upt_cols = [c for c in upt if c in agg.columns]
print("\n--- BIOSYNTHESIS ---")
print(agg[bio_cols].to_string())
print("\n--- UPTAKE ---")
print(agg[upt_cols].to_string())
PYEOF

echo "=== DONE ==="
