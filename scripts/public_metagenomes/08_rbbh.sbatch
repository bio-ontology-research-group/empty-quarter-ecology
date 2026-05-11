#!/bin/bash
#SBATCH --job-name=betA_rbbh
#SBATCH --partition=debug
#SBATCH --exclude=node001,node004,node006,node007
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/rbbh_%j.out
#SBATCH --error=/data/emptyquarter/ecology-paper-runs/public_metagenomes/logs/rbbh_%j.err
set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
WORK=/data/emptyquarter/ecology-paper-runs/public_metagenomes
MSA=$WORK/msa
cd $WORK
echo [concat proteomes] $(date)
# Concatenate all MAG .faa with genome_id prefix in headers for attribution
if [ ! -s $MSA/all_gem_proteins.faa ]; then
  > $MSA/all_gem_proteins.faa
  for f in orfs/*.faa; do
    gid=$(basename $f .faa)
    awk -v g=$gid '/^>/ {sub(/^>/, ">"g"|"); print; next} {print}' $f >> $MSA/all_gem_proteins.faa
  done
fi
grep -c '^>' $MSA/all_gem_proteins.faa
echo [diamond makedb] $(date)
diamond makedb --in $MSA/all_gem_proteins.faa --db $MSA/all_gem.dmnd --quiet
echo [query set] $(date)
# Build query: 5 EQ CSP1-2 betA + E. coli betA
python3 <<PY
import re
from pathlib import Path
MSA=Path('/data/emptyquarter/ecology-paper-runs/public_metagenomes/msa')
out=[]
# EQ CSP1-2 betA sequences from candidate_betA.faa
with open(MSA/'candidate_betA.faa') as fh:
    name=None; buf=[]
    for line in fh:
        if line.startswith('>'):
            if name and name.startswith('EQ_CSP12'): out.append(('>'+name+'\n'+''.join(buf)))
            name=line[1:].strip(); buf=[]
        else:
            buf.append(line.strip())
    if name and name.startswith('EQ_CSP12'): out.append(('>'+name+'\n'+''.join(buf)))
# E. coli betA
with open(MSA/'ref_betA_swissprot.faa') as fh:
    content=fh.read().split('>')
    for r in content:
        if 'BETA_ECOLI' in r:
            out.append('>REF_Ecoli_betA\n'+r.partition('\n')[2].replace('\n','').strip()); break
open(MSA/'rbbh_query.faa','w').write('\n'.join(out)+'\n')
print('Query sequences:', len(out))
PY
grep -c '^>' $MSA/rbbh_query.faa
echo [diamond blastp] $(date)
diamond blastp --query $MSA/rbbh_query.faa --db $MSA/all_gem.dmnd   --out $MSA/rbbh.tsv --outfmt 6 qseqid sseqid pident length evalue bitscore qcovhsp scovhsp   --threads 16 --evalue 1e-10 --max-target-seqs 5000 --quiet
wc -l $MSA/rbbh.tsv
echo [reverse BLAST vs SwissProt refs] $(date)
diamond makedb --in $MSA/ref_betA_swissprot.faa --db $MSA/ref_betA.dmnd --quiet
# For each EQ CSP1-2 query, what's the top SwissProt hit?
diamond blastp --query $MSA/rbbh_query.faa --db $MSA/ref_betA.dmnd   --out $MSA/eq_vs_swissprot.tsv --outfmt 6 qseqid sseqid pident length evalue bitscore qcovhsp scovhsp   --threads 16 --evalue 1e-5 --quiet
echo [DONE] $(date)
head -5 $MSA/rbbh.tsv
echo; echo ==EQ CSP1-2 betA vs SwissProt betA==; head -30 $MSA/eq_vs_swissprot.tsv
