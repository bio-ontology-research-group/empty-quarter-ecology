#!/bin/bash
set -euo pipefail
source /storage/miniforge3/etc/profile.d/conda.sh
conda activate metagenomics
WORK=/data/emptyquarter/ecology-paper-runs/public_metagenomes
MSA=$WORK/msa
mkdir -p $MSA $MSA/eq_faa
# 1) Re-run Prodigal on the 5 EQ CSP1-2 MAGs (4 single + 1 co-assembly re-done for consistency)
EQ_DIR=/data/emptyquarter/ecology-paper-runs/mags/csp_checkm2/input
COAS=/data/emptyquarter/ecology-paper-runs/csp_mag/refine_370_v3/SemiBin_370_HQ_final.fa
for f in $EQ_DIR/V27Dr2__SemiBin_73.fna $EQ_DIR/V30PRr1__SemiBin_7.fna $EQ_DIR/V32PRr1__SemiBin_26.fna $EQ_DIR/V38PRr3__SemiBin_38.fna; do
  id=$(basename $f .fna)
  out=$MSA/eq_faa/CSP12_$id.faa
  if [ ! -s $out ]; then
    prodigal -p meta -q -i $f -a $out -o /dev/null
  fi
done
# Co-assembly HQ MAG
if [ -s $COAS ] && [ ! -s $MSA/eq_faa/coassembly_CSP12_SemiBin_370.faa ]; then
  prodigal -p meta -q -i $COAS -a $MSA/eq_faa/coassembly_CSP12_SemiBin_370.faa -o /dev/null
fi
ls -la $MSA/eq_faa/
# 2) hmmsearch K00108 on all 5 EQ MAGs
for f in $MSA/eq_faa/*.faa; do
  id=$(basename $f .faa)
  tbl=$MSA/eq_faa/$id.K00108.tbl
  if [ ! -s $tbl ]; then
    # Build a mini-HMM with K00108 only
    if [ ! -s $MSA/K00108.hmm ]; then
      python3 -c "
import re
with open('/data/emptyquarter/ecology-paper-runs/t1_mechanism/kofam/target.hmm') as fh:
    hmm=fh.read()
blocks=re.split(r'(?=^HMMER3', hmm, flags=re.M)
for b in blocks:
    if 'NAME  K00108' in b:
        open('$MSA/K00108.hmm','w').write(b); break"
    fi
    hmmsearch --tblout $tbl $MSA/K00108.hmm $f > /dev/null
  fi
done
echo ==EQ K00108 hits==
for t in $MSA/eq_faa/*.K00108.tbl; do
  id=$(basename $t .K00108.tbl)
  best=$(grep -v '^#' $t | sort -k6,6gr | head -1)
  echo "$id: $best"
done
