# PICRUSt2 Functional Analysis Data

This directory contains the results of the PICRUSt2 functional analysis performed on the 16S amplicon sequencing data from the Rub al-Khali expedition.

## Source Information
- **Remote Path:** `dragon:/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/amplicon_16S/novaseq_14_07_25/final_analysis/functional_analysis/PICRUSt2/`
- **Date Downloaded:** 2026-02-10

## Files
- **path_abun_unstrat.tsv**: MetaCyc pathway abundances (unstratified).
- **path_abun_unstrat_descriptions.tsv**: Descriptions for the MetaCyc pathways.
- **path_abun_unstrat_relative_pct.tsv**: Relative abundance (%) of MetaCyc pathways.
- **EC_predicted.tsv**: Enzyme Commission (EC) number predicted abundances (at sequence level).
- **KO_predicted.tsv**: KEGG Ortholog (KO) predicted abundances (at sequence level).
- **metagenome_pred_metagenome_unstrat.tsv**: Unstratified metagenome functional predictions (EC abundances per sample).
- **ko_pred_metagenome_unstrat.tsv**: Unstratified metagenome functional predictions (KO abundances per sample).
- **sample_mapping.tsv**: Mapping between PICRUSt2 column IDs and standard project Sample IDs.

## Sample Mapping
The PICRUSt2 column headers follow the format `eNNNN_SampleID` (e.g., `e0325_10Dr2`). 
The `sample_mapping.tsv` file provides a direct link between these IDs and the parsed `sample_id`.

| Field | Description |
|-------|-------------|
| picrust2_id | The column header in the TSV files |
| sequence_index | The numeric prefix after 'e' |
| sample_id | The standard Sample Name used in the project |