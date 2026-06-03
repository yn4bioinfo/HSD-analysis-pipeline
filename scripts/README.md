# Scripts


This directory contains the four main scripts used in the HSD analysis pipeline.
The scripts are arranged in the approximate order of the current analysis workflow:

| Step |       Script        | Purpose |
|------|---------------------|------------------------------------------------------------------|
|  1.  | 1_4GTF.py           | Convert BLAST results into GTF-like duplicated-locus annotations |
|  2.  | 2_retroclean.sh     | Remove repeat-mediated retrocopy-like or fragmentary loci |
|  3.  | 3_DupFamMaker.py    | Build broad duplication families |
|  4.  | 4_TopoExonMapper.py | Analyse exon/station topology within each family |




## 1. `1_4GTF.py`

### Purpose

`1_4GTF.py` converts transcript-to-genome BLAST results into GTF-like duplicated-locus annotations.

This script is designed to process **one transcript query BLAST result at a time**. In large-scale analyses, many transcript-level BLAST result files are processed in parallel using a Slurm array job.

For each query transcript, the script identifies the native locus, duplicated loci, single-exon-like copies, retrocopy-like copies, and repeat-associated cases, and writes GTF-like annotations and summary tables for downstream analyses.

### Main inputs

| Argument                  | Required | Description                                                                                               |
| ------------------------- | -------: | --------------------------------------------------------------------------------------------------------- |
| `-b`, `--blast_result`    |      yes | BLAST result table for one transcript query                                                               |
| `-c`, `--chr_converter`   |      yes | Two-column chromosome label converter table                                                               |
| `-s`, `--species_label`   |      yes | Species/genome label used in output summaries, e.g. `HS_GRCh38`, `HS_CHM13`, `PTR_panTro6`, `GGO_gorGor1` |
| `-n`, `--native_align`    |       no | Native-locus exon alignment table used in cross-species analysis                                          |
| `-t`, `--TEdup_threshold` |       no | Threshold for TE-associated duplicated regions; default: `100`                                            |

### BLAST result format

The `--blast_result` file is expected to be a tabular BLAST output (BLASTN -outfmt 6), without a header.

The expected columns are:

| Column | Description                              |
| ------ | ---------------------------------------- |
| 1      | query accession / transcript information |
| 2      | subject accession                        |
| 3      | percent identity                         |
| 4      | alignment length                         |
| 5      | mismatches                               |
| 6      | gap opens                                |
| 7      | query start                              |
| 8      | query end                                |
| 9      | subject start                            |
| 10     | subject end                              |
| 11     | e-value                                  |
| 12     | bit score                                |

This corresponds to a BLAST tabular output similar to:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

### Chromosome converter format

The `--chr_converter` file is a whitespace-delimited two-column table that converts subject accession IDs to chromosome labels.

Example:

```text
NC_000001.11  1
NC_000002.12  2
NC_000003.12  3
NC_000023.11  X
NC_000024.10  Y
```

The first column is the subject accession ID used in the BLAST database.
The second column is the chromosome label used in the downstream GTF-like outputs.

### Minimal command for one transcript query

```bash
python 1_4GTF.py \
  --blast_result path/to/ENSTxxxx_blast_result.tsv \
  --chr_converter path/to/chr_label_converter.tsv \
  --species_label HS_GRCh38
```

Equivalent short-option form:

```bash
python 1_4GTF.py \
  -b path/to/ENSTxxxx_blast_result.tsv \
  -c path/to/chr_label_converter.tsv \
  -s HS_GRCh38
```

### Large-scale execution with a Slurm array job

In the actual genome-wide analysis, this script was executed as a Slurm array job.

The array job reads a manifest file containing paths to individual transcript-level BLAST result files. Each array task selects one line from the manifest and passes that file to `1_4GTF.py` via `--blast_result`.

Example manifest file:

```text
path/to/blast_result_ENST000001.tsv
path/to/blast_result_ENST000002.tsv
path/to/blast_result_ENST000003.tsv
```

Example Slurm array usage:

```bash
#!/bin/bash
#SBATCH --array=1-144
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=5
#SBATCH --output=logs/4GTF_%A_%a.out
#SBATCH --error=logs/4GTF_%A_%a.err

set -euo pipefail
shopt -s nullglob

manifest="path/to/blast_result_file_list.tsv"
chr_converter="path/to/chr_label_converter_HS_GRCh38.tsv"
species_label="HS_GRCh38"

mapfile -t blast_files < "${manifest}"
blast_file="${blast_files[$SLURM_ARRAY_TASK_ID-1]}"

python 1_4GTF.py \
  -b "${blast_file}" \
  -c "${chr_converter}" \
  -s "${species_label}"
```

### Main outputs

The script writes output files to the current working directory.

Major output files include:

| Output pattern                       | Description                                                         |
| ------------------------------------ | ------------------------------------------------------------------- |
| `GTF_from_BLASTresult_*.gtf`         | GTF-like annotation of duplicated loci                              |
| `single_GTF_from_BLASTresult_*.gtf`  | GTF-like annotation for single-exon-like cases                      |
| `dups_*.tsv`                         | Table of duplicated hits/loci for the query transcript              |
| `dups_area_*.tsv`                    | Query-exon/duplicated-area mapping table used by downstream scripts |
| `statistics_*.tsv`                   | Copy-number and classification summary for the query transcript     |
| `Repeat_*.tsv`                       | Repeat-associated case report, when applicable                      |
| `TE_*.tsv`                           | TE-associated duplicated-hit report, when applicable                |
| `native_locus_exon_alignments_*.tsv` | Native-locus exon alignment table, when applicable                  |

### Notes

* This script processes one transcript-level BLAST result at a time.
* For genome-wide runs, use a Slurm array job or an equivalent workflow manager.
* The manifest file used in the Slurm array job is not an input to `1_4GTF.py` itself; it is used by the batch script to select the `--blast_result` file for each task.
* Because output files are written to the current working directory, it is recommended to run the script in a dedicated output directory or to move output files after each batch run.
* Absolute paths used in local or supercomputer environments should be replaced with project-relative paths in public documentation.
