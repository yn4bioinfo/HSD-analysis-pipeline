# HSD-analysis-pipeline

This repository contains core scripts used for analysing duplicated gene copies associated with human-specific segmental duplication genes (HSDs).

The pipeline is currently under active development. It consists of four main modules:

1. 4GTF: converts BLAST-derived transcript/genome alignments into GTF-like duplicated-locus annotations.
2. retroclean: filters repeat-mediated retrocopy-like or fragmentary loci using RepeatMasker annotation.
3. DFM / DupFamMaker: constructs broad duplication families from transcript-sharing evidence.
4. TEM / TopoExonMapper: resolves within-family transcript/exon architecture and copy-level structural relationships.

Large genome-wide input files and full intermediate outputs are not included in this repository because of file size and data management constraints.
