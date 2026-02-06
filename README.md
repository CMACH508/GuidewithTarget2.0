# GuidewithTarget2.0

GuideWithTarget 2.0 is an optimized bioinformatics pipeline designed to identify and validate PIWIL3-associated piRNA cleavage events. Significant upgrades from version 1.0 include high-throughput batch processing for biological replicates and enhanced handling of sequencing duplicates, moving beyond one-by-one execution to a robust, integrated workflow. It bridges small RNA-seq (IP vs. Input) and degradome-seq (WT vs. Mutant) to pinpoint authentic cleavage signatures across the transcriptome and transposable elements with higher statistical power.

---

## Usage

To run the pipeline, ensure your environment is configured (e.g., Conda environment with `bowtie`, `bedtools`, and `samtools`) and execute the main bash script:

Configuration Parameters
The following parameters are defined within the script to control the behavior of the batch analysis:





```bash
# Edit PATHWAY DECLARATIONS inside the script first

Bowtie 1.0.0 – Short read aligner for hierarchical mapping.
Samtools 1.5 – Utilities for manipulating alignments.
Bedtools 2.27.1 – Genome arithmetic for cleavage site extraction.
Blast 2.12.0 – BLAST+ suite for target binding analysis.
Cutadapt 2.6 – Adapter trimming and quality control.
Seqtk 1.3 – Tool for processing sequences in FASTA/Q formats.
Python 3.13.9 – Core scripting language (includes Biopython, Pandas).
R 4.2.0 – Statistical analysis (includes ggplot2, dplyr).
Perl 5.26.2 – Legacy script support and text processing.

# Parameters Declarations
cleavage / nocleavage – Group prefixes for experimental conditions (e.g., WT vs MUT).
num_of_rep_cleav – Number of biological replicates for the cleavage group.
logfc_thres – Minimum Log2 Fold Change to consider a site enriched.
pval_thres – P-value cutoff (e.g., 0.05) for statistical validation across replicates.
set_range – Nucleotide window size (e.g., 20) flanking the cleavage site.
keep_zeroMUT – Set to Y to focus on sites with zero background in control groups.
