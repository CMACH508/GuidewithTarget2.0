#!/bin/bash
#Go Go Go!
#😄😄😄😄😄😄 

# Source conda initialization script
source ~/miniconda3/etc/profile.d/conda.sh  # Adjust path if needed

# Activate your environment
conda activate zlin3


: << '-------------------------------------------->PREPARATIONS<--------------------------------------------'
Prepare for the sequencing data, 

cd $transcriptome_dir
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/10036/103/GCF_017639785.1_BCM_Maur_2.0/GCF_017639785.1_BCM_Maur_2.0_genomic.gff.gz
gunzip GCF_017639785.1_BCM_Maur_2.0_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/10036/103/GCF_017639785.1_BCM_Maur_2.0/GCF_017639785.1_BCM_Maur_2.0_rna.fna.gz
gunzip GCF_017639785.1_BCM_Maur_2.0_rna.fna
python "${script_dir}clean_fasta.py" GCF_017639785.1_BCM_Maur_2.0_rna.fna > GCF_017639785.1_BCM_Maur_2.0_rna_clean.fna
python "${script_dir}clean_transposon.py" -i GCF_017639785.1_BCM_Maur_2.0_DB-families.fa -o GCF_017639785.1_BCM_Maur_2.0_DB-families_clean.fa
python "${script_dir}transform_gff_into_table.py" GCF_017639785.1_BCM_Maur_2.0_genomic.gff > transcript_table_anno_2.txt
bowtie-build --threads 40 GCF_017639785.1_BCM_Maur_2.0_rna_clean.fna transcript
bowtie-build --threads 40 GCF_017639785.1_BCM_Maur_2.0_DB-families_clean.fa transposon
python "${script_dir}seqlength_stat.py"  GCF_017639785.1_BCM_Maur_2.0_rna.fna BCM_Macur_seq_stat.txt
python "${script_dir}seqlength_stat.py"  GCF_017639785.1_BCM_Maur_2.0_DB-families_clean.fa BCM_Macur_seq_tp_stat.txt

cd $genome_dir
python "${script_dir}clean_fasta.py" BCM_Maur_2.0.fasta > BCM_Maur_2.0_clean.fasta
bowtie-build --threads 40 BCM_Maur_2.0_clean.fasta genome
python "${script_dir}seqlength_stat.py"  BCM_Maur_2.0_clean.fasta BCM_Macur_genome_stat.txt

cd $transcriptome_dir
RepeatMasker -pa 100 -lib GCF_017639785.1_BCM_Maur_2.0_DB-families_clean.fa -a -excln ../genome/BCM_Maur_2.0_clean.fasta
cd $genome_dir
nohup perl /home/user/Downloads/RepeatMasker/util/calcDivergenceFromAlign.pl \
    -s BCM_Maur_2.0_clean.summary \
    -noCpGMod \
    BCM_Maur_2.0_clean.fasta.align &
-------------------------------------------->PREPARATIONS<--------------------------------------------







: << '-------------------------------------------->PATHWAY DECLARATIONS<--------------------------------------------'
Adjust Below Codes for Established Your own working directory
Setup the experiment setting accordingly
-------------------------------------------->PATHWAY DECLARATIONS<--------------------------------------------

parent_dir="/data/zenglin/goldenham/degradome/"
download_dir="${parent_dir}download/"
script_dir="${parent_dir}script/"
rawdata_dir="/data/zenglin/P3code/rawdata_submit/goldenhamster/"
backup_dir="${parent_dir}backup/"
wkdir="${parent_dir}wkdir/"
transcriptome_dir="${parent_dir}transcriptome/"
genome_dir="${parent_dir}genome/"
guide="PIWIL3_associated"
species="goldenhamster"
#genome_file="mm10" # in this pipeline, needless to define
protein="PIWIL3"
risc="${protein}_${guide}"
reference="transcriptome"
cleavage="WT"
nocleavage="MUT"
data_dir="${wkdir}${protein}_${guide}_${cleavage}_${nocleavage}/"
set_range="20"
num_of_rep_cleav="2"
num_of_rep_nocleav="2"
piRNA_fc_thres=1 # IP enrichement threshold
logfc_thres=0 # cleavage products to align
pval_thres=1 # cleavage products to align
keep_zeroMUT="N" # cleavage products to align



: << '-------------------------------------------->DATA PREPARATIONS<--------------------------------------------'
We have done two-batch data sequencing, we need to combine them into one whole
So decompress data, double check completeness of data, and then concatenate the data so that they are complete!
-------------------------------------------->DATA PREPARATIONS<--------------------------------------------


: << '-------------------------------------------->PROCESS_FILES<--------------------------------------------'
Untar files
Double check completeness
Rename libraries so that names could be concisely and precisely short
Remove adaptors
-------------------------------------------->PROCESS_FILES<--------------------------------------------



echo "Your analysis starts with here..."
echo $rawdata_dir

echo "     _________"
echo "    //   ||  \\"
echo "   //____||___\\"
echo "   |   Start  |"
echo "  _|__________|_"
echo " |              |"
echo " |______________|"
echo "   (O)       (O)"

echo "UnZip, Check md5, Rename, Cut # Adaptors..."
sleep 10
echo "Too many trivial things to do, but necessarily important!"


cd $rawdata_dir
adaptR1="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
adaptR2="GTTCAGAGTTCTACAGTCCGACGATC"
adaptR3="GTTCAGAGTTCTACAGTCC"
for i in $(ls *R1.fastq.gz);
do
    cutadapt -m 15 --length 25 --max-n 0 -g $adaptR2 \
            -b $adaptR3 -a $adaptR1 -j 16  \
            -o ${i/.fastq/.trimmed} $i 1> $i.txt
done

if [ -d "$data_dir" ]; then
    echo "Your directory exists, removing it..."
    rm -r "$data_dir"
fi
mkdir -p "$data_dir"
echo "Directory created successfully!"

mv "${rawdata_dir}${cleavage}"*.trimmed.gz $data_dir
mv "${rawdata_dir}${nocleavage}"*.trimmed.gz $data_dir
mv *.txt $data_dir


# This step implements a rigorous preprocessing pipeline for degradome-seq data to isolate authentic piRNA-mediated cleavage signals from raw sequencing reads. The process begins by navigating to the raw data directory and defining three specific adapter sequences commonly used in degradome library preparation, ensuring comprehensive removal of both 5' and 3' adapter remnants that could interfere with downstream alignment. Each R1 fastq file undergoes systematic adapter trimming using cutadapt with stringent parameters: discarding reads shorter than 15 nucleotides, standardizing read lengths to 25 bases, and eliminating sequences containing ambiguous bases (N) to ensure mapping quality. The pipeline employs iterative adapter removal targeting multiple potential adapter configurations, leveraging 16 processing threads for computational efficiency while generating detailed log files for quality assessment. Following successful adapter trimming, the script performs directory management by either removing or creating a dedicated data storage location, then systematically organizes all processed files—including both cleavage and non-cleavage experimental conditions—into structured directories with accompanying documentation files, establishing a clean, organized dataset ready for subsequent alignment and cleavage site identification analyses.


: << '-------------------------------------------->MAPPING<--------------------------------------------'
Here we map reads to reference
Map to reference both strands are accepted,
Since maybe in some regions are bi-directionally transcribed~
Core part to do, you can of course modify script manaully,
e.g. 1 mismatch tolerated~
bowtie -p 40 -t -v 0 -x "${genome_dir}mm10" $i -S "$i.sam" 2>$i.map.log
# mapped rate so low!
you can use following comd to browse:
skills:
awk '$3 == "ENSMUST00000192833.1" && $4 == 833'  WT_1_R1.trimmed.gz.candidates.sam
# only keep forward strand hits!
-------------------------------------------->MAPPING<--------------------------------------------

cd $data_dir
gzip -dk *.gz


# for transcriptome
# --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)
# -n/--seedmms       max mismatches in seed (can be 0-3, default: -n 2)
# -t/--time          print wall-clock time taken by search phases
# HACK -k <int>           report up to <int> good alignments per read (default: 1) 
# --best             hits guaranteed best stratum; ties broken by quality
# --un <fname>       write unaligned reads/pairs to file(s) <fname>
echo "Mapping To transcriptome..."
for i in $(ls *_R1.trimmed);
do
        echo "Aligning this fastq file now --> $i";
        bowtie -nomaqround -n 1 -p 16 -t -k 1 --best --no-unal --un "$i.unaligned.fastq" \
                --best "${transcriptome_dir}transcript" $i -S "$i.sam" 2>$i.map.log 
        samtools flagstat "$i.sam" > "$i.flagstat.txt"
        sed '/^@/d' "$i.sam" | cut -f 1,2,3,4,10 | awk '($2 == "0")' > "$i.candidates.sam" # forward strand only
done
echo "Transcriptome Alignment Finished!"


echo "Mapping unaligned reads to transposon database..."
for i in $(ls *_R1.trimmed);
do
    base_name=$(basename "$i" _R1.trimmed)
    unaligned_fastq="${base_name}_R1.trimmed.unaligned.fastq"
    
    if [ -f "$unaligned_fastq" ]; then
        echo "Mapping unaligned reads from $base_name to transposon..."
        
        # Map to transposon with double-strand mapping and more mismatches
        # Key differences from transcriptome mapping:
        # - No strand restriction (removed filtering for forward strand only)
        # - More mismatches allowed (-n 2 instead of -n 1)
        # - We might want to allow more alignments per read since transposons can have repeats
        
        bowtie -nomaqround -n 2 -p 16 -t -k 1 --best --no-unal \
                "${transcriptome_dir}transposon" "$unaligned_fastq" \
                -S "${base_name}_transposon.sam" 2>${base_name}_transposon.map.log
        
        # Generate flagstat
        samtools flagstat "${base_name}_transposon.sam" > "${base_name}_transposon.flagstat.txt"
        
        # Extract alignment information (both strands, no strand filtering)
        # For transposon mapping, we might want to keep both forward and reverse alignments
        sed '/^@/d' "${base_name}_transposon.sam" | cut -f 1,2,3,4,10 > "${base_name}_transposon.candidates.sam"
        
        echo "Transposon mapping completed for $base_name"
        echo "Unaligned reads: $(grep -c '^@' "$unaligned_fastq" || echo "0")"
        echo "Transposon aligned: $(samtools view -F 4 "${base_name}_transposon.sam" | wc -l)"
    else
        echo "No unaligned reads file found for $base_name"
    fi
done

echo "Transposon mapping finished!"

# This step implements a tiered mapping strategy to comprehensively characterize degradome-seq reads by sequentially aligning them to distinct genomic features. The pipeline begins by mapping trimmed reads to the transcriptome reference using stringent Bowtie parameters: disabling Maq-style quality rounding, allowing only one mismatch in the seed region (-n 1), and requiring perfect strand orientation while tracking processing time and discarding unaligned reads to a separate file. Following transcriptome alignment, SAMtools flagstat generates alignment statistics, and a custom filtering step extracts only forward-strand alignments (flag 0) as candidate cleavage sites, focusing on conventional miRNA targeting patterns. The unaligned reads from this primary mapping then undergo secondary alignment to a transposon database using relaxed parameters that permit two mismatches (-n 2) and capture both DNA strands, acknowledging the repetitive nature and potential antisense transcription of transposable elements. This hierarchical approach ensures comprehensive read assignment while maintaining appropriate stringency for each target class, generating separate alignment files and statistics for transcriptome and transposon mappings that collectively provide a complete picture of degradome signal distribution across different genomic compartments.


: << '-------------------------------------------->FETCHING_POSITION<--------------------------------------------'
Prepare for Cutting edge of each splicing position!
Fetch the sequence of cutting edge
Para:
 cleav_suffix: cleavage lib e.g. Cleavage_1
 nocleav_suffix: non-cle lib e.g. NoCleavage_1
 data_add: where files located in
 reference: genome
 num_of_rep: number of biological repeats
 set_range: 1-based index, set 1st base cleavage products
as center, fetch region set_range-long right/left side,
e.g. 50, then 51 is center, 1-50 and 52-101 as flanking
Input:
 cleavage_candidate_samfile
Ouput:
 bed files of potential cleavage products
 e.g. chr1_10097054_+ (1-based coordiantes, then first base is 10097054)
-------------------------------------------->FETCHING_POSITION<--------------------------------------------

for file in {WT,MUT}[0-9]*; do mv -- "$file" "$(echo "$file" | sed -E 's/^(WT|MUT)([0-9])/\1_\2/')"; done

# transcriptome relevant cleavage events
echo "Fetching position of transcripts cleaged..."
Rscript "${script_dir}fetch_DegraPosition_rep.R" \
        $cleavage \
        $nocleavage \
        $data_dir \
        $num_of_rep_cleav \
        $num_of_rep_nocleav \
        $set_range \
        $risc \
        "${transcriptome_dir}BCM_Macur_seq_stat.txt"


# transcriptome relevant cleavage events
echo "Fetching position of transcripts cleaged..."
Rscript "${script_dir}fetch_DegraPosition_rep_tp.R" \
        $cleavage \
        $nocleavage \
        $data_dir \
        $num_of_rep_cleav \
        $num_of_rep_nocleav \
        $set_range \
        $risc \
        "${transcriptome_dir}BCM_Macur_seq_tp_stat.txt"


# Filter the products so that the mapping can focus on only these products
# Rscript script_name.R input.txt output.bed cleavage_group non_cleavage_group num_cleavage_reps num_non_cleavage_reps 2 0.05
echo "Filtering transcripts cleaged..."
cd $wkdir
Rscript "${script_dir}filter_products.R" \
        "${data_dir}${risc}_insert_expr.txt" \
        "${data_dir}${risc}_insert_filter_pval${pval_thres}_logfc${logfc_thres}_coding.bed" \
        $cleavage \
        $nocleavage \
        $num_of_rep_cleav \
        $num_of_rep_nocleav \
        "nofilter" \
        $pval_thres \
        $keep_zeroMUT

# Rscript script_name.R input.txt output.bed cleavage_group non_cleavage_group num_cleavage_reps num_non_cleavage_reps 2 0.05
echo "Filtering transcripts cleaged..."
cd $wkdir
Rscript "${script_dir}filter_products.R" \
        "${data_dir}${risc}_insert_expr.tp.txt" \
        "${data_dir}${risc}_insert_filter_pval${pval_thres}_logfc${logfc_thres}_tp.bed" \
        $cleavage \
        $nocleavage \
        $num_of_rep_cleav \
        $num_of_rep_nocleav \
        "nofilter" \
        $pval_thres \
        $keep_zeroMUT


# Following successful alignment of degradome-seq reads to transcriptomic and transposonic references, this pipeline implements a sophisticated two-tiered analysis to identify and validate authentic small RNA-directed cleavage events through comprehensive positional mapping and statistical filtering. The process begins with systematic file reorganization to standardize naming conventions across biological replicates, then executes specialized R scripts that extract precise cleavage site coordinates from forward-strand aligned reads while accounting for sequence length constraints and generating normalized 20-nucleotide windows centered on each putative cleavage position. These scripts process both transcriptome-derived and transposon-derived cleavage events in parallel, implementing rigorous normalization procedures that convert raw read counts to reads-per-million (RPM) while applying appropriate pseudocounts for low-abundance sites, and compiling comprehensive expression matrices that integrate data across multiple biological replicates of both cleavage (WT) and non-cleavage (MUT) conditions. Subsequently, a rigorous differential analysis framework identifies biologically significant PIWIL3-associated cleavage events by comparing normalized degradome signals between conditions, employing stratified statistical testing that accounts for experimental design: for single-replicate conditions, it applies straightforward log₂ fold-change filtering, while for multi-replicate datasets it calculates pooled variances, computes t-statistics, and derives empirical p-values through two-tailed t-distribution testing. The analysis specifically targets down-regulated events in the non-cleavage condition while incorporating an optional stringent filter that retains only cleavage sites completely absent in the mutant background, ensuring detection of PIWIL3-dependent cleavage signatures. Finally, the pipeline extracts genomic coordinates of statistically validated events, formats them into standard BED files, and generates filtered datasets ready for sequence extraction and subsequent alignment with PIWIL3-bound small RNAs for downstream motif discovery and functional validation analyses.



: << '-------------------------------------------->FETCHING_SEQUENCE<--------------------------------------------'
Prepare for Cutting edge of each splicing position!
Fetch the sequence of cutting edge
Para:
 cleav_suffix: cleavage lib e.g. Cleavage_1
 nocleav_suffix: non-cle lib e.g. NoCleavage_1
 data_add: where files located in
 reference: genome
 num_of_rep: number of biological repeats
 set_range: 1-based index, set 1st base cleavage products
as center, fetch region set_range-long right/left side,
e.g. 50, then 51 is center, 1-50 and 52-101 as flanking
Input:
 cleavage_candidate_samfile
Ouput:
 bed files of potential cleavage products
 e.g. chr1_10097054_+ (1-based coordiantes, then first base is 10097054)
-------------------------------------------->FETCHING_SEQUENCE<--------------------------------------------

# extract sequence from transcriptome file
# filter.bed comes from Cleavage In Vivo.ipynb
echo "Extracting Sequence from Transcriptome File:"
bedtools getfasta -s -fi "${transcriptome_dir}GCF_017639785.1_BCM_Maur_2.0_rna.fna" \
        -bed "${data_dir}${risc}_insert_filter_pval${pval_thres}_logfc${logfc_thres}_coding.bed" \
        -name > "${data_dir}${risc}_insert_filter_pval${pval_thres}_logfc${logfc_thres}_coding.fa"

echo "Extracting Sequence from Transcriptome File:"
#  -fi     Input FASTA file
#  -s      Force strandedness. If the feature occupies the antisense,
bedtools getfasta -s -fi "${transcriptome_dir}GCF_017639785.1_BCM_Maur_2.0_DB-families_clean.fa" \
        -bed "${data_dir}${risc}_insert_filter_pval${pval_thres}_logfc${logfc_thres}_tp.bed" \
        -name > "${data_dir}${risc}_insert_filter_pval${pval_thres}_logfc${logfc_thres}_tp.fa"

cat ${data_dir}${risc}_insert_filter_pval${pval_thres}_logfc${logfc_thres}_tp.fa \
    ${data_dir}${risc}_insert_filter_pval${pval_thres}_logfc${logfc_thres}_coding.fa \
    > ${data_dir}${risc}_insert_filter_pval${pval_thres}_logfc${logfc_thres}.fa

# This stage performs targeted sequence extraction from genomic reference databases to retrieve the specific nucleotide contexts surrounding validated cleavage sites for subsequent small RNA binding analysis. Using BEDtools, the pipeline extracts stranded genomic sequences from two distinct reference sources: first from the golden hamster transcriptome FASTA file to capture protein-coding cleavage contexts, and then from a curated transposon database to obtain repetitive element-associated cleavage regions, ensuring strand-specificity through the -s parameter to maintain proper orientation of cleavage products. These extractions are guided by filtered BED files containing statistically significant cleavage events that meet both fold-change and p-value thresholds, generating separate FASTA files for coding and transposon-derived sequences that include comprehensive header information linking each extracted sequence to its original genomic coordinates. Finally, the pipeline concatenates both sequence collections into a unified FASTA file, creating a comprehensive nucleotide dataset encompassing all validated cleavage events across different genomic compartments for downstream alignment with PIWIL3-associated small RNAs and motif discovery analyses.


: << '-------------------------------------------->ESTABLISH_LIBRARY<--------------------------------------------'
Prepare for Cutting edge of each splicing position!
Fetch the sequence of cutting edge
-------------------------------------------->ESTABLISH_LIBRARY<--------------------------------------------
makeblastdb -in "${data_dir}${risc}_insert_filter_pval${pval_thres}_logfc${logfc_thres}.fa" \
  -dbtype nucl \
  -max_file_sz 2GB \
  -out "${data_dir}${risc}_insert_db_filter_pval${pval_thres}_logfc${logfc_thres}"


: << '-------------------------------------------->MAPPING POSITIONS<--------------------------------------------'
Prepare for Cutting edge of each splicing position!
Fetch the sequence of cutting edge
Add align the small to the splicing positions
-------------------------------------------->MAPPING POSITIONS<--------------------------------------------

# align POI 
cd $data_dir
smallRNA_file="${download_dir}PIWIL3_associated_piRNA.fa"
cp $smallRNA_file "${data_dir}${guide}.fa"
echo $smallRNA_file
echo "I have copy this file as smallRNA file"
directory="${data_dir}blast_POI_filter"

if [ -d "$directory" ]; then
        echo "Directory exists. Removing it..."
        rm -rf "$directory"
fi

echo "Creating directory..."
mkdir $directory
cd $directory

blastn -query "${data_dir}${guide}.fa" \
        -num_threads 16 \
        -evalue 10000000000 \
        -db "${data_dir}${risc}_insert_db_filter_pval${pval_thres}_logfc${logfc_thres}" \
        -strand minus \
        -task blastn \
        -word_size 10 \
        -outfmt 0 \
        -ungapped > "${risc}_POI_BP.txt"


# Following cleavage site identification, this pipeline establishes a target sequence database and performs reverse-complement alignment to identify potential PIWIL3-bound piRNA binding sites using a specialized BLAST protocol. First, a BLAST nucleotide database is constructed from filtered cleavage site sequences encompassing both coding transcripts and transposable elements, creating a searchable reference for rapid sequence matching. Subsequently, PIWIL3-associated piRNA sequences are aligned against this database using stringent parameters optimized for small RNA recognition: enforcing reverse-strand orientation (-strand minus) to identify complementary base-pairing patterns, employing a 10-nucleotide word size for sensitive short sequence detection, and utilizing permissive e-value thresholds to capture potential binding interactions. This alignment strategy specifically targets the characteristic antisense binding mechanism of PIWI proteins, generating comprehensive pairing data between enriched piRNAs and statistically validated cleavage sites for subsequent motif analysis and functional validation.


Rscript "${script_dir}Parse_POI.R" \
        "${risc}_POI_BP.txt" \
        $risc
find . -name '*Parsed*' -print0 | \
        xargs -0 cat > \
        "${risc}_parsed_POI_filter.txt"
mv "${risc}_parsed_POI_filter.txt" $data_dir
cd $data_dir

# This step implements a sophisticated parallelized parsing algorithm to extract and format detailed alignment metrics from BLAST results between PIWIL3-associated piRNAs and statistically validated cleavage sites. The custom R script processes BLAST output files in 100 parallelized batches using 16 computing cores, systematically extracting alignment coordinates, identity percentages, match/mismatch positions, and strand information for each piRNA query against corresponding cleavage site sequences. Through robust error handling and batch processing, the algorithm efficiently parses complex BLAST formats to generate comprehensive interaction tables that include both absolute genomic coordinates and relative positions within cleavage windows, ultimately consolidating all parsed results into a unified dataset for downstream motif analysis and functional validation.

# align background
# Iterate from 1 to 10
for i in {1..10}; do
    echo "Processing iteration $i"
    
    # Create directory and change to it
    cd "$data_dir"
    mkdir "blast_NPOI_filter_$i"
    cd "blast_NPOI_filter_$i"
    
    # Copy the background file
    cp "${download_dir}/bg_random_shuffle${i}.fasta" "bg_${i}.fa"
    
    # Run blastn
    blastn -query "bg_${i}.fa" \
            -num_threads 40 \
            -evalue 10000000000 \
            -db "${data_dir}${risc}_insert_db_filter_pval${pval_thres}_logfc${logfc_thres}" \
            -strand minus \
            -task blastn \
            -word_size 10 \
            -outfmt 0 \
            -ungapped > "${risc}_NPOI_BP_${i}_filter.txt"
    
    # Run the R script
    Rscript "${script_dir}Parse_NPOI.R" \
            "${risc}_NPOI_BP_${i}_filter.txt" \
            "$risc"
    
    # Combine parsed files and move result
    find . -name '*Parsed*' -print0 | \
            xargs -0 cat > \
            "${risc}_parsed_NPOI_${i}_filter.txt"
    mv "${risc}_parsed_NPOI_${i}_filter.txt" "$data_dir"
    
    # Return to data directory
    cd "$data_dir"
    
    echo "Completed iteration $i"
done

echo "All iterations completed!"

# # Manage Products
cd $script_dir
Rscript ManageProducts.R