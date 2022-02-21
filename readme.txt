
#==========================================================================================================================
# General information

# Brief introduction:
# These in-house scripts were wrote for RNA-seq analysis of adherent cells from Salmon distal intestine and head kidney. It includes two parts, part1 "salmonDI_mRNA_part1.sh" for trimming adapters, aligning to genome, controlling the quality, and quantifying read counts; and part2 "salmonDI_mRNA_part2.r" for comparing differential transcripts and generating plots. For more details please see below.

# Author:
# Qirui Zhang (qirui.zhang@med.lu.se)

# Citation:
# Y Park, Q Zhang, J M.O. Fernandes, and V Kiron. Phagocytic activity and marker gene expression of the adherent cells from intestine and head kidney of Atlantic salmon. XXX, 2020, XXX


#==========================================================================================================================
# Experimental design

# Leukocytes were isolated from Atlantic salmon distal intestine (DI) and head kidney (HK) tissues, and were shortly cultured to obtain adherent cells. Six fish were used for generating DI and HK samples, resulting in 6 DI replicates and 6 HK replicates. Adherent cells were harvested and subjected to extract total RNA, construct libraries, and run RNA-seq on an Illumina NextSeq 500 sequencing platform. Gene expression profiles were compared between DI and HK replicates. More details can be found in "Y Park, Q Zhang, J M.O. Fernandes, and V Kiron. XXX, 2020, XXX"


#==========================================================================================================================
# Running environment

# The script "salmonDI_mRNA_part1.sh" was originally run in author's private work directory "/home2/park/salmonDI_mrna", raw fastq files were stored in directory "/home2/park/rawData/flowcell2_rna/fastqFiles", STAR index was build in directory "/home2/park/genome/star_index". Other users have to modify these THREE directory variables, "WorkDir", "RawReads", and "Index", in line5-7 of script to their own directory for custom use. Adapter and barcode sequence are stored in file "config.txt", other users also have to modify this file to their own adapter and barcode sequence.

# The script "salmonDI_mRNA_part2.r" can be run directly in other users' own directory, please keep the input files ("salmonDI_mRNA_readcount.mx" and "sample.info") in the same directory. The file "salmonDI_mRNA_readcount.mx" has to be firstly generated from script "salmonDI_mRNA_part1.sh".


#==========================================================================================================================
# Part1 script: "salmonDI_mRNA_part1.sh"


############## Usage of the script1 ################

# Usage: bash salmonDI_mRNA_part1.sh -c config.txt -s3,4,5
#
#	Parameters
#	-c|--config	(required)	Configure file ("config.txt")
#	-s|--step	(optional)	Which step to run (see "Steps" and "Examples" below, by default run all steps)
#	-h|--help	(optional)	Return help information
#
#	Steps
#	1	Quality control of raw fastq
#	2	Trim adapter
#	3	Quality control of trimmed fastq
#	4	Align to salmon genome (ICSASG_v2)
#	5	Statistics of alignments
#	6	Generate readcounts matrix
#
#	Examples
#	bash salmonDI_mRNA_part1.sh --config=config.txt --step=3,4,5 (or -c config.txt -s3,4,5)
#	bash salmonDI_mRNA_part1.sh -c config.txt (if not mentioning steps, run all steps by default)
#	bash salmonDI_mRNA_part1.sh -h (for help information)
#
#	Function
#	Analyze salmon distal intestine and head kidney RNA-seq data, including quality control of raw fastq; trim adapter; quality control of trimmed fastq; align to genome; and generate read counts matrix.


########### Input and output of the script1 #############

# Input:
# File "config.txt" is required. It is a manually-prepared configure file and should be kept in the same directory with the script.
# Step "-s/--step" is optional, users can specify step numbers to run. If ignored, all steps will be run automatically by default.

# Output:
# Folders of qualityControl/, alignment/, cleanReads/, qualityControl/rawReads/, qualityControl/cleanReads/, and qualityControl/alignment/ are created under ${WorkDir}
# Clean reads are stored in ${WorkDir}/cleanReads/
# Alignments are stored in ${WorkDir}/alignment/
# Quality control results of raw reads, clean reads, and alignments are stored in the corresponding subdirectories under ${WorkDir}/qualityControl/.
# "salmonDI_readcounts.matrix" is the newly generated raw read counts matrix of all 12 samples (6 DI and 6 HK replicates), which is to be used in the script "salmonDI_mRNA_part2.r".


#==========================================================================================================================
# Part2 script: "salmonDI_mRNA_part2.r"


############## Usage of the script2 ################

# Usage: Rscript salmonDI_mRNA_part2.r <readCounts.mx> <sampleInfo> <readCounts/gene> <nonZeroLib/gene> <controlGroup> <adjPvalue> <foldChange> <ssaGeneName>
#
#	Parameters (all required)
#	<readCounts.mx>		Read counts matrix file ("salmonDI_mRNA_readcount.mx", generated from "salmonDI_mRNA_part1.sh")
#	<sampleInfo>		Sample information file ("sample.info")
#	<readCounts/gene>	Minimal total read counts in all 12 libraries per gene, genes with read counts less than this number will be filtered out (suggest: 100)
#	<nonZeroLib/gene>	Minimal non-zero libraries per gene, genes with non-zero libraries less than this number will be filtered out (suggest: 5)
#	<controlGroup>		Control group, "DI" or "HK" (DI: distal intestine, HK: head kidney)
#	<adjPvalue>			BH adjusted p value for running DESeq2 (suggest: 0.05)
#	<foldChange>		Minimal fold change threshold for determining DEGs (suggest: 2)
#	<ssaGeneName>		"ssa_geneName.list"
#
#	Example
#	Rscript salmonDI_mRNA_part2.r salmonDI_mRNA_readcount.mx sample.info 100 5 DI 0.05 2 ssa_geneName.list
#
#	Function
#	Run DESeq2, to find differentially expressed genes and generate plots.



########### Input and output of the script2 #############

# Input:
# "salmonDI_mRNA_readcount.mx" is the RAW read counts matrix of mRNAs generated from "salmonDI_mRNA_part1.sh".
# "sample.info" is a manually-prepared sample information file.

# Output:
# "salmonDI_mRNA_baseMean.tsv" is the DESeq2-normalized read counts of all 12 samples (6 DI and 6 HK replicates).
# "salmonDI_mRNA_DEGlist.tsv" is the full list of DEGs with an absolute fold change > foldChange, and BH adjusted < adjPvalue.
# "salmonDI_mRNA_DEGplots.pdf" is the plot collection file, including sample correlation heatmap, PCA, MA plot, normalization dispersion estimate plot, DEG volcano plot, and DEG heatmap.


#==========================================================================================================================

