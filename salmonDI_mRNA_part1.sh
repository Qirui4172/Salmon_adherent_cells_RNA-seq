#!/bin/bash


Time=`date "+%Y-%m-%d %H:%M:%S"`
WorkDir="/home2/park/salmonDI_mrna"
RawReads="/home2/park/rawData/flowcell2_rna/fastqFiles"
Index="/home2/park/genome/star_index"

#-------------------------------------------------------------------------------------------------
Args=`getopt -o c:s::h --long config:,step::,help -- "$@"`
eval set -- "$Args"

Usage()
{
	printf "\n\tUsage: bash salmonDI_mRNA_part1.sh -c config.txt -s3,4,5\n"

	printf "\n\tParameters\n"
	printf "\t-c|--config	(required)	Configure file (\"config.txt\")\n"
	printf "\t-s|--step	(optional)	Which step to run (see \"Steps\" and \"Examples\" below, by default run all steps)\n"
	printf "\t-h|--help	(optional)	Return help information\n\n"

	printf "\tSteps\n"
	printf "\t1	Quality control of raw fastq\n"
	printf "\t2	Trim adapter\n"
	printf "\t3	Quality control of trimmed fastq\n"
	printf "\t4	Align to salmon genome (ICSASG_v2)\n"
	printf "\t5	Statistics of alignments\n"
	printf "\t6	Generate readcounts matrix\n\n"

	printf "\tExamples\n"
	printf "\tbash salmonDI_mRNA_part1.sh --config=config.txt --step=3,4,5 (or -c config.txt -s3,4,5)\n"
	printf "\tbash salmonDI_mRNA_part1.sh -c config.txt (if not mentioning steps, run all steps by default)\n"
	printf "\tbash salmonDI_mRNA_part1.sh -h (for help information)\n\n"

	printf "\tFunction\n"
	printf "\tAnalyze salmon DI RNA-seq data: quality control of raw fastq; trim adapter; quality control of trimmed fastq; align to genome; generate read counts matrix.\n\n"

	printf "\tContact: Qirui Zhang (qirui.zhang@med.lu.se)\n"
	printf "\tDate: 16-11-2019\n\n"
}

#-------------------------------------------------------------------------------------------------
# Inspect work folder
#if [[ $PWD != $WorkDir ]];then printf "\n\t*** WARNING: You're attempting to run the script in the wrong folder, please run the script in the folder: \"$WorkDir/\"! and keep \"config.txt\" file in the SAME folder!!";Usage;exit;fi

# Inspect input parameters
printf "Inspecting parameters ...\n"
if ([[ $@ =~ "-s" ]] && [[ $# -ne 5 ]]) || ([[ $@ != *"-s"* ]] && [[ $# -ne 3 ]]);then
	if [[ $# -ne 1 ]] && [[ $# -ne 2 ]];then	# no parameter or "-h" parameter
		printf "\n\t*** WARNING: Incorrect parameters, please check again!\n"
	fi
	Usage; exit
fi

# Read config file and prepare steps
printf "Reading config file & checking steps ...\n"
if [[ $@ != *"-s"* ]];then step="all";fi
while [[ $# -ne 0 ]]
do
	case $1 in
		-c|--config)
		Samples=($(grep '^Samples' $2|cut -f2|sed 's/,/\t/g'))
		A1=$(grep '^A1' $2|cut -d "," -f2)
		A2=$(grep '^A2' $2|cut -d "," -f2)
		G1=$(grep '^G1' $2|cut -d "," -f2)
		G2=$(grep '^G2' $2|cut -d "," -f2)
		Barcodes=($(grep '^Adapter' $2|cut -d "," -f2,3))
		shift 2
		;;
		-s|--step)
		case $2 in
			"")
			printf "\n\t*** WARNING: Step numbers are missing! please specify step numbers after \"-s\" or \"--step=\"!\n"
			Usage
			shift
			exit
			;;
			[1-9]*)
			Steps=($(echo $2| sed 's/,/\t/g'))
			shift 2
			;;
		esac
		;;
		-h|--help)
		Usage
		shift
		exit
		;;
		--)
		shift
		break
		;;
	esac
done

mkdir -p qualityControl alignment cleanReads qualityControl/rawReads qualityControl/cleanReads qualityControl/alignment

#-------------------------------------------------------------------------------------------------
printf "\n=======================================================================================\n"
printf "Step1: Quality control of raw fastq files\n\n"

# Step1: Quality control of raw fastq files
if [[ ${Steps[@]} =~ 1 ]] || [[ $step == "all" ]];then
	printf "$Time\tStep1 starts\n"

	for sample in ${Samples[@]}
	do
		printf "Controlling quality of raw fastq of sample ${sample} ...\n"
		fastqc -o qualityControl/rawReads -q $RawReads/${sample}*fastq.gz
	done

	printf "$Time\tStep1 done!\n"
else
	printf "Step1 skipped\n"
fi

#-------------------------------------------------------------------------------------------------
printf "\n=======================================================================================\n"
printf "Step2: Trim adapter\n\n"

# Step2: Trim adapter
if [[ ${Steps[@]} =~ 2 ]] || [[ $step == "all" ]];then
	printf "$Time\tStep2 starts\n"

	for sampBarc in ${Barcodes[@]}
	do
		sample=${sampBarc/%,*/}
		printf "Trimming sample ${sample} ...\n"
		barcode=${sampBarc/#*,/}
		A=${A1}${barcode}${A2}
		A2=${barcode}${A2}
		G=${G1}${G2}

		cutadapt -a ${A} -q 20 --quality-base=33 -m 20 --trim-n -o cleanReads/${sample}.trim1.fq $RawReads/${sample}*fastq.gz
		if [[ $? -ne 0 ]];then echo "error in cutadapt 1";exit;fi
		cutadapt -a ${A1} -q 20 --quality-base=33 -m 20 --trim-n -o cleanReads/${sample}.trim2.fq cleanReads/${sample}.trim1.fq
		if [[ $? -ne 0 ]];then echo "error in cutadapt 2";exit;fi
		cutadapt -a ${A2} -q 20 --quality-base=33 -m 20 --trim-n -o cleanReads/${sample}.trim3.fq cleanReads/${sample}.trim2.fq
		if [[ $? -ne 0 ]];then echo "error in cutadapt 3";exit;fi
		cutadapt -g ${G} -q 20 --quality-base=33 -m 20 --trim-n -o cleanReads/${sample}.trim4.fq cleanReads/${sample}.trim3.fq
		if [[ $? -ne 0 ]];then echo "error in cutadapt 4";exit;fi
		cutadapt -g ${G2} -q 20 --quality-base=33 -m 20 --trim-n -o cleanReads/${sample}.trim5.fq cleanReads/${sample}.trim4.fq
		if [[ $? -ne 0 ]];then echo "error in cutadapt 5";exit;fi
		fastq_quality_filter -q 20 -p 80 -i cleanReads/${sample}.trim5.fq -o cleanReads/${sample}.clean.fq
		if [[ $? -ne 0 ]];then echo "error in cutadapt 6";exit;fi
		rm cleanReads/${sample}.trim*.fq

	done

	printf "$Time\tStep2 done!\n"
else
	printf "Step2 skipped\n"
fi

#-------------------------------------------------------------------------------------------------
printf "\n=======================================================================================\n"
printf "Step3: Quality control of trimmed fastq\n"

# Step3: Quality control of trimmed fastq
if [[ ${Steps[@]} =~ 3 ]] || [[ $step == "all" ]];then
	printf "$Time\tStep3 starts\n"

	for sample in ${Samples[@]}
	do
		printf "Controlling quality of trimmed sample ${sample} ...\n"
		fastqc -o qualityControl/cleanReads -q cleanReads/${sample}.clean.fq
	done

	printf "$Time\tStep3 done!\n"
else
	printf "Step3 skipped\n"
fi

#-------------------------------------------------------------------------------------------------
printf "\n=======================================================================================\n"
printf "Step4: Align to Atlantic salmon genome (ICSASG_v2)\n\n"

# Step4: Align to genome
if [[ ${Steps[@]} =~ 4 ]] || [[ $step == "all" ]];then
	printf "$Time\tStep4 starts\n"

	for sample in ${Samples[@]}
	do
		printf "Aligning sample ${sample} ...\n"
		mkdir -p alignment/${sample}
		STAR --runThreadN 16 --genomeDir $Index --readFilesIn cleanReads/${sample}.clean.fq --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix alignment/${sample}/
	done

	printf "$Time\tStep4 done!\n"
else
	printf "Step4 skipped\n"
fi

#-------------------------------------------------------------------------------------------------
printf "\n=======================================================================================\n"
printf "Step5: Statistics of alignments\n\n"

# Step5: Statistics of alignments
if [[ ${Steps[@]} =~ 5 ]] || [[ $step == "all" ]];then
	printf "$Time\tStep5 Starts\n"

	for sample in ${Samples[@]}
	do
		printf "Summarizing alignment statistics of sample ${sample} ...\n"
		mkdir -p alignment/${sample}
		samtools flagstat alignment/${sample}/Aligned.out.bam > qualityControl/alignment/${sample}.flagstat
	done

	printf "$Time\tStep5 done!\n"
else
	printf "Step5 skipped\n"
fi

#-------------------------------------------------------------------------------------------------
printf "\n=======================================================================================\n"
printf "Step6: Generate read counts matrix\n\n"

# Step6: Generate read counts matrix
if [[ ${Steps[@]} =~ 6 ]] || [[ $step == "all" ]];then
	printf "$Time\tStep6 starts\n"

	flag=1
	for sample in ${Samples[@]}
	do
		printf "Generating read counts of sample ${sample} ...\n"
		if [ $flag -eq 1 ];then
			cut -f1 alignment/${sample}/ReadsPerGene.out.tab > gene.id
			sed -i 's/:/:LOC/g' gene.id
			sed -i '1{x;p;x;}' gene.id
			flag=2
		fi
		cut -f4 alignment/${sample}/ReadsPerGene.out.tab > ${sample}.tmp
		sed -i "1i ${sample}" ${sample}.tmp
	done

	paste gene.id *tmp > salmonDI_readcounts.matrix
	sed -i '2,5d' salmonDI_readcounts.matrix
	sed -i 's/^\t//' salmonDI_readcounts.matrix
	sed -i '/^rna/d' salmonDI_readcounts.matrix
	rm gene.id *tmp

	printf "$Time\tStep6 done!\n"
else
	printf "Step6 skipped\n"
fi

#-------------------------------------------------------------------------------------------------
printf "\n=======================================================================================\n"
printf "$Time\tThe whole analysis done!\n\n"
printf "Congrats Oppa, your script has successfully finished!\nWish you a good day/Ha en fin dag/祝你有愉快的一天/좋은 하루 되세요 ^^\n\n"

