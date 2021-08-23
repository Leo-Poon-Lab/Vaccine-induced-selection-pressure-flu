#!/bin/bash

mkdir ../results
mkdir ../results/ivar_consensus
mkdir ../results/ivar_snvs
mkdir ../results/lofreq_snvs
mkdir ../FastQC
mkdir ../archive
mkdir ../tmp

cd ../data

data_suffix1=_1.fastq.gz
data_suffix2=_2.fastq.gz

~/softwares/samtools-1.11/bin/samtools faidx reference.fasta
~/softwares/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index reference.fasta
~/softwares/bwa/bwa index reference.fasta
java -jar ~/softwares/picard.jar CreateSequenceDictionary -R reference.fasta -O reference.dict
~/softwares/bioawk/bioawk -c fastx '{print $name"\t1\t"length($seq)"\t"$name}' reference.fasta > reference.bed

bedtools makewindows -w 3408 -b reference.bed > ref_slide_window.bed

while true
do			

	for fwdread in `ls | grep $data_suffix1 | head -n 12`
	do
		{
			sample=$(echo $fwdread | cut -d"_" -f 1)
			rwsread=$sample$data_suffix2
			# make sure file1 and file2 exist and are not opening
			if [[ -f "$fwdread" ]] && [[ -f "$rwsread" ]] && ! [[ `find $fwdread -mmin -1` ]] && ! [[ `find $rwsread -mmin -1` ]]; then
			
				# begin ngs
				# de novo assembly	
				# ~/softwares/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -t 46 -1 $fwdread -2 $rwsread -o ../results/MEGAHIT/$sample

				# Align reads to reference
				~/softwares/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 4 reference.fasta $fwdread $rwsread | ~/softwares/samtools-1.11/bin/samtools view -b --threads 4 - | ~/softwares/samtools-1.11/bin/samtools sort -o $sample"_sorted.bam" --threads 4 -m 4G - 
				~/softwares/samtools-1.11/bin/samtools index $sample"_sorted.bam"

				# Primer trimming and quality trimming
				ivar trim -q 20 -m 30 -e -b ./primers.bed -i $sample"_sorted.bam" -p $sample".trim"
				~/softwares/samtools-1.11/bin/samtools sort $sample".trim.bam" -o $sample"-trim.sorted.bam"
				~/softwares/samtools-1.11/bin/samtools index $sample"-trim.sorted.bam"

				# call consensus
				cat ref_slide_window.bed | awk '{print $1":"$2"-"$3}' | parallel -j 4 -k "~/softwares/samtools-1.11/bin/samtools mpileup -A -aa -B -d 0 -Q 0 --reference reference.fasta --region {} $sample"-trim.sorted.bam" " | awk '!seen[$2]++' > $sample"_mpileup.txt"

				# awk '!seen[$2]++' $sample"_mpileup_ori.txt" > $sample"_mpileup.txt"
				# rm $sample"_mpileup_ori.txt"

				cat $sample"_mpileup.txt" | ivar consensus -p $sample"_consensus" -n N -m 10 -q 20
				mv $sample"_consensus.fa" ../results/ivar_consensus/$sample"_consensus.fa"
				
				## readcount
				cat $sample"_mpileup.txt" | ~/softwares/mpileup2readcounts/build/mpileup2readcounts > ../FastQC/$sample"_bam.readcount.mpileup.txt"	

				## FastQC
				fastqc $sample"-trim.sorted.bam" -t 4 -o ../FastQC/					

				# call variants
				cat $sample"_mpileup.txt" |  ivar variants -p $sample -q 20 -t 0.01 -m 100 -r reference.fasta

				# ivar getmasked -i $sample".tsv" -b primers.bed -f df_pair_info.tsv -p $sample".masked.txt"

				# ivar removereads -i $sample"-trim.sorted.bam" -p $sample"-us.trimmed.masked.bam" -t $sample".masked.txt" -b primers.bed 
				# ~/softwares/samtools-1.11/bin/samtools sort $sample"-us.trimmed.masked.bam" -o $sample"-trimmed.masked.bam"
				# ~/softwares/samtools-1.11/bin/samtools index $sample"-trimmed.masked.bam"

				# ~/softwares/samtools-1.11/bin/samtools mpileup -A -aa -B -d 0 -Q 0 --reference reference.fasta $sample"-trim.sorted.bam" | ivar variants -p $sample"_ivar" -q 20 -t 0.01 -m 100 -r reference.fasta 

				# mv $sample"_ivar.tsv" ../results/ivar_snvs/"ivar_"$sample.tsv
				mv $sample".tsv" ../results/ivar_snvs/"ivar_"$sample".tsv"

				# /home/hggu/softwares/lofreq/src/lofreq/lofreq indelqual --dindel --ref reference.fasta $sample"-trim.sorted.bam" | /home/hggu/softwares/lofreq/src/lofreq/lofreq call --call-indels -f reference.fasta --force-overwrite -o ../results/lofreq_snvs/"lofreq_"$sample".vcf" -

				# mv data to achive
				touch $sample"-trim.sorted.bam"
				mv -f $sample"-trim.sorted.bam" ../tmp/$sample"-trim.sorted.bam"
				touch $sample"-trim.sorted.bam.bai"
				mv -f $sample"-trim.sorted.bam.bai" ../tmp/$sample"-trim.sorted.bam.bai"
				mv -f $fwdread ../archive/$fwdread
				mv -f $rwsread ../archive/$rwsread
				rm $sample"_consensus"*
				rm $sample"_sorted"*
				rm $sample".trim"*
				# rm $sample"-us.trimmed.masked.bam"
				rm $sample"_mpileup.txt" 
				# # rm $sample".tsv" 
				# rm $sample".masked.txt" 
				# rm $sample"-primertrim"* 
			else
				if [[ -f "$fwdread" ]]; then 
					echo $fwdread exist
					if [[ `find $fwdread -mmin -1` ]]; then 
						echo but $fwdread is still busy
					fi
				fi
				if [[ -f "$rwsread" ]]; then 
					echo $rwsread exist
					if [[ `find $rwsread -mmin -1` ]]; then 
						echo but $rwsread is still busy
					fi
				fi
			fi
		} &
	done
	echo "wait for the batch:"
	echo `ls | grep $data_suffix1 | head -n 10`
	wait
	echo "success!"
	
	if `! ls | grep "fastq.gz"`; then
		echo No new data
	fi

	sleep 6

done