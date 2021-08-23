#!/bin/bash
cd ../tmp/

cp ../data/reference.fasta ./
~/softwares/samtools-1.11/bin/samtools faidx reference.fasta
~/softwares/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index reference.fasta
~/softwares/bwa/bwa index reference.fasta
java -jar ~/softwares/picard.jar CreateSequenceDictionary -R reference.fasta -O reference.dict
~/softwares/bioawk/bioawk -c fastx '{print $name"\t1\t"length($seq)"\t"$name}' reference.fasta > reference.bed

while true
do	

	for bai in `ls | grep 'trim.sorted.bam.bai' | head -n 24`
		do
			{
				sample=$(echo $bai | awk -F '-trim' '{print $1}')
				bamfile=$sample"-trim.sorted.bam"
				# make sure file1 and file2 exist and are not opening
				if [[ -f "$bamfile" ]] && ! [[ `find $bamfile -mmin -1` ]]; then

					## lofreq
					/home/hggu/softwares/lofreq/src/lofreq/lofreq indelqual --dindel --ref reference.fasta $bamfile | /home/hggu/softwares/lofreq/src/lofreq/lofreq call --call-indels -f reference.fasta --force-overwrite -o ../results/lofreq_snvs/"lofreq_"$sample".vcf" -

					mv -f $bamfile ../archive/$bamfile
					mv -f $bai ../archive/$bai

				else
					if [[ -f "$bamfile" ]]; then 
						echo $bamfile exist
						if [[ `find $bamfile -mmin -1` ]]; then 
							echo but $bamfile is still busy
						fi
					fi
					
				fi

			} &
		done
		wait
	echo "success!"

	sleep 6

done
