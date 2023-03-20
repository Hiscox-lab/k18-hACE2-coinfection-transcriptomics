
for i in $(cat /home/rebee/projects/cambridge-upr-project/bin/path-to-samples.txt)
do
mkdir -p ../star_out/${i}

	STAR --genomeDir /home/rebee/references/mouse/star \
	--readFilesIn "${i}/"*R1.fastq.gz "${i}/"*R2.fastq.gz \
       	--readFilesCommand zcat \
	--runThreadN 12 \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ../star_out/${i}
done
