
mkdir -p ../featureCounts/
for i in $(ls *Aligned.sortedByCoord.out.bam | sed 's/Aligned.sortedByCoord.out.bam//' | uniq)
do
featureCounts -T 56 -p -a /home/rebee/references/mouse_2021/gencode.vM27.annotation.gtf -o ../featureCounts/${i}.txt ${i}Aligned.sortedByCoord.out.bam
done
