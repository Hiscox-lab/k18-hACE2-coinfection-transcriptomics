
#!/bin/bash
for fn in $(ls)
do
echo "Processing sample ${fn}"
salmon quant -i /home/hlrpenri/references/mouse_gencode_salmon -l A \
         -1 "${fn}/"*R1.fastq.gz \
         -2 "${fn}/"*R2.fastq.gz \
         -p 100 --validateMappings -o ../quants/${fn} --seqBias --gcBias
done 

