reads=$1
assembly=$2
thread=$3


sample=$(echo "$reads" | sed 's/\./\t/' | cut -f 1 )

echo "$sample"

seqkit fq2fa ${reads} -o ${sample}.fa.gz
pigz -p ${thread} -d ${sample}.fa.gz
seqkit fx2tab -j ${thread} -l -n -i -H ${sample}.fa | sort -k2 -n | tail -n 1000 | cut -f 1 | seqkit grep -f /dev/stdin ${sample}.fa > ${sample}.top1000.fa
seqkit stat ${sample}.top1000.fa > ${sample}.top1000.stat

seqkit sliding -g -s 5000 -W 5000 -j 60 ${sample}.top1000.fa -o ${sample}.top1000.win5k.fa

makeblastdb -in ${assembly} -dbtype nucl -out db/$(basename ${assembly})/$(basename ${assembly}).db

mkdir blastout
blastn -query ${sample}.top1000.win5k.fa -out blastout/$(basename ${assembly}).blastout -db db/$(basename ${assembly})/$(basename ${assembly}).db -outfmt 6 -max_target_seqs 2 -max_hsps 2 -evalue 1e-5 -num_threads 48

python switch_error.py -b blastout/$(basename ${assembly}).blastout > $(basename ${assembly}).switch_error.txt

