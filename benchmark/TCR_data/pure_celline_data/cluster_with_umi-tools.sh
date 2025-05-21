code=$(readlink -f $0)
code_dir=`dirname $code`
gene_lst="A B D G"
read -ra genes <<< "$gene_lst"

for sample in `cat runList.tsv | cut -f1`; do
    mkdir -p $sample/log
    cp /data/$sample/umi_clustering.consensus.input $sample/
    cp /data/$sample/$sample.primer $sample/
    cp /data/$sample/$sample.UMI $sample/
    cat $sample/umi_clustering.consensus.input | awk -v s=$sample '{print "@"$2"_"$6"\n"$3"\n+\n"$4 > s"/"s"."$1".R1.fastq"; print "@"$2"_"$6"\n"$8"\n+\n"$9 > s"/"s"."$1".R2.fastq"}'
    p=0
    for i in "${genes[@]}"; do
        p=`echo "$p+1" | bc` 
        bwa mem -t 72 -k 10 -L 0 -B 1 -O 1 -T 18 $code_dir/refgenome/ref.fa $sample/$sample.$p.R1.fastq $sample/$sample.$p.R2.fastq 2> $sample/log/bwa.TR$i.log | samtools view -Sb - 2>> $sample/log/bwa.TR$i.log | samtools sort - -o $sample/TR$i.srt.bam 2>> $sample/log/bwa.TR$i.log
        samtools index $sample/TR$i.srt.bam
        umi_tools group -I $sample/TR$i.srt.bam --edit-distance-threshold=1 --group-out=$sample/TR$i.grouped.tsv --log=$sample/log/umi_tools.TR$i.log --method=adjacency
        touch $sample/umi.clustered.TR$i
        join -1 2 -2 1 <(cat $sample/umi_clustering.consensus.input | awk -v p=$p '$1==p' | cut -f1-5,8- | sort -k2,2) <(tail -n+2 $sample/TR$i.grouped.tsv | cut -f1,7 | sed 's/_/\t/' | cut -f1,3 | sort -k1,1) | awk -v OFS="\t" '{print $2,$1,$3,$4,$5,$8,$2,$6,$7}' | sort -k6,6 > $sample/umi.clustered.consensus.TR$i
        rm -rf $sample/$sample.$p.*.fastq
    done
done