gt=$1
rpu_cutoff=$2
code=$(readlink -f $0)
code_dir=`dirname $code`
gene_lst="A B D G"
read -ra genes <<< "$gene_lst"

find_threshold() {
    infile=$1
    cutoff=`echo "16*0.02" | bc -l`
    base_sim=`echo "16*0.5" | bc`
    i=0
    s0=0
    t0=0
    gt=0
    cat $infile | grep "Threshold" | cut -d" " -f 2,5 | sed 's/://g' | while read a; do
        t=`echo $a | cut -d" " -f1`
        s=`echo $a | cut -d" " -f2`
        if [ $i -gt 0 ]; then
            c=`echo "$s-$s0" | bc -l`
            if (( $(echo "$c < $cutoff" | bc -l) )) && (( $(echo "$t > $base_sim" | bc -l) )) && (( $(echo "$t > 15" | bc -l) )); then
                gt=$t0
                return $gt
            fi
        fi
        s0=$s
        t0=$t
        i=`echo "$i+1" | bc`
    done
}

for run in "M06463_0121new" "230117_VH01211_6_AAC7YLMM5"; do
    mkdir -p $run/
    cd $run
    cp /data/$run/runList.tsv .
    for sample in `cat runList.tsv | cut -f1`; do
        mkdir -p $sample/log
        if [ ! -f $sample/umi_clustering.consensus.input ]; then
            cp /data/$run/$sample/umi_clustering.consensus.input $sample/
            cp /data/$run/$sample/$sample.primer $sample/
            cp /data/$run/$sample/$sample.UMI $sample/
        fi
        cat $sample/umi_clustering.consensus.input | awk -v s=$sample '$6!~/N/{print "@"$2"_"$6"\n"$3"\n+\n"$4 > s"/"s"."$1".fastq";print ">"$2"_"$6"\n"$6 > s"/"s"."$1".fa"}'
        p=0
        for i in "${genes[@]}"; do
            p=`echo "$p+1" | bc`
            if [ $rpu_cutoff -eq 0 ]; then
                rpu_cutoff=`cat /data/$run/$sample/umi.clustered.TR$i.estimate | grep "rpu_cutoff" | awk '{print $2}'`
            fi
            python /Download/UMIC-seq/UMIC-seq.py -T 48 clustertest -i $sample/$sample.$p.fa -o $sample/$sample.$p --steps 5 25 1 > $sample/log/UMIC-seq.$p.clustertest.log 2>&1
            if [ $gt -eq 0 ]; then
                find_threshold $sample/$sample.$p.clustertest
                gt=$?
            fi
            touch $sample/umi.clustered.TR$i
            if [ ! -f $sample/$sample.$p.rpu$rpu_cutoff.cluster.reads ]; then
                python /Download/UMIC-seq/UMIC-seq.py -T 48 clusterfull -i $sample/$sample.$p.fa -o $sample/$sample.$p --reads $sample/$sample.$p.fastq --aln_thresh $gt --size_thresh 1 --stop_thresh 1 >> $sample/log/UMIC-seq.$p.log 2>&1
                for clt in `ls -1q $sample/$sample.$p/*.fasta`; do
                    rpu=`cat $clt | grep ">" | wc -l`
                    if [ $rpu -ge $rpu_cutoff ]; then
                        cat $clt | grep ">" | sed 's/>/@/' | sed 's/_/ /' | awk -v s="" '{if(NR==1){r=$2};print $1,r}' >> $sample/$sample.$p.rpu$rpu_cutoff.cluster.reads
                    fi
                done
            fi
            join -1 2 -2 1 <(cat $sample/umi_clustering.consensus.input | awk -v p=$p '$1==p' | cut -f1-5,8- | sort -k2,2) <(cat $sample/$sample.$p.rpu$rpu_cutoff.cluster.reads | sed 's/@//' | sort -k1,1) | awk -v OFS="\t" '{print $2,$1,$3,$4,$5,$8,$2,$6,$7}' | sort -k6,6 > $sample/umi.clustered.consensus.TR$i
            rm -rf $sample/$sample.$p.*.fastq
            rm -rf $sample/$sample.$p/
        done
    done
    cd ..
done
