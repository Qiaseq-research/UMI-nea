gt=$1
rpu_cutoff=$2
code=$(readlink -f $0)
code_dir=`dirname $code`
gene_lst="A B D G"
read -ra genes <<< "$gene_lst"
umi_len=18
seq="CTGCTATACAGACTTGTCACTGGATTTAGAGTCTCTCAGCTGGTACACGGCAGGGTCAGGGTTCTGGATATTTGCTCTTACAGTTACTGTGGTTCCGGCTCCAAAGCTGAGCTTGTAGTCGTTAGAACCCTCTCTCATTGCACAGAAATACATTGCTGAGTCCCCCAGTTGTGAAGCGGAGATGACAAGGTTGGCGGATTTTCTTGCCTTCTGGAAATTCAATGAGTAGCGACCTTCTGTTGCATTTTGCTCGTCATAAGA"
qlt="BBBBBFFFFFFFGGGGGGFFFGHHHHHHHGHDEHBGHHHHHHHCGGHHGGGEEEG2GFHHA2FGHHHEFHHGFHHHGHDFHGEH5DFHB@@E2FDCEEEFGCHHFFHGHHGHHFHHGFHHGEFAEHH3?/?FHFGFHHHFGD2DFF2FHHH2GHGHHHFHFGFHFFGGHFHFD11FFDCG?<GGHGHEHHB<G?BA@C-.F0BF99C/0B9CFBEFFGFGG00BFEFBFG?--BAFFFFFFF/9BFFFFFFFFFFFFEFFF"

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
        cp /data/$run/$sample/umi_clustering.consensus.input $sample/
        cp /data/$run/$sample/$sample.primer $sample/
        cp /data/$run/$sample/$sample.UMI $sample/
        cat $sample/umi_clustering.consensus.input | awk -v s=$sample -v ul=$umi_len -v seq=$seq -v qlt=$qlt -v OFS="\n" '{if(length($6)<ul){$6=$6""substr($8,1,ul-length($6))}else{$6=substr($6,1,ul)};print "@"$2"_"$6,$6""substr(seq,1,150-ul),"+",substr(qlt,1,150) > s"/"s"."$1".fastq"; print ">"$2"_"$6,$6 > s"/"s"."$1".fa"}'

        p=0
        for i in "${genes[@]}"; do
            p=`echo "$p+1" | bc` 
            python /Download/UMIC-seq/UMIC-seq.py -T 48 clustertest -i $sample/$sample.$p.fa -o $sample/$sample.$p --steps 5 25 1 > $sample/log/UMIC-seq.$p.clustertest.log 2>&1
            if [ $gt -eq 0 ]; then
                find_threshold $sample/$sample.$p.clustertest
                gt=$?
            fi
            python /Download/UMIC-seq/UMIC-seq.py -T 48 clusterfull -i $sample/$sample.$p.fa -o $sample/$sample.$p --reads $sample/$sample.$p.fastq --aln_thresh $gt --size_thresh 1 --stop_thresh 1 >> $sample/log/UMIC-seq.$p.log 2>&1
        done
    done
    cd ..
done
