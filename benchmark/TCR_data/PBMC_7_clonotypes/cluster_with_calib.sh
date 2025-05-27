rpu_cutoff=$1
code=$(readlink -f $0)
code_dir=`dirname $code`
gene_lst="A B D G"
read -ra genes <<< "$gene_lst"
umi_len=18
seq="CTGCTATACAGACTTGTCACTGGATTTAGAGTCTCTCAGCTGGTACACGGCAGGGTCAGGGTTCTGGATATTTGCTCTTACAGTTACTGTGGTTCCGGCTCCAAAGCTGAGCTTGTAGTCGTTAGAACCCTCTCTCATTGCACAGAAATACATTGCTGAGTCCCCCAGTTGTGAAGCGGAGATGACAAGGTTGGCGGATTTTCTTGCCTTCTGGAAATTCAATGAGTAGCGACCTTCTGTTGCATTTTGCTCGTCATAAGA"
qlt="BBBBBFFFFFFFGGGGGGFFFGHHHHHHHGHDEHBGHHHHHHHCGGHHGGGEEEG2GFHHA2FGHHHEFHHGFHHHGHDFHGEH5DFHB@@E2FDCEEEFGCHHFFHGHHGHHFHHGFHHGEFAEHH3?/?FHFGFHHHFGD2DFF2FHHH2GHGHHHFHFGFHFFGGHFHFD11FFDCG?<GGHGHEHHB<G?BA@C-.F0BF99C/0B9CFBEFFGFGG00BFEFBFG?--BAFFFFFFF/9BFFFFFFFFFFFFEFFF"

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
        cat $sample/umi_clustering.consensus.input | awk -v s=$sample -v ul=$umi_len -v seq=$seq -v qlt=$qlt -v OFS="\n" '{if(length($6)<ul){$6=$6""substr($8,1,ul-length($6))}else{$6=substr($6,1,ul)};print "@"$2"-"$6,substr(seq,1,150),"+",substr(qlt,1,150) > s"/"s"."$1".R1.fastq"; print "@"$2"-"$6,$6""substr(seq,1,150-ul),"+",substr(qlt,1,150) > s"/"s"."$1".R2.fastq"}'
        p=0
        for i in "${genes[@]}"; do
            p=`echo "$p+1" | bc`
            if [ $rpu_cutoff -eq 0 ]; then
                rpu_cutoff=`cat /data/$run/$sample/umi.clustered.TR$i.estimate | grep "rpu_cutoff" | awk '{print $2}'`
            fi
            touch $sample/umi.clustered.TR$i
            if [ ! -f $sample/$sample.$p.cluster ]; then
                /Download/calib/calib -f $sample/$sample.$p.R1.fastq -r $sample/$sample.$p.R2.fastq -l1 0 -l2 $umi_len -o $sample/$sample.$p. -c 8 1> $sample/log/calib.$p.log
            fi
            if [ ! -f $sample/umi.clustered.consensus.TR$i ]; then
                join -1 2 -2 1 <(cat $sample/$sample.$p.cluster | cut -f1,4 | awk '{l=split($2,n,"-");rn=n[1];for(i=2;i<l;i++){rn=rn"-"n[i]};print rn,$1,n[l]}' | sort -k2,2) <(cat $sample/$sample.$p.cluster | cut -f1 | sort | uniq -c | awk -v r=$rpu_cutoff '$1>=r{print $2}' | sort) | awk -v c=-1 '{if($1!=c){u=$3;c=$1};print $2,u}' > $sample/$sample.$p.rpu$rpu_cutoff.cluster.reads
                join -1 2 -2 1 <(cat $sample/umi_clustering.consensus.input | awk -v p=$p '$1==p' | cut -f1-5,8- | sort -k2,2) <(cat $sample/$sample.$p.rpu$rpu_cutoff.cluster.reads | sed 's/@//' | sort -k1,1) | awk -v OFS="\t" '{print $2,$1,$3,$4,$5,$8,$2,$6,$7}' | sort -k6,6 > $sample/umi.clustered.consensus.TR$i
            fi
            rm -rf $sample/$sample.$p.*.fastq
        done
    done
    cd ../
done