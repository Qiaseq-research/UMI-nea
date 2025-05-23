# bash /code/compare_3_tools.sh 20 0.995 1000 20 0 1 0
umi_len=$1 # umi length
acc=$2 # accuracy
pN=$3 # number of founder
oN=$4 # number of offspring
num_rep=$5 # number of replicates
dist=$6 # use UMI-nea -e to infer dist set dist to 0 else input calculated dist
umic_threshold=$7 # threhold number for UMIC-seq clustering; set to 0 for auto finding threhold
mut_ratio=$8 # ins-del-sub
p1_ratio=$9 # var/mean ratio for number of children simulation with negative bionmial distribution
run_umi_nea_only=${10} # only run UMI-nea
code=$(readlink -f $0)
code_dir=`dirname $code`
name=sim_${pN}_${oN}_ul${umi_len}_acc${acc}
err_rate=`echo $acc | awk '{printf "%1.3f\n",(1-$1)}'`
time_lim=86400s # time limit for each tools
mkdir -p $name/log
cd $name

seq="CTGCTATACAGACTTGTCACTGGATTTAGAGTCTCTCAGCTGGTACACGGCAGGGTCAGGGTTCTGGATATTTGCTCTTACAGTTACTGTGGTTCCGGCTCCAAAGCTGAGCTTGTAGTCGTTAGAACCCTCTCTCATTGCACAGAAATACATTGCTGAGTCCCCCAGTTGTGAAGCGGAGATGACAAGGTTGGCGGATTTTCTTGCCTTCTGGAAATTCAATGAGTAGCGACCTTCTGTTGCATTTTGCTCGTCATAAGA"
qlt="BBBBBFFFFFFFGGGGGGFFFGHHHHHHHGHDEHBGHHHHHHHCGGHHGGGEEEG2GFHHA2FGHHHEFHHGFHHHGHDFHGEH5DFHB@@E2FDCEEEFGCHHFFHGHHGHHFHHGFHHGEFAEHH3?/?FHFGFHHHFGD2DFF2FHHH2GHGHHHFHFGFHFFGGHFHFD11FFDCG?<GGHGHEHHB<G?BA@C-.F0BF99C/0B9CFBEFFGFGG00BFEFBFG?--BAFFFFFFF/9BFFFFFFFFFFFFEFFF"
rdn="@M01750:63:000000000-KLHVC:1:1101:18381:1596"

simulate_umi() {
    rep=$1
    python $code_dir/simulate_UMI_indel.py sim${rep} $umi_len $acc $pN $oN 1 $mut_ratio $p1_ratio > log/simulate.sim${rep}.log
    cat sim${rep}.truth.labels | sort -k1,1 -k2,2n > sim.l && mv sim.l sim${rep}.truth.labels
}

get_clustering_score() {
    n_complete=`cat $1 | wc -l`
    n_truth=`cat $2 | wc -l`
    if [ $n_complete -eq $n_truth ]; then
        python $code_dir/clustering_score.py $2 $1 > $3
    else
        echo "$1 homogeneity_score is 0 completeness_score is 0 V is 0" > $3
        rm -f $1
    fi
}

run_UMI-nea() {
    rep=$1
    mkdir -p UMI-nea
    if [ ! -f UMI-nea/sim${rep}.labels ]; then
        cat sim${rep}.out | cut -f1 | sort | uniq -c | awk '{print "1\t"$2"\t"$1}' | sort -k3,3nr > UMI-nea/sim${rep}.input
        maxl=`cat UMI-nea/sim${rep}.input | cut -f2 | awk '{print length()}' | sort -nr | head -1`
        echo "$name $rep UMI-nea" >> UMI-nea.time
        if [ $dist -gt 0 ]; then
            echo "err is $dist"
            { time timeout ${time_lim} bash -c "/Download/UMI-nea/UMI-nea/UMI-nea -i UMI-nea/sim${rep}.input -o UMI-nea/sim${rep}.clustered -l $maxl -t 48 -m $dist -a >> log/UMI-nea.sim${rep}.log"; } 2>> UMI-nea.time
        else
            { time timeout ${time_lim} bash -c "/Download/UMI-nea/UMI-nea/UMI-nea -i UMI-nea/sim${rep}.input -o UMI-nea/sim${rep}.clustered -l $maxl -t 48 -e $err_rate -a >> log/UMI-nea.sim${rep}.log"; } 2>> UMI-nea.time
        fi
        join <(cat UMI-nea/sim${rep}.clustered | awk '{print $2,$3}' | sort -k1,1) <(cat UMI-nea/sim${rep}.input | awk '{print $2,$3}' | sort -k1,1) | sort -k2,2 | awk -v n=0 -v p="" '{if(p=="" || $2==p){p=$2;print $0,n}else{n+=1;p=$2;print $0,n}}' | sort -k1,1 | awk '{for(i=1;i<=$3;i++){print $1,$NF}}' > UMI-nea/sim${rep}.labels
        get_clustering_score UMI-nea/sim${rep}.labels sim${rep}.truth.labels UMI-nea.sim${rep}.score
    fi
}

run_umi-tools() {
    rep=$1
    mkdir -p umi-tools
    if [ ! -f umi-tools/sim${rep}.labels ]; then
        cat sim${rep}.out | cut -f1 | awk -v r="$rdn" -v s="$seq" -v q="$qlt" '{print r":c"NR"_"$1"\n"s"\n""+""\n"q}' > umi-tools/sim${rep}.fastq
        bwa mem -t 72 -k 10 -L 0 -B 1 -O 1 -T 18 $code_dir/refgenome/ref.fa umi-tools/sim${rep}.fastq 2> log/umi-tools.bwa.log | samtools view -Sb - 1> umi-tools/sim${rep}.bam
        samtools sort umi-tools/sim${rep}.bam -o umi-tools/sim${rep}.srt.bam
        samtools index umi-tools/sim${rep}.srt.bam
        echo "$name $rep umi-tools" >> umi-tools.time
    if [ $dist -eq 0 ]; then
        maxdist=`cat log/UMI-nea.sim${rep}.log | grep "maxdist" | awk '{print $NF}'`
    else
        maxdist=$dist
    fi
        { time timeout ${time_lim} bash -c "umi_tools group -I umi-tools/sim${rep}.srt.bam --edit-distance-threshold=$maxdist --group-out=umi-tools/sim${rep}.grouped.tsv --log=log/umi-tools.sim${rep}.log --method=adjacency"; } 2>> umi-tools.time
        if [ -s umi-tools/sim${rep}.grouped.tsv ]; then
            join -1 2 -2 1 -o 1.1 1.2 1.3 1.4 2.2 <(join <(cat umi-tools/sim${rep}.grouped.tsv | cut -f5,7 | sort | uniq -c | awk '{print $2,$3,$1}' | sort -k1,1) <(cat sim${rep}.out | cut -f1 | sort | uniq -c | awk '{print $2,$1,NR-1}' | sort -k1,1) | sort -k2,2) <(cat sim${rep}.out | cut -f1 | sort -u | awk '{print $1,NR-1}' | sort -k1,1) | sort -k1,1 | awk '{print $1,$3,$5}' | awk '{for(i=1;i<=$2;i++){print $1,$3}}' > umi-tools/sim${rep}.labels
        fi
        get_clustering_score umi-tools/sim${rep}.labels sim${rep}.truth.labels umi-tools.sim${rep}.score
    fi
}

find_threshold() {
    rep=$1
    cutoff=`echo "$umi_len*0.02" | bc -l`
    base_sim=`echo "$umi_len*0.5" | bc`
    i=0
    s0=0
    t0=0
    gt=0
    cat UMIC-seq/sim${rep}.clustertest | grep "Threshold" | cut -d" " -f 2,5 | sed 's/://g' | while read a; do
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

run_UMIC-seq() {
    rep=$1
    mkdir -p UMIC-seq
    if [ ! -f UMIC-seq/sim${rep}.clustertest ]; then
        cat sim${rep}.out | cut -f1 | awk -v r="$rdn" '{print r":c"NR"_"$1"\n"$1}' | sed 's/^@/>/g' > UMIC-seq/sim${rep}.fa
        cat sim${rep}.out | cut -f1 | awk -v r="$rdn" -v s="$seq" -v q="$qlt" '{print r":c"NR"_"$1"\n"s"\n""+""\n"q}' > UMIC-seq/sim${rep}.fastq
        s=20
        e=40
        if [ $umi_len -lt 20 ]; then
            s=5
            e=25
        fi
        gt=$s
        echo "$name $rep UMIC-seq" >> UMIC-seq.time
        { time timeout ${time_lim} bash -c "python /Download/UMIC-seq/UMIC-seq.py -T 48 clustertest -i UMIC-seq/sim${rep}.fa -o UMIC-seq/sim${rep} --steps $s $e 1 > log/UMIC-seq.sim${rep}.clustertest.log 2>&1"; } 2>> UMIC-seq.time
    fi
    if [ ! -f UMIC-seq/sim${rep}.labels ]; then
        if [ $umic_threshold -eq 0 ]; then
            find_threshold $rep
            gt=$?
        else
            gt=$umic_threshold
        fi
        echo "Threshold is $gt" > log/UMIC-seq.sim${rep}.log
        { time timeout ${time_lim} bash -c "python /Download/UMIC-seq/UMIC-seq.py -T 48 clusterfull -i UMIC-seq/sim${rep}.fa -o UMIC-seq/sim${rep} --reads UMIC-seq/sim${rep}.fastq --aln_thresh $gt --size_thresh 1 --stop_thresh 1 >> log/UMIC-seq.sim${rep}.log 2>&1"; } 2>> UMIC-seq.time
        for i in UMIC-seq/sim${rep}/cluster_*.fasta; do
            c=`echo $i | cut -d/ -f3 | cut -d_ -f2 | cut -d. -f1`
            cat $i | grep ">" | awk -v clst=$c '{split($1,n,"_"); print n[2],clst}' >> UMIC-seq/sim${rep}.labels
        done
        cat UMIC-seq/sim${rep}.labels | sort -k1,1 > t && mv t UMIC-seq/sim${rep}.labels
        get_clustering_score UMIC-seq/sim${rep}.labels sim${rep}.truth.labels UMIC-seq.sim${rep}.score
    fi
}

run_calib() {
    rep=$1
    mkdir -p calib/
    if [ ! -f calib/sim${rep}.labels ]; then
        seq150=${seq:0:150}
        seq150_rc=`echo $seq150 | rev | tr ACGT TGCA`
        echo "$name $rep umi-tools" >> calib.time
        cat sim${rep}.out | awk -v ul=$umi_len -v seq=$seq150 -v rdn=$rdn -v OFS="\n" '{if(length($1)<ul){ul=length($1)};print rdn":"NR"-"$1,substr($1,1,ul)substr(seq,1,length(seq)-ul),"+",substr($1,1,ul)substr(seq,1,length(seq)-ul)}' > calib/sim${rep}.R1.fastq
        cat sim${rep}.out | awk -v seq=$seq150_rc -v rdn=$rdn -v OFS="\n" '{print rdn":"NR"-"$1,seq,"+",seq}' > calib/sim${rep}.R2.fastq
        { time timeout ${time_lim} bash -c "/Download/calib/calib -f calib/sim${rep}.R1.fastq -r calib/sim${rep}.R2.fastq -l1 $umi_len -l2 0 -o calib/sim${rep}. -c 8 1> log/calib.sim${rep}.log 2>&1"; } 2>> calib.time
        cat calib/sim${rep}.cluster | cut -f1,4 | awk '{l=split($2,n,"-");print n[l],$1}' | sort -k1,1 > calib/sim${rep}.labels
        get_clustering_score calib/sim${rep}.labels sim${rep}.truth.labels calib.sim${rep}.score
    fi
}

echo "num_founder mean_children_num variance_children_num umi_len err_rate insertion:deletion:substitution replicate total_umi simulated_mean_children_num simulated_variance_children_num substitution_base indel_base simulated_insertion:deletion:substitution substitution_only_umi indel_umi uniq_umi tool dist thread runtime_in_sec dedup_umi_cluster V-measure homogeneity_score completeness_score RPU_cutoff RPU_cutoff_model estimated_molecule" > performance.txt
for rep in `seq 1 $num_rep`; do
    if [ ! -f sim${rep}.truth.labels ]; then
        simulate_umi $rep
    fi
    input_ratio=`echo $mut_ratio | sed 's/-/:/g'`
    var_oN=`echo "$oN*$p1_ratio" | bc`
    total_umi=`cat sim${rep}.out | wc -l`
    s_mu=`cat sim${rep}.truth.labels | awk '{print $2}' | sort | uniq -c | awk '{print $1}' | awk -v n=0 '{n+=$1}END{print n/NR}'`
    s_var=`cat sim${rep}.truth.labels | awk '{print $2}' | sort | uniq -c | awk '{print $1}' | awk -v a=0 -v m=$s_mu '{s=($1-m)*($1-m);a+=s}END{print int(a/(NR-1))}'`
    umi_err=`cat sim${rep}.out | awk '$3!=""' | wc -l`
    umi_sub=`cat sim${rep}.out | awk '$3!=""{gsub("[0-9]","",$3);gsub(":","",$3);split($3,n,"_");idl=0;for(i in n){m=n[i];split(m,n1,"/");if(length(n1[2])!=1){idl=1;break}};if(idl==0){print}}' | wc -l`
    umi_idl=`echo "$umi_err-$umi_sub" | bc`
    elist=`cat sim${rep}.out | awk -v s=0 -v i=0 -v d=0 '$3!=""{gsub("[0-9]","",$3);gsub(":","",$3);split($3,n,"_");for(x in n){m=n[x];split(m,n1,"/");if(length(n1[2])==1){s+=1}else if(length(n1[2])>1){i+=1}else{d+=1}}}END{print i,d,s}'`
    bp_ins=`echo $elist | awk '{print $1}'`
    bp_del=`echo $elist | awk '{print $2}'`
    bp_sub=`echo $elist | awk '{print $3}'`
    bp_idl=`echo "$bp_ins+$bp_del" | bc`
    bp_ratio=`echo "$bp_ins $bp_del $bp_sub" | awk '{a=$1;for(i=1;i<=3;i++){a=($i<a?$i:a)};printf "%3.3f:",$1/a;printf "%3.3f:",$2/a;printf "%3.3f",$3/a}'`
    uniq_umi=`cat sim${rep}.out | cut -f1 | sort | uniq | wc -l`
    
: <<'END'
END
    if [ ! -f UMI-nea.sim${rep}.score ]; then
        run_UMI-nea $rep
    fi
    if [ $run_umi_nea_only -eq 0 ]; then
        if [ ! -f umi-tools.sim${rep}.score ]; then
            run_umi-tools $rep
        fi

        if [ ! -f UMIC-seq.sim${rep}.score ]; then
            run_UMIC-seq $rep
        fi
        if [ ! -f calib.sim${rep}.score ]; then
            run_calib $rep
        fi
    fi
    for eval_t in "UMI-nea" "umi-tools" "UMIC-seq" "calib"; do
        if [ ! -f $eval_t.sim${rep}.score ]; then
            runtime_t="NA"
            score_v="NA"
            score_h="NA"
            score_c="NA"
            n_cluster="NA"
            rpu_cutoff="NA"
            rpu_model="NA"
            est_mol="NA"
            thread=48
        else
            if [ $eval_t != "UMIC-seq" ]; then
                runtime_t=`cat $eval_t.time | grep -A 2 "$name $rep" | tail -1 | awk '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1];print int(n)}'`
                maxdist=`cat log/UMI-nea.sim${rep}.log | grep "maxdist" | awk '{print $NF}'`
            else
                runtime_t=`cat $eval_t.time | grep -A 6 "$name $rep" | awk 'NR==3 || NR==7{print $2}' | awk -v n=0 '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1]}END{print int(n)}'`
                maxdist=`cat log/$eval_t.sim${rep}.log | head -1 | awk '{print $3}'`
            fi
            score_v=`cat $eval_t.sim${rep}.score | awk '{printf "%.4f\n",$NF}'`
            score_h=`cat $eval_t.sim${rep}.score | awk '{printf "%.4f\n",$4}'`
            score_c=`cat $eval_t.sim${rep}.score | awk '{printf "%.4f\n",$7}'`
            n_cluster=`cat $eval_t/sim${rep}.labels | cut -d" " -f2 | sort -u | wc -l`
            if [ $eval_t == "UMI-nea" ]; then
                rpu_cutoff=`cat $eval_t/sim${rep}.clustered.estimate | head -3 | tail -1 | awk '{print $NF}'`
                rpu_model=`cat $eval_t/sim${rep}.clustered.estimate | head -1 | awk '{print ($1~/NB/?"negbinom":"kneeplot")}'`
                est_mol=`cat $eval_t/sim${rep}.clustered.estimate | tail -1 | awk '{print $NF}'`
            else
                rpu_cutoff="NA"
                rpu_model="NA"
                est_mol=$n_cluster
            fi
            if [ $eval_t == "umi-tools" ]; then
                thread=1
            fi
            if [ $eval_t == "calib" ]; then
                thread=8
            fi
        fi
        echo "$pN $oN $var_oN $umi_len $err_rate $input_ratio $rep $total_umi $s_mu $s_var $bp_sub $bp_idl $bp_ratio $umi_sub $umi_idl $uniq_umi $eval_t $maxdist 48 $runtime_t $n_cluster $score_v $score_h $score_c $rpu_cutoff $rpu_model $est_mol" >> performance.txt
    done
done

