umi_len=$1 # umi length
err_rate=$2 # error rate
pN=$3 # number of founder
oN=$4 # number of offspring
num_rep=$5 # number of replicates
do_indel=$6 # 1 for simulate umi seq with indels 0 for no indel
mut_ratio=$7 # ratio of ins-del-sub
dist=$8 # edit distance for UMI-nea and umi-tools *
umic_threshold=$9 # threhold number for UMIC-seq clustering; set to 0 for auto finding threhold
thread=${10} # number of thread
p1_ratio=${11} # var/mean ratio for number of children simulation with negative bionmial distribution
tools_to_run=${12} # tools to run UMI-nea,umi-tools,UMIC-seq,calib
code=$(readlink -f $0)
code_dir=`dirname $code`
name=sim_${pN}_${oN}_ul${umi_len}_err${err_rate}
time_lim=86400s # time limit for each tools

echo -e "simulation\numi_len=$umi_len\nerr_rate=$err_rate\nnum_founder=$pN\nnum_children=$oN\nnum_replicate=$num_rep\ndo_indel=$do_indel\nratio=$mut_ratio\n UMI-nea&umi-tools_threshold=$dist\nUMIC-seq_threshold=$umic_threshold\n num_thread=$thread\nvar/mean_ratio=$p1_ratio\nevaluate_tools=$tools_to_run"

mkdir -p $name/log
cd $name

seq="GTGGAGCGCGCCGCCACGGACCACGGGCGGGCTGGCGGGCGAGCGGCGAGCGCGCGGCGATCCGAGCCCCTAGGGCGGATCCCGGCTCCAGGCCCGCGCGCGCCTCAGGCCGTTTCCCTATTTAAGGCCTCGCCGCCGCGGGGTGTGTGA"
aln="0,chr1,1000001,60,150M,*,0,0,GTGGAGCGCGCCGCCACGGACCACGGGCGGGCTGGCGGGCGAGCGGCGAGCGCGCGGCGATCCGAGCCCCTAGGGCGGATCCCGGCTCCAGGCCCGCGCGCGCCTCAGGCCGTTTCCCTATTTAAGGCCTCGCCGCCGCGGGGTGTGTGA,GTGGAGCGCGCCGCCACGGACCACGGGCGGGCTGGCGGGCGAGCGGCGAGCGCGCGGCGATCCGAGCCCCTAGGGCGGATCCCGGCTCCAGGCCCGCGCGCGCCTCAGGCCGTTTCCCTATTTAAGGCCTCGCCGCCGCGGGGTGTGTGA,NM:i:0,MD:Z:150,AS:i:150,XS:i:0"

simulate_umi() {
    rep=$1
    python $code_dir/simulate_UMI_indel.py sim${rep} $umi_len $err_rate $pN $oN $do_indel $mut_ratio $p1_ratio > log/simulate.sim${rep}.log
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
    td=$2
    mkdir -p UMI-nea
    if [ ! -f UMI-nea/sim${rep}.t$td.labels ]; then
        if [ ! -f UMI-nea/sim${rep}.input ]; then
            cat sim${rep}.out | cut -f1 | sort | uniq -c | awk '{print "1\t"$2"\t"$1}' | sort -k3,3nr > UMI-nea/sim${rep}.input
        fi
        maxl=`cat UMI-nea/sim${rep}.input | cut -f2 | awk '{print length()}' | sort -nr | head -1`
        echo "$name r=$rep t=$td UMI-nea" >> UMI-nea.time
        if [ $dist -eq 0 ]; then
            { time timeout ${time_lim} bash -c "/Download/UMI-nea/UMI-nea/UMI-nea -i UMI-nea/sim${rep}.input -o UMI-nea/sim${rep}.t$td.clustered -l $maxl -t $td -e $err_rate >> log/UMI-nea.sim${rep}.t$td.log"; } 2>> UMI-nea.time
            dist=`cat log/UMI-nea.sim${rep}.t$td.log | grep "maxdist" | awk '{print $NF}'`
        else
            { time timeout ${time_lim} bash -c "/Download/UMI-nea/UMI-nea/UMI-nea -i UMI-nea/sim${rep}.input -o UMI-nea/sim${rep}.t$td.clustered -l $maxl -t $td -m $dist >> log/UMI-nea.sim${rep}.t$td.log"; } 2>> UMI-nea.time
        fi
        join <(cat UMI-nea/sim${rep}.t$td.clustered | awk '{print $2,$3}' | sort -k1,1) <(cat UMI-nea/sim${rep}.input | awk '{print $2,$3}' | sort -k1,1) | sort -k2,2 | awk -v n=0 -v p="" '{if(p=="" || $2==p){p=$2;print $0,n}else{n+=1;p=$2;print $0,n}}' | sort -k1,1 | awk '{for(i=1;i<=$3;i++){print $1,$NF}}' > UMI-nea/sim${rep}.t$td.labels
        get_clustering_score UMI-nea/sim${rep}.t$td.labels sim${rep}.truth.labels UMI-nea.sim${rep}.t$td.score
    fi
}

run_umi-tools() {
    rep=$1
    td=$2
    mkdir -p umi-tools
    if [ ! -f umi-tools/sim${rep}.t$td.labels ]; then
        if [ ! -f umi-tools/sim${rep}.srt.bam.bai ]; then
            echo -e "@SQ\tSN:chr1\tLN:248956422" > umi-tools/sim${rep}.sam
            cat sim${rep}.out | cut -f1 | awk -v a=$aln '{gsub(",","\t",a);print "read-"NR"_"$1"\t"a}' >> umi-tools/sim${rep}.sam
            samtools sort umi-tools/sim${rep}.sam | samtools view - -Sb -o umi-tools/sim${rep}.srt.bam
            samtools index umi-tools/sim${rep}.srt.bam
        fi
        if [ -f log/UMI-nea.sim${rep}.*.log ]; then
            dist=`cat log/UMI-nea.sim${rep}.*.log | grep "maxdist" | head -1 | awk '{print $NF}'`
        fi
        if [ $dist -gt 0 ]; then
            echo "$name r=$rep t=$td umi-tools" >> umi-tools.time
            { time timeout ${time_lim} bash -c "umi_tools group -I umi-tools/sim${rep}.srt.bam --edit-distance-threshold=$dist --group-out=umi-tools/sim${rep}.grouped.tsv --log=log/umi-tools.sim${rep}.t$td.log"; } 2>> umi-tools.time
            if [ -s umi-tools/sim${rep}.grouped.tsv ]; then
                tail -n+2 umi-tools/sim${rep}.grouped.tsv | cut -f5,9 | awk '{print $1,$2}' | sort -k1,1 > umi-tools/sim${rep}.t$td.labels
            fi
        else
            echo "Distance value cannot be 0 !"
        fi
        get_clustering_score umi-tools/sim${rep}.t$td.labels sim${rep}.truth.labels umi-tools.sim${rep}.t$td.score
    fi
}

find_threshold() {
    rep=$1
    td=$2
    cutoff=`echo "$umi_len*0.02" | bc -l`
    base_sim=`echo "$umi_len*0.5" | bc`
    i=0
    s0=0
    t0=0
    gt=0
    cat UMIC-seq/sim${rep}.t$td.clustertest | grep "Threshold" | cut -d" " -f 2,5 | sed 's/://g' | while read a; do
        t=`echo $a | cut -d" " -f1`
        s=`echo $a | cut -d" " -f2`
        if [ $i -gt 0 ]; then
            c=`echo "$s-$s0" | bc -l`
            if (( $(echo "$c < $cutoff" | bc -l) )) && (( $(echo "$s > $base_sim" | bc -l) )) && (( $(echo "$t0 >= 15" | bc -l) )); then
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
    td=$2
    mkdir -p UMIC-seq
    if [ ! -f UMIC-seq/sim${rep}.t$td.clustertest ]; then
        cat sim${rep}.out | cut -f1 | awk '{print ">read-"NR"_"$1"\n"$1}' > UMIC-seq/sim${rep}.fa
        cat sim${rep}.out | cut -f1 | awk -v s="$seq" '{print "@read-"NR"_"$1"\n"s"\n""+""\n"s}' > UMIC-seq/sim${rep}.fastq
        s=20
        e=40
        if [ $umi_len -lt 20 ]; then
            s=5
            e=25
        fi
        gt=$s
        echo "$name r=$rep t=$td UMIC-seq" >> UMIC-seq.time
        { time timeout ${time_lim} bash -c "python /Download/UMIC-seq/UMIC-seq.py -T $td clustertest -i UMIC-seq/sim${rep}.fa -o UMIC-seq/sim${rep}.t$td --steps $s $e 1 > log/UMIC-seq.sim${rep}.t$td.clustertest.log 2>&1"; } 2>> UMIC-seq.time
    fi
    if [ ! -f UMIC-seq/sim${rep}.t$td.labels ]; then
        if [ $umic_threshold -eq 0 ]; then
            find_threshold $rep $td
            gt=$?
        else
            gt=$umic_threshold
        fi
        echo "Threshold is $gt" > log/UMIC-seq.sim${rep}.t$td.log
        { time timeout ${time_lim} bash -c "python /Download/UMIC-seq/UMIC-seq.py -T $td clusterfull -i UMIC-seq/sim${rep}.fa -o UMIC-seq/sim${rep}.t$td --reads UMIC-seq/sim${rep}.fastq --aln_thresh $gt --size_thresh 1 --stop_thresh 1 >> log/UMIC-seq.sim${rep}.t$td.log 2>&1"; } 2>> UMIC-seq.time
        for i in UMIC-seq/sim${rep}.t$td/cluster_*.fasta; do
            c=`echo $i | cut -d/ -f3 | cut -d_ -f2 | cut -d. -f1`
            cat $i | grep ">" | awk -v clst=$c '{split($1,n,"_"); print n[2],clst}' >> UMIC-seq/sim${rep}.t$td.labels
        done
        cat UMIC-seq/sim${rep}.t$td.labels | sort -k1,1 > t && mv t UMIC-seq/sim${rep}.t$td.labels
        get_clustering_score UMIC-seq/sim${rep}.t$td.labels sim${rep}.truth.labels UMIC-seq.sim${rep}.t$td.score
    fi
}

run_calib() {
    rep=$1
    td=$2
    mkdir -p calib/
    if [ ! -f calib/sim${rep}.t$td.labels ]; then
        seq_rc=`echo $seq | rev | tr ACGT TGCA`
        echo "$name r=$rep t=$td calib" >> calib.time
        if [ ! -f calib/sim${rep}.R1.fastq ] || [ ! -f calib/sim${rep}.R2.fastq ]; then
            cat sim${rep}.out | awk -v ul=$umi_len -v seq=$seq -v OFS="\n" '{if(length($1)<ul){ul=length($1)};print "@read:"NR"-"$1,substr($1,1,ul)substr(seq,1,length(seq)-ul),"+",substr($1,1,ul)substr(seq,1,length(seq)-ul)}' > calib/sim${rep}.R1.fastq
            cat sim${rep}.out | awk -v seq=$seq_rc -v OFS="\n" '{print "@read:"NR"-"$1,seq,"+",seq}' > calib/sim${rep}.R2.fastq
        fi
        { time timeout ${time_lim} bash -c "/Download/calib/calib -f calib/sim${rep}.R1.fastq -r calib/sim${rep}.R2.fastq -l1 $umi_len -l2 0 -o calib/sim${rep}.t$td. -c $td 1> log/calib.sim${rep}.t$td.log 2>&1"; } 2>> calib.time
        cat calib/sim${rep}.t$td.cluster | cut -f1,4 | awk '{l=split($2,n,"-");print n[l],$1}' | sort -k1,1 > calib/sim${rep}.t$td.labels
        get_clustering_score calib/sim${rep}.t$td.labels sim${rep}.truth.labels calib.sim${rep}.t$td.score
    fi
}

if [ ! -f performance.txt ]; then
    echo "num_founder mean_children_num variance_children_num umi_len err_rate insertion-deletion-substitution replicate total_umi simulated_mean_children_num simulated_variance_children_num substitution_base indel_base simulated_insertion-deletion-substitution substitution_only_umi indel_umi uniq_umi tool clustering_threshold thread runtime_in_sec dedup_umi_cluster V-measure homogeneity_score completeness_score RPU_cutoff RPU_cutoff_model estimated_molecule" > performance.txt
fi
for rep in `seq 1 $num_rep`; do
    if [ ! -f sim${rep}.truth.labels ]; then
        simulate_umi $rep
    fi
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
    bp_ratio=`echo "$bp_ins $bp_del $bp_sub" | awk '{a=$1;for(i=1;i<=3;i++){a=($i<a?$i:a)};printf "%3.3f-",$1/a;printf "%3.3f-",$2/a;printf "%3.3f",$3/a}'`
    uniq_umi=`cat sim${rep}.out | cut -f1 | sort | uniq | wc -l`

: <<'END'
END
    for eval_t in `echo $tools_to_run | sed 's/,/\n/g'`; do
        if [ $eval_t == "umi-tools" ]; then
            td=1
        else
            td=$thread
        fi
        scoref=$eval_t.sim${rep}.t$td.score
        if [ ! -f $scoref ]; then
            run_${eval_t} $rep $td
        fi
        runtime_t="NA"
        rpu_cutoff="NA"
        rpu_model="NA"
        score_v=`if [ -f $scoref ]; then cat $scoref | awk '{printf "%.4f\n",$NF}'; else echo "NA"; fi`
        score_h=`if [ -f $scoref ]; then cat $scoref | awk '{printf "%.4f\n",$4}'; else echo "NA"; fi`
        score_c=`if [ -f $scoref ]; then cat $scoref | awk '{printf "%.4f\n",$7}'; else echo "NA"; fi`
        n_cluster=`cat $eval_t/sim${rep}.t$td.labels | cut -d" " -f2 | sort -u | wc -l`
        est_mol=$n_cluster
        case $eval_t in
        UMIC-seq)
            runtime_t=`cat $eval_t.time | grep -A 6 "$name r=$rep t=$td" | awk 'NR==3 || NR==7{print $2}' | awk -v n=0 '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1]}END{printf "%.2f\n",n}'`
            maxdist=`cat log/$eval_t.sim${rep}.t$td.log | head -1 | awk '{print $3}'`
            ;;
        calib)
            runtime_t=`cat $eval_t.time | grep -A 2 "$name r=$rep t=$td" | tail -1 | awk '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1];printf "%.2f\n",n}'`
            maxdist=`cat log/$eval_t.sim${rep}.t$td.log | grep "error_tolerance" | head -1 | awk '{print $2}'`
            ;;
        umi-tools)
            runtime_t=`cat $eval_t.time | grep -A 2 "$name r=$rep t=$td" | tail -1 | awk '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1];printf "%.2f\n",n}'`
            maxdist=`cat log/$eval_t.sim${rep}.t$td.log | grep "# threshold" | awk '{print $4}'`
            ;;
        UMI-nea)
            runtime_t=`cat $eval_t.time | grep -A 2 "$name r=$rep t=$td" | tail -1 | awk '{split($NF,a,"m");n+=a[1]*60;split(a[2],b,"s");n+=b[1];printf "%.2f\n",n}'`
            maxdist=`cat log/$eval_t.sim${rep}.t$td.log | grep "maxdist" | awk '{print $NF}'`
            rpu_cutoff=`cat $eval_t/sim${rep}.t$td.clustered.estimate | head -3 | tail -1 | awk '{print $NF}'`
            rpu_model=`cat $eval_t/sim${rep}.t$td.clustered.estimate | head -1 | awk '{print ($1~/NB/?"negbinom":"kneeplot")}'`
            est_mol=`cat $eval_t/sim${rep}.t$td.clustered.estimate | tail -1 | awk '{print $NF}'`
            ;;
        *)
            echo "invalid tool"
        esac
        echo "$pN $oN $var_oN $umi_len $err_rate $mut_ratio $rep $total_umi $s_mu $s_var $bp_sub $bp_idl $bp_ratio $umi_sub $umi_idl $uniq_umi $eval_t $maxdist $td $runtime_t $n_cluster $score_v $score_h $score_c $rpu_cutoff $rpu_model $est_mol" >> performance.txt
    done
done
