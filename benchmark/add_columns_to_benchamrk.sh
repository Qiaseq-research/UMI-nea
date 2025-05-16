for f in `ls -1q /datadrive/UMI-nea-publication/benchmark_oN100_20250506/`; do
    mkdir -p $f
    cp /datadrive/UMI-nea-publication/benchmark_oN100_20250506/$f/performance.txt $f/
    echo "replicate uniq_umi umi_num umi_with_error_num umi_with_only_substitution_num umi_with_indel_num umi_with_error_ratio umi_with_only_substitution_ratio umi_with_indel_ratio" > $f/sim.txt
    echo "uniq_umi tool estimated_umi_num" > $f/umi.txt
    for rep in `seq 1 3`; do
        simf=/datadrive/UMI-nea-publication/benchmark_oN100_20250506/$f/sim${rep}.out
        uniq_umi=`cat $simf | cut -f1 | sort -u | wc -l`
        n_read=`cat $simf | wc -l`
        n_err=`cat $simf | awk '$3!=""' | wc -l`
        n_sub=`cat $simf | awk '$3!=""{gsub("[0-9]","",$3);gsub(":","",$3);split($3,n,"_");idl=0;for(i in n){m=n[i];split(m,n1,"/");if(length(n1[2])!=1){idl=1;break}};if(idl==0){print}}' | wc -l`
        n_idl=`echo "$n_err-$n_sub" | bc`
        r_err=`echo "$n_err/$n_read" | bc -l`
        r_sub=`echo "$n_sub/$n_read" | bc -l`
        r_idl=`echo "$n_idl/$n_read" | bc -l`
        echo "$rep $uniq_umi $n_read $n_err $n_sub $n_idl $r_err $r_sub $r_idl" >> $f/sim.txt
        for tools in "UMI-nea" "umi-tools" "UMIC-seq"; do
            labf=/datadrive/UMI-nea-publication/benchmark_oN100_20250506/$f/$tools/sim${rep}.labels
            n_mol=`cat $labf | cut -d" " -f2 | sort -u | wc -l`
            echo "$uniq_umi $tools $n_mol" >> $f/umi.txt
        done
    done
done
