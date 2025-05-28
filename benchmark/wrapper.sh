for f in "1000" "10000"; do
    for a in "0.999" "0.995"; do
        for l in "12" "18"; do
            bash /Download/UMI-nea/benchmark/compare_tools.sh $l $a $f 100 3 1 1-1-40 0 0 48 2 UMI-nea,umi-tools,UMIC-seq,calib
        done
    done
done
for f in "1000" "10000"; do
    for a in "0.99" "0.97"; do
        for l in "25" "50"; do
            bash /Download/UMI-nea/benchmark/compare_tools.sh $l $a $f 100 3 1 1-1-1 0 0 48 2 UMI-nea,umi-tools,UMIC-seq,calib
        done
    done
done
