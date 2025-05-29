for f in "1000" "10000"; do
    for e in "0.001" "0.005"; do
        for l in "12" "18"; do
            bash /Download/UMI-nea/benchmark/compare_tools.sh $l $e $f 100 3 1 1-1-40 0 0 48 2 UMI-nea,umi-tools,UMIC-seq,calib
        done
    done
done
for f in "1000" "10000"; do
    for e in "0.01" "0.03"; do
        for l in "25" "50"; do
            bash /Download/UMI-nea/benchmark/compare_tools.sh $l $e $f 100 3 1 1-1-1 0 0 48 2 UMI-nea,umi-tools,UMIC-seq,calib
        done
    done
done
