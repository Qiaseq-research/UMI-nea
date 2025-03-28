for f in "10000" "100000" "200000"; do
    for l in "12" "24" "48"; do
        for a in "0.999" "0.995" "0.97"; do
            bash /Download/UMI-nea/benchmark/compare_3_tools.sh $l $a $f 20 0 0
        done
    done
done
