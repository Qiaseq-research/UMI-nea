for num_t in "12" "24" "36" "48"; do
    bash /Download/UMI-nea/benchmark/compare_3_tools.sh 12 0.995 1000 100 3 0 0 $num_t 1-1-40 2 UMI-nea,UMIC-seq,calib
    bash /Download/UMI-nea/benchmark/compare_3_tools.sh 50 0.99 1000 100 3 0 0 $num_t 1-1-1 2 UMI-nea,UMIC-seq,calib
done