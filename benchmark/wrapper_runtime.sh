for num_t in "12" "24" "36" "48"; do
    bash /Download/UMI-nea/benchmark/compare_tools.sh 12 0.995 10000 100 3 1 1-1-40 0 0 $num_t 2 UMI-nea,UMIC-seq,calib
    bash /Download/UMI-nea/benchmark/compare_tools.sh 50 0.99 10000 100 3 1 1-1-1 0 0 $num_t 2 UMI-nea,UMIC-seq,calib
done