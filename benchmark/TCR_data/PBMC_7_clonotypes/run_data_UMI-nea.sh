for sample in `cat runList.tsv | cut -f1`; do
    mkdir -p $sample/log
    cp ~/support/Song/TCR_V3/PBMC_reproducibility/230117_VH01211_6_AAC7YLMM5/$sample/umi_clustering.consensus.input $sample/
    cp ~/support/Song/TCR_V3/PBMC_reproducibility/230117_VH01211_6_AAC7YLMM5/$sample/$sample.primer $sample/
    cp ~/support/Song/TCR_V3/PBMC_reproducibility/230117_VH01211_6_AAC7YLMM5/$sample/$sample.UMI $sample/
    cp ~/support/Song/TCR_V3/PBMC_reproducibility/230117_VH01211_6_AAC7YLMM5/$sample/umi_clustering.input.TR* $sample/

    for g in "A" "B" "D" "G"; do
        if [ -f $sample/umi.clustered.TR$g ]; then
            continue
        fi
        maxlen=`cat $sample/umi_clustering.input.TR$g | cut -f2 | awk '{print length($1)}' | sort -nr | head -1`
        ~/code/UMI-nea/UMI-nea/UMI-nea -i $sample/umi_clustering.input.TR$g -o $sample/umi.clustered.TR$g -l $maxlen -t 64 -m 2 -a > $sample/log/UMI_clustering.TR${g}.log
    done
done
docker run --name tcr_v3 -v ${PWD}:/home/qiauser -w /home/qiauser -v /mnt/fdkbio10/home/zhangj/code/TCR_V3/:/srv/qgen/code/ -v /mnt/fdknas01/raw-files/230117_VH01211_6_AAC7YLMM5/Analysis/1/Data/fastq/:/data tcr:v3.1 bash -c "python /srv/qgen/code/core/run.py runList.tsv"
docker rm tcr_v3

for sample in `cat runList.tsv | cut -f1`; do
    mkdir -p $sample/log
    cp ~/support/Song/TCR_V3/PBMC_reproducibility/M06463_0121new/$sample/umi_clustering.consensus.input $sample/
    cp ~/support/Song/TCR_V3/PBMC_reproducibility/M06463_0121new/$sample/$sample.primer $sample/
    cp ~/support/Song/TCR_V3/PBMC_reproducibility/M06463_0121new/$sample/$sample.UMI $sample/
    cp ~/support/Song/TCR_V3/PBMC_reproducibility/M06463_0121new/$sample/umi_clustering.input.TR* $sample/

    for g in "A" "B" "D" "G"; do
        if [ -f $sample/umi.clustered.TR$g ]; then
            continue
        fi
        maxlen=`cat $sample/umi_clustering.input.TR$g | cut -f2 | awk '{print length($1)}' | sort -nr | head -1`
        ~/code/UMI-nea/UMI-nea/UMI-nea -i $sample/umi_clustering.input.TR$g -o $sample/umi.clustered.TR$g -l $maxlen -t 64 -m 2 -a > $sample/log/UMI_clustering.TR${g}.log
    done
done
docker run --name tcr_v3 -v ${PWD}:/home/qiauser -w /home/qiauser -v /mnt/fdkbio10/home/zhangj/code/TCR_V3/:/srv/qgen/code/ -v /mnt/fdknas01/seq-reads/M06463_0121new/2022_M6463_121-376880952/FASTQ_Generation_2022-12-23_18_04_22Z-642543903/:/data tcr:v3.1 bash -c "python /srv/qgen/code/core/run.py runList.tsv"
docker rm tcr_v3
