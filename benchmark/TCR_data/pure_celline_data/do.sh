sample=$1

mkdir -p $sample/log

cp ~/support/Song/TCR_V3/TCR_3_cellline_data/$sample/umi_clustering.consensus.input $sample/
cp ~/support/Song/TCR_V3/TCR_3_cellline_data/$sample/$sample.primer ~/support/Song/TCR_V3/TCR_3_cellline_data/$sample/$sample.UMI $sample/
cp ~/support/Song/TCR_V3/TCR_3_cellline_data/$sample/umi_clustering.input.TR* $sample/

for g in "A" "B" "D" "G"; do
    if [ -f $sample/umi.clustered.TR$g ]; then
        continue
    fi
    maxlen=`cat $sample/umi_clustering.input.TR$g | cut -f2 | awk '{print length($1)}' | sort -nr | head -1`
    ~/code/UMI-nea/UMI-nea/UMI-nea -i $sample/umi_clustering.input.TR$g -o $sample/umi.clustered.TR$g -l $maxlen -t 64 -m 2 -a > $sample/log/UMI_clustering.TR${g}.log
done

cat allfiles.tsv | grep $sample > runList.tsv
docker run --name tcr_v3 -v ${PWD}:/home/qiauser -w /home/qiauser -v /mnt/fdkbio10/home/zhangj/code/TCR_V3/:/srv/qgen/code/ -v /mnt/fdknas01/users/JixinDeng/TCR/Old_Three_Cell_Line_Data/:/data tcr:v3.1 bash -c "python /srv/qgen/code/core/run.py runList.tsv"
docker rm tcr_v3
