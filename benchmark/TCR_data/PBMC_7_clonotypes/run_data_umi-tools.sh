docker run --name clono7 -v ${PWD}:/home/qiauser -v ~/support/UMI-nea_paper/PBMC_reproducibility_UMI-nea/:/data -w /home/qiauser umi-nea:latest bash -c "bash /Download/UMI-nea/benchmark/TCR_data/PBMC_7_clonotypes/cluster_with_umi-tools.sh"
docker rm clono7
cd M06463_0121new/
docker run --name tcr_v3 -v ${PWD}:/home/qiauser -w /home/qiauser -v /mnt/fdkbio10/home/zhangj/code/TCR_V3/:/srv/qgen/code/ -v /mnt/fdknas01/seq-reads/M06463_0121new/2022_M6463_121-376880952/FASTQ_Generation_2022-12-23_18_04_22Z-642543903/:/data tcr:v3.1 bash -c "python /srv/qgen/code/core/run.py runList.tsv"
docker rm tcr_v3
cd ..
cd 230117_VH01211_6_AAC7YLMM5
docker run --name tcr_v3 -v ${PWD}:/home/qiauser -w /home/qiauser -v /mnt/fdkbio10/home/zhangj/code/TCR_V3/:/srv/qgen/code/ -v /mnt/fdknas01/raw-files/230117_VH01211_6_AAC7YLMM5/Analysis/1/Data/fastq/:/data tcr:v3.1 bash -c "python /srv/qgen/code/core/run.py runList.tsv"
docker rm tcr_v3
cd ../