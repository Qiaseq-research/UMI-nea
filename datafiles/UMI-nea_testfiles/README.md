## Test files for UMI-nea

### Extract UMI output

docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser \
qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/UMI-nea/UMI-nea_helper.sh \
-f sim1.R1.fastq.gz -r sim1.R2.fastq.gz -a 1:18 -n sim1"

#### input files

`sim1.R1.fastq.gz`: paired-end reads, only R1 contains 18bp UMI sequences
`sim1.R2.fastq.gz`: paired-end reads

#### output files

`sim1.input`: A tab separated file with number of read for each unique UMI sequence, UMI-nea input
`sim1.umi.R1.fastq`: UMI extracted R1 reads with UMI sequences attached to the end of readname with "_"
`sim1.umi.R2.fastq`: R2 reads with UMI sequences attached to the end of readname with "_"

### Run UMI-nea clustering

docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser \
qiaseqresearch/umi-nea:latest bash -c "/Download/UMI-nea/UMI-nea/UMI-nea \
-i sim1.input -l 18 -e 0.01 -t 48 -o sim1.t48.clustered"

#### input files

`sim1.input`: generated from helper script

#### output files

`sim1.t48.clustered`: A tab separated file with error corrected founder for each UMI sequences
`sim1.t48.clustered.estimate`: Model estimated number of molecules from UMI clustering reuslt
`sim1.t48.clustered.umi.reads.count`: A tab separated file with cluster size of each error corrected founder

### Generate output fastq

docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser \
qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/UMI-nea/convert_to_fastq.sh \
-f sim1.umi.R1.fastq -r sim1.umi.R2.fastq -u sim1.t48.clustered -n sim1"

#### input files

`sim1.umi.R1.fastq`: UMI extracted R1 reads generated from helper script
`sim1.umi.R2.fastq`: R2 reads generated from helper script
`sim1.t48.clustered`: output file from UMI-nea

#### output files

`sim1.grouped.R1.fastq`: R1 reads with error corrected UMI sequences attached to the end of readname with "_"
`sim1.grouped.R2.fastq`: R2 reads with error corrected UMI sequences attached to the end of readname with "_"
