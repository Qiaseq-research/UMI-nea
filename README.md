# UMI-nea

UMI-nea is a reference-free UMI deduplication tool designed to accurately and robustly quantify the absolute count of molecules in UMI-tagged sequencing libraries. 

UMI-nea, which relies solely on UMI sequences and is optimized for computational efficiency, can be broadly applied to UMI-based workflows in both DNA-seq and RNA-seq


## Directories

```
├── benchmark                                  # benchmark and simulation module
│   ├── clustering_score.py
│   ├── compare_tools.sh
│   ├── plot_clonotype_UMI_count.py
│   ├── plot_PCA.py
│   ├── plot_runtime.py
│   ├── plot_V-measure.py
│   ├── README.md                              
│   ├── simulate_UMI_indel.py
│   ├── wrapper_runtime.sh
│   └── wrapper.sh
├── docker                                     # Dockerfile
│   ├── docker-compose.yaml
│   ├── Dockerfile
│   ├── network.py
│   ├── README.md
│   ├── requirements.txt
│   └── UMIC-seq.py
├── README.md
└── UMI-nea                                    # UMI-nea and helper scripts
    ├── convert_to_fastq.sh
    ├── edlib
    ├── include
    ├── LICENSE
    ├── Makefile
    ├── src
    ├── UMI-nea_helper.sh
    └── UMI-nea.main.cpp
```

## Installation

UMI-nea has two ways to install. You can install UMI-nea directly through source code or use our docker image

### From source code

prerequisite
1. c++14
2. boost library

```bash
git clone https://github.com/Qiaseq-research/UMI-nea.git
cd UMI-nea/UMI-nea
make
```

### From docker image

To pull docker image of umi-nea
```bash
docker pull qiaseqresearch/umi-nea:latest
```
To run UMI-nea in docker container:
```bash
docker run --name <umi-nea-container-name> -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest /Download/UMI-nea/UMI-nea/UMI-nea
```

## Run UMI-nea

There are two main usage of UMI-nea, one is do UMI clustering and quantification, the other is do quantification only based on input file

### Extract UMI

We provide a helper script to extract UMI sequence and generate UMI-nea input.

To extract UMI with helper script
```bash
bash UMI-nea_helper.sh -f <read1-file> -a <position> -r <read2-file> -b <position>
```

##### required
1. `-f <str>`: forward read fastq file, end with .fq/.fastq/.fq.gz/.fastq.gz
2. `-a <int:int>`: 1-based umi start and end positions at forward reads e.g. a 12bp umi at 1:12

##### optional
1. `-r <str>`: reverse read fastq file, end with .fq/.fastq/.fq.gz/.fastq.gz
2. `-b <int:int>`: 1-based umi start and end positions at reverse reads
3. `-n <str>`: output file prefix
4. `-h`: Show help

#### Example run

Single end
```bash
docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/UMI-nea/UMI-nea_helper.sh -f test.fastq -a 1:12"
```

Pair end
both R1 and R2 has UMI sequence
```bash
docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/UMI-nea/UMI-nea_helper.sh -f test.R1.fastq -a 1:8 -r test.R2.fastq -b 1:8"
```
Only R2 has UMI sequence
```bash
docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/UMI-nea/UMI-nea_helper.sh -f test.R1.fastq -r test.R2.fastq -b 1:18"
```

#### Extract UMI output

Extract UMI on pair-end reads with UMI at 1:18 on R2
1. UMI-nea input(no header)

| Amplicon_ID | Unique UMI sequence | number of reads |
|:-----------:|:-------------------:|:---------------:|
|      1      |  TTTACAAGCACCCCATCA |        146      |


2. UMI extracted fastqs with UMI sequences attached to the end of readname with "_"
R1
```
@read-1_TTTACAAGCACCCCATCA
TCACACACCCCGCGGCGGCGAGGCCTTAAATAGGGAAACGGCCTGAGGCGCGCGCGGGCCTGGAGCCGGGATCCGCCCTAGGGGCTCGGATCGCCGCGCGCTCGCCGCTCGCCCGCCAGCCCGCCCGTGGTCCGTGGCGGCGCGCTCCAC
+
ABBBAFBA@B@BGGGGGGFDGGFHF5G53AFGGFBGHH22AA0AAEGC10F?001EF2?EGHHGGFB43FE3FDGB333F@/F//0?B/E/F2??G?/0BF?/<//2?G2FHD2GH2F222?<<110//--;@F9BBFD;..;-9BB/9.
```

R2
```
@read-1_TTTACAAGCACCCCATCA
GTGGAGCGCGCCGCCACGGACCACGGGCGGGCTGGCGGGCGAGCGGCGAGCGCGCGGCGATCCGAGCCCCTAGGGCGGATCCCGGCTCCAGGCCCGCGCGCGCCTCAGGCCGTTTCCCTATTTAAGGCCTCG
+
FDGGFHF5G53AFGGFBGHH22AA0AAEGC10F?001EF2?EGHHGGFB43FE3FDGB333F@/F//0?B/E/F2??G?/0BF?/<//2?G2FHD2GH2F222?<<110//--;@F9BBFD;..;-9BB/9.
```

### Clustering

To run UMI-nea clustering with quantification:
```bash
./UMI-nea -i <input-file> -l <max-length> -o <output-file> -e <error-rate>
```

#### Clustering parameters

##### required

1. `--maxL -l <int|default:-1>`: Max length of UMI, longer UMIs will be trimmed to the length set! Required to set for UMI clustering
2. `--in -i <fname>`: Input file, a tab delim file with three cols: GroupID, UMI, NumofReads, sorted in descending by NumofReads
3. `--out -o <fname>`: Output file, output has three cols, GroupID, UMI, FounderUMI

##### optional

4. `--maxC -m <int|default:1>`: Max levenshtein distance cutoff within UMI cluster; if set -m do not set -e
5. `--errorR -e <float (0,1]|default:0.001>`: If set -e, do not set -m and maxC -m will be calculated based on binomial CI
6. `--minF -f <int|default:1>`: Min reads count for a UMI to be founder, overruled by -a or -n or -k
7. `--nb -n <default:false>`: Apply negative binomial model to decide min reads for a founder, override -f
8. `--kp -k <default:false>`: Apply knee plot to decide min reads for a founder, override -f
9. `--auto -a <default:true>`: Combine knee plot strategy and negative binomial model to decide min reads for a founder, override -f
10. `--just -j <default:false>`: Just estimate molecule number and rpu cutoff
11. `--prob -q <float|default:0.001>`: probability for nb lower tail quantile cutoff in quantification! Not reommended to change
12. `--angle -b <default:120>`: minimal angle for knee point
13. `--greedy -g <default:false>`: Greedy mode, first founder below cutoff will be selected once found, which speed up computation but affect reprouciblity. Default is false which enforce to find the best founder! Not recommended!
14. `--thread -t <int|default:10>`: Num of thread to use, minimal 2
15. `--pool -p <int|default:1000>`: Total UMIs to process in each thread at one time
16. `--help -h`: Show help

#### Clustering input

`<fname>.input`: A tab separated file with number of read for each unique UMI sequence(no header)

| Amplicon_ID | Unique UMI sequence | number of reads |
|:-----------:|:-------------------:|:---------------:|
|      1      |  TTTACAAGCACCCCATCA |        146      |

#### Clustering output

`<fname>`: A tab separated file with error corrected founder for each UMI sequences with below columns(no header)

| Amplicon_ID | original UMI sequence | error corrected UMI sequence |
|:-----------:|:---------------------:|:----------------------------:|
|      1      |  TTTACAAGCACCCCATCA   |      TTTACAAGCACCCCATCA      |

`<fname>.estimate`: Model estimated number of molecules from UMI clustering reuslt

```
NB_estimate     ON
median_rpu      103
rpu_cutoff      50
estimated_molecules     1000

```

`<fname>.umi.reads.count`: A tab separated file with cluster size of each error corrected founder(no header)

| Amplicon_ID | error corrected UMI sequence | number of reads |
|:-----------:|:----------------------------:|:---------------:|
|      1      |      TTTACAAGCACCCCATCA      |        154      |

#### Example run

```bash
docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest /Download/UMI-nea/UMI-nea/UMI-nea -i sim1.input -o sim1.clustered -l 19 -e 0.005
```
```
********************Input Parameters:************************
inFile=sim1.input
outFile=sim1.t48.clustered
maxlenUMI = 19
errorRate = 0.005
maxdist = 2
auto-Estimate = ON
min_ReadsPerUmi_founder  = 2
threads = 48
poolSize = 1000
********************************************

All done!
```

### Quantification only

To run UMI-nea quantification only:
```bash
./UMI-nea -i <input-file> -j
```

#### Example run

```bash
docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest /Download/UMI-nea/UMI-nea/UMI-nea -i sim1.input -j
```
#### Quantification Output

```
KP_estimate     ON
knee_angle      103
median_rpu      93
rpu_cutoff      2
estimated_molecules     1506

```

### Generate output fastq

We provide a helper script to convert UMI-nea output file to fastq files.

To generate fastq files
```bash
bash convert_to_fastq.sh -f <UMI-extracted read1-file> -r <UMI-extracted read2-file> -u <UMI-nea output file>
```

##### required
1. `-f <str>`: UMI sequences extracted forward read fastq file, end with .fq/.fastq
2. `-u <str>`: output file from UMI-nea with error corrected founder for each UMI sequences

##### optional
1. `-r <str>`: UMI sequences extracted reverse read fastq file, end with .fq/.fastq
2. `-n <str>`: output file prefix
3. `-h`: Show help

#### Input files

`UMI extracted fastq files`: UMI sequences extracted and attached to the end of readname with "_"
```
@read-1_TTTACAAGCACCCCATCA
GTGGAGCGCGCCGCCACGGACCACGGGCGGGCTGGCGGGCGAGCGGCGAGCGCGCGGCGATCCGAGCCCCTAGGGCGGATCCCGGCTCCAGGCCCGCGCGCGCCTCAGGCCGTTTCCCTATTTAAGGCCTCG
+
FDGGFHF5G53AFGGFBGHH22AA0AAEGC10F?001EF2?EGHHGGFB43FE3FDGB333F@/F//0?B/E/F2??G?/0BF?/<//2?G2FHD2GH2F222?<<110//--;@F9BBFD;..;-9BB/9.
```

`UMI-nea output file`: A tab separated file with error corrected founder for each UMI sequences with below columns(no header)

| Amplicon_ID | original UMI sequence | error corrected UMI sequence |
|:-----------:|:---------------------:|:----------------------------:|
|      1      |  TTTACAAGCACCCCATCA   |      TTTACAAGCAACCCATCA      |

#### Example run

Single end
```bash
docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/UMI-nea/convert_to_fastq.sh -f out.umi.R1.fastq -u out.clustered"
```

Pair end
```bash
docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/UMI-nea/UMI-nea_helper.sh -f out.umi.R1.fastq -r out.umi.R2.fastq -u out.clustered"
```
#### Output files

Pair end
R1
```
@read-100_TTTACAAGCAACCCATCA
TCACACACCCCGCGGCGGCGAGGCCTTAAATAGGGAAACGGCCTGAGGCGCGCGCGGGCCTGGAGCCGGGATCCGCCCTAGGGGCTCGGATCGCCGCGCGCTCGCCGCTCGCCCGCCAGCCCGCCCGTGGTCCGTGGCGGCGCGCTCCAC
+
ABBBAFBA@B@BGGGGGGFDGGFHF5G53AFGGFBGHH22AA0AAEGC10F?001EF2?EGHHGGFB43FE3FDGB333F@/F//0?B/E/F2??G?/0BF?/<//2?G2FHD2GH2F222?<<110//--;@F9BBFD;..;-9BB/9.
```
R2
```
@read-100_TTTACAAGCAACCCATCA
GTGGAGCGCGCCGCCACGGACCACGGGCGGGCTGGCGGGCGAGCGGCGAGCGCGCGGCGATCCGAGCCCCTAGGGCGGATCCCGGCTCCAGGCCCGCGCGCGCCTCAGGCCGTTTCCCTATTTAAGGCCTCG
+
FDGGFHF5G53AFGGFBGHH22AA0AAEGC10F?001EF2?EGHHGGFB43FE3FDGB333F@/F//0?B/E/F2??G?/0BF?/<//2?G2FHD2GH2F222?<<110//--;@F9BBFD;..;-9BB/9.
```