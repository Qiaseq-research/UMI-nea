# UMI-nea

UMI-nea is a alignment-free UMI deduplication tool designed to accurately and robustly quantify the absolute count of molecules in UMI-tagged sequencing libraries. 

UMI-nea relies solely on UMI sequences and is optimized for computational efficiency, making it particularly suitable for reference free genomic applications such as TCR/BCR sequencing and large-scale UMI datasets.

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
docker run --name <umi-nea-container-name> qiaseqresearch/umi-nea:latest /Download/UMI-nea/UMI-nea/UMI-nea
```

## Run UMI-nea

There are two main usage of UMI-nea, one is do UMI clustering and quantification, the other is do quantification only based on input file

### Clustering

To run UMI-nea clustering with quantification:
```bash
UMI-nea -i <input-file> -l <max-length> -o <output-file>
```

#### Clustering parameters

##### required

1. `--maxL -l <int|default:-1>`: Max length of UMI, longer UMIs will be trimmed to the length set! Required to set for UMI clustering
2. `--in -i <fname>`: Input file, a tab delim file with three cols: GroupID, UMI, NumofReads, sorted in descending by NumofReads
3. `--out -o <fname>`: Output file, output has three cols, GroupID, UMI, FounderUMI

##### optional

4. `--minC -m <int|default:2>`: Min levenshtein distance cutoff, overruled by -e
5. `--errorR -e <float|default:none>`: If set, -m will be calculated based on binomial CI, override -m
6. `--minF -f <int|default:1>`: Min reads count for a UMI to be founder, overruled by -n or -k
7. `--nb -n <default:false>`: Apply negative binomial model to decide min reads for a founder, override -f
8. `--kp -k <default:false>`: Apply knee plot to decide min reads for a founder, override -f
9. `--auto -a <default:false>`: Combine knee plot strategy and negative binomial model to decide min reads for a founder, override -f
10. `--just -j <default:false>`: Just estimate molecule number and rpu cutoff
11. `--prob -q <default:0.001>`: probability for nb lower tail quantile cutoff in quantification
12. `--first -d <default:false>`: First founder mode, first founder below cutoff will be selected once found, which speed up computation but affect reprouciblity. Default is false which enforce to find the best founder
13. `--thread -t <int|default:10>`: Num of thread to use, minimal 2
14. `--pool -p <int|default:1000>`: Total UMIs to process in each thread at one time
15. `--help -h`: Show help

#### Clustering input

`<fname>.input`: A tab separated file with number of read for each unique UMI sequence

| Amplicon_ID | Unique UMI sequence | number of reads |
|:-----------:|:-------------------:|:---------------:|

#### Clustering output

`<fname>`: A tab separated file with error corrected founder for each UMI sequences with below columns

| Amplicon_ID | original UMI sequence | error corrected UMI sequence |
|:-----------:|:---------------------:|:----------------------------:|

`<fname>.estimate`: Model estimated number of molecules from UMI clustering reuslt

`<fname>.umi.reads.count`: A tab separated file with cluster size of each error corrected founder

| Amplicon_ID | error corrected UMI sequence | number of reads |
|:-----------:|:----------------------------:|:---------------:|

#### Example run

```bash
docker run --name umi_nea qiaseqresearch/umi-nea:latest /Download/UMI-nea/UMI-nea/UMI-nea -i sim1.input -o sim1.clustered -l 19 -e 0.001
```
```

********************Input Parameters:************************
inFile=../test/pN1000_oN100/sim_1000_100_ul18_err0.001/UMI-nea/sim1.input
outFile=sim1.clustered
maxlenUMI = 19
errorRate = 0.001
maxdist = 1
auto-Estimate = ON
threads = 10
poolSize = 1000
********************************************

All done!
```


### Quantification only

To run UMI-nea quantification only:
```bash
UMI-nea -i <input-file> -j
```

#### Example run

```bash
docker run --name umi_nea qiaseqresearch/umi-nea:latest /Download/UMI-nea/UMI-nea/UMI-nea -i sim1.input -j
```
Output

```
KP_estimate     ON
knee_angle      113
median_rpu      101
rpu_cutoff      2
estimated_molecules     1091
after_rpu-cutoff_molecules      1022

```