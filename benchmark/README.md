# Benchmark with simulated datasets

We compared the performance of UMI-nea against both alignment based and alignment-free UMI clustering tools with simulated data. The alignment-free tools we tested are calib (Orabi et al., 2019) and UMIc-seq (Zurek et al., 2020). The alignment-based tool is UMI-tools (Smith et al., 2017). 

## Run simulation

UMI-nea includes a module that performed simulations of UMI sequencing under different read types and library sizes with know clustering truth. We start this module from UMIc-seq simulation approach (https://github.com/fhlab/UMIC-seq/blob/master/figures/SuppFig1/simulate_UMI-length-errorrate.py). This module is used for comparing UMI-nea performance against other tools. To use this module, we provide a docker image.

To pull docker image of umi-nea
```bash
docker pull qiaseqresearch/umi-nea:latest
```

### Generating simulation datasets

To run generate a simuation dataset using docker image
```bash
docker run --name <umi-nea-container-name> -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "python /Download/UMI-nea/benchmark/simulate_UMI_indel.py <output-file> <umi_len> <err_rate> <num_founder> <num_children> <include_indel> <mutation_ratio> <dispersion>"
```

#### Simulation parameters required

1. `output-file <str>`: Output file name
2. `umi_len <int>`: UMI length
3. `err_rate <float>`: per UMI base error rate raneg from 0-1
4. `num_founder <int>`: Number of founder UMI to simulate
5. `num_children <int>`: Number of progeny UMI to simulate
6. `include_indel <int>`: 1 for including insertion and deletion errors, 0 for substitution only
7. `mutation_ratio <int-int-int>`: Ratio of number of insertion-deletion-substitution
8. `dispersion <float>`: Ratio of variance/mean of negative bionomial distribution of number of progeny UMI

#### Simulation output

`<fname>.out`: A tab separated file with founder UMI sequence and simulated error for each UMI sequence simulated

| simulated UMI sequence | founder UMI sequence | error simulated |
|:----------------------:|:--------------------:|:---------------:|

`<fname>.truth.labels`: A space separated file with founder UMI ID for each UMI sequence simulated

| simulated UMI sequence | founder UMI ID |
|:----------------------:|:--------------:|


#### Example run

```bash
docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "python /Download/UMI-nea/benchmark/simulate_UMI_indel.py sim 18 0.005 10000 100 1 1-1-40 2"
```

### Benchmark tools

To comapre V-measure and runtime of the 4 tools, we include a benchmark script that output metrics of comparisons

To run benchmark using docker image
```bash
docker run --name <umi-nea-container-name> -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/benchmark/compare_tools.sh <umi_len> <err_rate> <num_founder> <num_children> <num_replicates> <include_indel> <mutation_ratio> <edit_distance> <umic_threshold> <thread> <dispersion> <tools_to_compare>"
``` 

#### Benchmark parameters required

1. `umi_len <int>`: UMI length
2. `err_rate <float>`: per UMI base error rate raneg from 0-1
3. `num_founder <int>`: Number of founder UMI to simulate
4. `num_children <int>`: Number of progeny UMI to simulate
5. `num_replicates <int>`: Number of replicates to simulate
6. `include_indel <int>`: 1 for including insertion and deletion errors, 0 for substitution only
7. `mutation_ratio <int-int-int>`: Ratio of number of insertion-deletion-substitution
8. `edit_distance <int>`: Edit distance for UMI-nea and umi-tools, set to 0 for UMI-nea to determine by itself
9. `umic_threshold <int>`: Alignment threshold for UMIC-seq, set to 0 for UMIC-seq to determine by itself
10. `thread <int>`: Number of thread to use, minimal 2
11. `dispersion <float>`: Ratio of variance/mean of negative bionomial distribution of number of progeny UMI

#### Benchmark output

`performance.txt`: metrics that summarize the simulation dataset and tools performance
* num_founder
* mean_children_num
* variance_children_num
* umi_len
* err_rate
* insertion-deletion-substitution: input ratio
* replicate
* total_umi
* simulated_mean_children_num
* simulated_variance_children_num
* substitution_base
* indel_base
* simulated_insertion-deletion-substitution: ratio calculated from simulation data
* substitution_only_umi: number of UMI contains only substitution errors
* indel_umi: number of UMI contains indel errors
* uniq_umi
* tool
* clustering_threshold: edit_distance for UMI-nea/umi-tools, alignment threshold for UMIc-seq, error_tolerance for calib
* thread
* runtime_in_sec
* dedup_umi_cluster: number of UMI clusters after deduplication
* V-measure
* homogeneity_score
* completeness_score
* RPU_cutoff: UMI-nea only, reads per UMI cutoff for quantification
* RPU_cutoff_model: UMI-nea only, model used for estimating reads per UMI cutoff for quantification
* estimated_molecule: UMI-nea only, quantification of number of molecules

#### Example run

```bash
docker run --name umi_nea -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/benchmark/compare_tools.sh 18 0.005 10000 100 3 1 1-1-40 0 0 48 2 UMI-nea,umi-tools,UMIC-seq,calib"
```

