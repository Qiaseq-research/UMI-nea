To start a detached docker container running evaluation, see example
docker run --name umi_nea -d -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/benchmark/compare_3_tools.sh 20 0.995 1000 20 0 0"

To benchmark UMI-nea runtime, see example
docker run --name umi_nea -d -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/benchmark/wrapper_runtime.sh"

To plot UMI-nea runtime linechart, see example
python /Download/UMI-nea-publication/benchmark/plot_runtime.py

To comapre the V-measure of three tools, see example
docker run --name umi_nea -d -v ${PWD}:/home/qiauser -w /home/qiauser qiaseqresearch/umi-nea:latest bash -c "bash /Download/UMI-nea/benchmark/wrapper.sh"

To plot V-measure barplot, see example
python /Download/UMI-nea-publication/benchmark/plot_toolcompare.py