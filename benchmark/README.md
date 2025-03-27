To start a detached docker container running evaluation
docker run --name umi_nea -d -v ${PWD}:/home/qiauser -w /home/qiauser jingxiaozhang/umi-nea:latest bash -c "bash /Download/UMI-nea/benchmark/compare_3_tools.sh 20 0.995 1000 20 0 0"
