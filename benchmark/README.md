To start a detached docker container running UMI-nea
docker run --name umi_nea -it -v ${PWD}:/home/qiauser -w /home/qiauser -v <data>:/data jingxiaozhang/umi-nea-publication:v1.0.0 bash -c "bash /code/compare_3_tools.sh 20 0.995 1000 20 0 0"
