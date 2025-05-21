To pull docker image of umi-nea, run:
docker pull qiaseqresearch/umi-nea:latest

To start a interactive docker container:
docker run -it  --name <umi-nea-container-name> umi-nea

To start a detached docker container running UMI-nea

docker run -d -v <data>:/home/UMI-nea/data --name <umi-nea-container-name>  umi-nea bash -c "/home/UMI-nea/UMI-nea/UMI-nea -i /home/UMI-nea/data/<input-file> -l <max-length>  -e <error_rate> -t <threads> -o  /home/UMI-nea/data/<output-file> 2> /home/UMI-nea/data/umi-nea.error.txt"

To run UMI-nea for UMI clustering and quantification, ie. for most users :
./UMI-nea -l <max-umi-len> -i <input-file> -e <error_rate> -t <threads>

To run UMI-nea only for quantification based on input file, ie. for most users :
./UMI-nea -i  <input-file> -j

Parameter for UMI-nea:

--maxL -l <int|default:-1>          Max length of UMI, longer UMIs will be trimmed to the length set! Required to set for UMI clustering!

--in -i <fname>:                    Input file, a tab delim file with three cols: GroupID, UMI, NumofReads, sorted in descending by NumofReads

--out -o <fname>:                   Output file, output has three cols, GroupID, UMI, FounderUMI

--minC -m <int|default:2>:          Min levenshtein distance cutoff, overruled by -e

--errorR -e <float|default:none>:   If set, -m will be calculated based on binomial CI, override -m

--minF -f <int|default:1>           Min reads count for a UMI to be founder, overruled by -n or -k

--nb -n <default:false>:            Apply negative binomial model to decide min reads for a founder, override -f

--kp -k <default:false>:            Apply knee plot to decide min reads for a founder, override -f

--auto -a <default:false>:          Combine knee plot strategy and negative binomial model to decide min reads for a founder, override -f

--just -j <default:false>:          Just estimate molecule number and rpu cutoff

--prob -q <default:0.001>:          probability for nb lower tail quantile cutoff in quantification

--first -d <default:false>:         First founder mode, first founder below cutoff will be selected once found, which speed up computation but affect reprouciblity. Default is false which enforce to find the best founder

--thread -t <int|default:10>:       Num of thread to use, minimal 2

--pool -p <int|default:1000>:       Total UMIs to process in each thread at one time

--help -h:                          Show help
