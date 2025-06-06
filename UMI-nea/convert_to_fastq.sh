code=$(readlink -f $0)
code_dir=`dirname $code`
r1="NA"
r2="NA"
umifile="NA"
outname="out"
while getopts "hf:r:u:n:" opt; do
    case ${opt} in
        h) echo -e "usage [-h] [-f umi trimmed forward read REQUIRED]\n[-r umi trimmed reverse read used for parir-end reads]\n[-u UMI-nea output file REQUIRED]\n[-n output file prefix]"; exit 0 ;;
        f) r1=$OPTARG ;;
        r) r2=$OPTARG ;;
        u) umifile=$OPTARG ;;
        n) outname=$OPTARG ;;
        :) echo "Error: -$OPTARG requires an argument"; exit 1 ;;
        ?) echo "Error: Invalid option -${OPTARG}"; exit 1 ;;
    esac
done
if [ $r1 == "NA" ]; then echo "Please input forward reads"; exit 1; fi 
echo -e "forward read file is $r1\nreverse read file is $r2\nUMI-nea output is $umifile\noutfile is $outname"
join <(cat $r1 | paste - - - - | awk '{l=split($1,n,"_");rn=n[1];for(i=2;i<=(l-1);i++){rn=rn"_"n[i]};print n[l],rn,$2,$3,$4}' | sort -k1,1) <(cat $umifile | cut -f2,3 | sort -k1,1) | sort -k6,6 | awk -v OFS="\n" '{print $2"_"$6,$3,$4,$5}' > $outname.grouped.R1.fastq

if [ $r2 != "NA" ]; then
    join <(cat $r2 | paste - - - - | awk '{l=split($1,n,"_");rn=n[1];for(i=2;i<=(l-1);i++){rn=rn"_"n[i]};print n[l],rn,$2,$3,$4}' | sort -k1,1) <(cat $umifile | cut -f2,3 | sort -k1,1) | sort -k6,6 | awk -v OFS="\n" '{print $2"_"$6,$3,$4,$5}' > $outname.grouped.R2.fastq
fi