code=$(readlink -f $0)
code_dir=`dirname $code`
r1="NA"
l1="NA"
r2="NA"
l2="NA"
outname="out"
while getopts "hf:a:r:b:n:" opt; do
    case ${opt} in
        h) echo -e "usage [-h] [-f forward read REQUIRED]\n[-r reverse read used for parir-end reads]\n[-a 1-based umi position at forward reads e.g. 1:12]\n[-b 1-based umi position at reversed reads]\n[-n output file prefix]"; exit 0 ;;
        f) r1=$OPTARG ;;
        r) r2=$OPTARG ;;
        a) l1=$OPTARG ;;
        b) l2=$OPTARG ;;
        n) outname=$OPTARG ;;
        :) echo "Error: -$OPTARG requires an argument"; exit 1 ;;
        ?) echo "Error: Invalid option -${OPTARG}"; exit 1 ;;
    esac
done
if [ $r1 == "NA" ]; then echo "Please input forward reads"; exit 1; fi 
if [ $l1 == "NA" ] && [ $l2 == "NA" ]; then echo "Please input umi position at forward reads or reverse reads"; exit 1; fi
if [ $r2 == "NA" ] && [ $l2 != "NA" ]; then echo "Please input reverse reads"; exit 1; fi
echo -e "forward read file is $r1\nreverse read file is $r2\numi at forward read $l1\numi at reverse read $l2\noutfile is $outname"
if [ $l1 != "NA" ]; then
    l1_s=`echo $l1 | cut -d: -f1`
    l1_e=`echo $l1 | cut -d: -f2 | awk '{print $1+1}'`
else
    l1_s=0
    l1_e=0
fi

if [ $l2 != "NA" ]; then
    l2_s=`echo $l2 | cut -d: -f1`
    l2_e=`echo $l2 | cut -d: -f2 | awk '{print $1+1}'`
else
    l2_s=0
    l2_e=0
fi

if [[ $r1 == *".gz" ]]; then
    r1_name=`basename $r1 | sed 's/fq.gz//' | sed 's/.fastq.gz//'`
    zcat $r1 | paste - - - - | awk -v s=$l1_s -v e=$l1_e -v OFS="\t" '{u=substr($2,s,(e-s));print $1,u,substr($2,e,length($2)),"+",substr($4,e,length($4))}' > $outname.umi.R1.txt
else
    r1_name=`basename $r1 | sed 's/fq//' | sed 's/.fastq//'`
    cat $r1 | paste - - - - | awk -v s=$l1_s -v e=$l1_e -v OFS="\t" '{u=substr($2,s,(e-s));print $1,u,substr($2,e,length($2)),"+",substr($4,e,length($4))}' > $outname.umi.R1.txt
fi

if [ $r2 != "NA" ]; then
    if [[ $r2 == *".gz" ]]; then
        r2_name=`basename $r2 | sed 's/fq.gz//' | sed 's/.fastq.gz//'`
        zcat $r2 | paste - - - - | awk -v s=$l2_s -v e=$l2_e -v OFS="\t" '{u=substr($2,s,(e-s));print $1,u,substr($2,e,length($2)),"+",substr($4,e,length($4))}' > $outname.umi.R2.txt
    else
        r2_name=`basename $r2 | sed 's/fq//' | sed 's/.fastq//'`
        cat $r2 | paste - - - - | awk -v s=$l2_s -v e=$l2_e -v OFS="\t" '{u=substr($2,s,(e-s));print $1,u,substr($2,e,length($2)),"+",substr($4,e,length($4))}' > $outname.umi.R2.txt
    fi
    paste <(cat $outname.umi.R1.txt | awk -F"\t" '{print $2}') <(cat $outname.umi.R2.txt | awk -F"\t" '{print $2}') -d "" | sort | uniq -c | awk '{print "1\t"$2"\t"$1}' | sort -k3,3nr > $outname.input
    paste <(cat $outname.umi.R1.txt) <(cat $outname.umi.R2.txt) | awk -F"\t" -v OFS="\n" -v o=$outname '{u=$2""$7; print $1"_"u,$3,$4,$5 > o".umi.R1.fastq";print $1"_"u,$8,$9,$10 > o".umi.R2.fastq"}'
    rm -f $outname.umi.R*.txt
else
    cat $outname.umi.R1.txt | awk -F"\t" '{print $2}' | sort | uniq -c | awk '{print "1\t"$2"\t"$1}' | sort -k3,3nr > $outname.input
    cat $outname.umi.R1.txt | awk -F"\t" -v OFS="\n" '{print $1"_"$2,$3,$4,$5}' > $outname.umi.R1.fastq
    rm -f $outname.umi.R1.txt
fi
