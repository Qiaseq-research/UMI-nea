code=$(readlink -f $0)
code_dir=`dirname $code`
echo -e "platform\tcondition\treplicate\tclonotype\tnumber of UMI" > umi_count_7_clonotypes.txt
for clono in `cat $code_dir/7_clonotypes_truth.txt`; do
    for p in "230117_VH01211_6_AAC7YLMM5" "M06463_0121new"; do
        for f in `ls -1q $p/*/*-TCR.clonotypes.new.txt`; do
            s=`basename $f | cut -d"-" -f1`
            cond=`cat $code_dir/library.txt | grep $s | cut -f2,3 | sed 's/ /_/'`
            gene=`echo $clono | cut -d":" -f1`
            v=`echo $clono | cut -d":" -f2`
            j=`echo $clono | cut -d":" -f4`
            cdr=`echo $clono | cut -d":" -f3`
            o=`cat $f | grep $v | grep $j | awk -F"\t" -v OFS="\t" -v g=$gene -v cdr=$cdr '$2==g && $5==cdr' | awk -F"\t" -v OFS=":" -v p="$p" -v c="$cond" '{print p"\t"c"\t"$2,$3,$5,$4"\t"$NF}' | head -1`
            nl=${#o}
            if [ $nl -le 1 ]; then
                echo "" | awk -v OFS="\t" -v p="$p" -v c="$cond" -v cl="$clono" '{print p,c,cl,0}'
            else
                echo $o | awk -v OFS="\t" -v p="$p" -v c="$cond" -v cl="$clono" '{print p,c,cl,$NF}'
            fi
        done
    done | sed 's/ /\t/g' | sed 's/%_/% /'
done >> umi_count_7_clonotypes.txt
