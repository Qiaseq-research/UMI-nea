mkdir -p log
for rep in 1 2 3; do
    for acc in "0.995" "0.95"; do
        if [ $acc == "0.995" ]; then
            umi_len=12
        else
            umi_len=50
        fi
        prefix=sim.rep$rep.acc$acc.pN100000.oN20.Len$umi_len
        python /Download/UMI-nea/benchmark/simulate_UMI_indel.py $prefix $umi_len $acc 100000 20 1 > log/$prefix.log
        join -1 2 -2 1 -o 1.1 1.2 1.3 1.4 2.2 <(join <(cat $prefix.out | cut -f1,2 | sort | uniq -c | awk '{print $2,$3,$1}' | sort -k1,1) <(cat $prefix.out | cut -f1 | sort | uniq -c | awk '{print $2,$1,NR-1}' | sort -k1,1) | sort -k2,2) <(cat $prefix.out | cut -f1 | sort -u | awk '{print $1,NR-1}' | sort -k1,1) | sort -k1,1 | awk '{print $1,$3,$5}' | awk '{for(i=1;i<=$2;i++){print $1,$3}}' > $prefix.truth.labels
    done
done

for i in sim*.out; do
    for thread in 12 24 48 60 72 84 96 ; do
        cat $i | cut -f1 | sort | uniq -c | awk '{print "1\t"$2"\t"$1}' | sort -k3,3nr > $i.input
        maxl=`cat $i | cut -f2 | awk '{print length($1)}' | sort -nr | head -1`
        acc=`basename $i | sed "s/.*acc//" | cut -f1,2 -d .`
        err_rate=`echo $acc | awk '{printf "%1.3f\n",(1-$1)}'`
        echo "UMI-nea $i thread=$thread len=$maxl err=$err_rate" >> UMI-nea.time
        { time -p bash -c "/Download/UMI-nea/UMI-nea/UMI-nea -i $i.input -o $i.thread$thread.clustered -l $maxl -t $thread -e $err_rate -a > log/$i.thread$thread.log"; }  2>> UMI-nea.time
    done
done