cat $1 | awk '{print $2}' | sort | uniq -c | awk '{print $1}' > sim.umi.count
fn=`wc -l sim.umi.count`
m=`cat sim.umi.count | awk -v n=0 '{n+=$1}END{print n/NR}'`
v=`cat sim.umi.count | awk -v a=0 -v m=$m '{s=($1-m)*($1-m);a+=s}END{print int(a/(NR-1))}'`
r=`echo "($m*$m)/($v-$m)" | bc -l`
p=`echo "$r/($r+$m)" | bc -l`
echo "$1 $fn $m $v $r $p"
