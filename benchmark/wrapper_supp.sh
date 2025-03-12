for f in "10000" "100000" "200000"; do
	for l in "12" "24" "48"; do
                for a in "0.999" "0.995" "0.97"; do
			n_t=2
			dist=0
			if [ $l -eq 12 ] && [ $a == "0.97" ]; then
				dist=1
			fi
			echo "$l $a $f $dist"
			bash /code/compare_3_tools.sh $l $a $f 20 $dist $n_t 0
		done
	done
done
