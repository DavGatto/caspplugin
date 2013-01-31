learningOptions=( deletion forward backward jumpforward cc wcc range )

echo -ne "instanceName" > result.csv

for learningOption in "${learningOptions[@]}"
do
	echo -ne ", " $learningOption >> result.csv
done

echo "" >>result.csv

for instance in simple.hex
do
	echo -ne $instance >>result.csv
	for learningOption in "${learningOptions[@]}"
	do
		echo "Running instance" $instance "with learning" $learningOption
		TIMEFORMAT=%R time1=$({ time (dlvhex2 --extlearn --heuristics=monolithic --csplearning=$learningOption --silent $instance 2> /dev/null) >/dev/null; } 2>&1 )
		echo -ne "," $time1 >>result.csv
	done

	echo "" >> result.csv

	#echo $instance "," $time1 "," $time2 "," $time3 "," $time4 "," $time5 "," $time6 "," $time7 >> result.csv
done
