for instance in learning_test_7.hex
do
	echo "Running instance" $instance " with no learning"
	echo $(timeout 3m /usr/bin/time --format "%e" dlvhex2 --extlearn --heuristics=monolithic --csplearning=none --silent $instance >>/dev/null)

	echo "Running instance" $instance " with deletion filtering learning"
	echo $(timeout 3m /usr/bin/time --format "%e" dlvhex2 --extlearn --heuristics=monolithic --csplearning=deletion --silent $instance >>/dev/null)

	echo "Running instance" $instance " with forward filtering learning"
	echo $(timeout 3m /usr/bin/time --format "%e" dlvhex2 --extlearn --heuristics=monolithic --csplearning=forward --silent $instance >>/dev/null)

	echo "Running instance" $instance " with backward filtering learning"
	echo $(timeout 3m /usr/bin/time --format "%e" dlvhex2 --extlearn --heuristics=monolithic --csplearning=backward --silent $instance >>/dev/null)

	echo "Running instance" $instance " with connected components filtering learning"
	echo $(timeout 3m /usr/bin/time --format "%e" dlvhex2 --extlearn --heuristics=monolithic --csplearning=cc --silent $instance >>/dev/null)

	echo "=========================================================================="
done
