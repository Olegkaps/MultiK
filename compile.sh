#!/bin/bash


for name in CreateMatrix DevideUni MergePriorProbs Multi2Uni PriorC; do
	cd $name
	./compile.sh
	echo "${name} Сompiled"
	cp $name ../exe/
	mv $name ../test/exe/
	cd ..
done

cd test
./test.sh
cd ..
