#!/bin/bash


mkdir exe
mkdir test/exe


for name in CreateMatrix DevideUni MergePriorProbs Multi2Uni PriorC; do
	cd $name
	chmod "u=rwx" ./compile.sh
	./compile.sh
	echo "${name} Сompiled"
	cp $name ../exe/
	mv $name ../test/exe/
	cd ..
done

cd test
./test.sh
cd ..
