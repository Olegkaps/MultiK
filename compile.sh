#!/bin/bash

for name in DevideUni MergePriorProbs Multi2Uni PriorC; do
	cd $name
	./compile.sh
	cp $name ../exe/
	mv $name ../test/exe/
	cd ..
done

cd test
./test.sh
