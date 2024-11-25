declare -A progs

progs["DevideUni"]="DevideUni"
progs["MergePriorProbs"]="MergePriorProbs"
progs["PriorC"]="PriorC"
progs["Multi2Uni"]="Multi2Uni"

for dir in DevideUni MergePriorProbs Multi2Uni PriorC; do
	cd $dir
	./compile.sh
	cp ${progs[${dir}]} ../exe/
	mv ${progs[${dir}]} ../test/exe/
	cd ..
done

cd test
./test.sh
