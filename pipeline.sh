
cd exe


echo "==== Create Matrix start ===="

./CreateMatrix --input_format contacts -i "../input_data/bam-to-contacts-test.txt" \
 -o "../results/multi_main_parsed.txt" "../results/uni_main_parsed.txt" \
 -r 100000 -a "../input_data/main_test/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff" \
 --chr_ids "../input_data/main_test/chromosomes.txt"

echo "==== Create Matrix files ===="

echo "Multi Parsed"
head "../results/multi_main_parsed.txt"
echo "Uni Parsed"
head "../results/uni_main_parsed.txt"

echo "==== Create Matrix end ===="



echo "==== Devide Uni start ===="
./DevideUni -i "../results/uni_main_parsed.txt" -1 "../results/uni_main_1.txt" -2 "../results/uni_main_2.txt"
echo "==== Devide Uni end ===="



echo "==== Prior start ===="

./PriorC -i "../results/uni_main_1.txt" -f "../input_data/frags.txt" -r 100000 -o "../results/spline_main_1.txt" -b 200 -d "../results/density.txt"
./PriorC -i "../results/uni_main_2.txt" -f "../input_data/frags.txt" -r 100000 -o "../results/spline_main_2.txt" -b 200

echo "==== Prior files ===="
echo "spline 1"
head "../results/spline_main_1.txt"
echo "spline 2"
head "../results/spline_main_2.txt"
echo "density"
head "../results/density.txt"

echo "==== Prior end ===="




echo "==== Merge Prior start ===="
./MergePriorProbs -1 "../results/spline_main_1.txt" -2 "../results/spline_main_2.txt" -o "../results/spline_main_merged.txt"
echo "==== Merge Prior files ===="
echo "spline"
head "../results/spline_main_merged.txt"
echo "==== Merge Prior end ===="



echo "==== Multi to Uni start ===="
./Multi2Uni -p "../results/spline_main_merged.txt" -m "../results/multi_main_parsed.txt" -f "../results/main_probs.txt" -d "../results/density.txt"
echo "==== Multi to Uni result ===="
head -n 200 ../results/main_probs.txt
echo "==== Multi to Uni end ===="
