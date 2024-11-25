
for dir in "devide_probs" "em" "merge_probs" "prior_probs"; do
	cd ${dir}
	rm -f res/*
	pytest
	rm -rf  __pycache__
	rm -f answ.*
	cd ..
done

rm -rf  __pycache__
