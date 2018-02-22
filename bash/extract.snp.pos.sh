echo "snp\tchr\tpos" > snp.pos.tsv
awk '{print $2"\t"$1"\t"$3"\t"$4"\t"$5}' full_sorted >> snp.pos.tsv