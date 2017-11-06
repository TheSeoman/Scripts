echo "snp\tchr\tpos" > snp.pos.tsv
awk '{print $2\t$1\t$3}' full_sorted > snp.pos.tsv        