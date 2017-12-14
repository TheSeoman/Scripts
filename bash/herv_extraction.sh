#in dir that contains rmsk.txt
awk '$11 ~ /ERV/ { print $6"\t"$7"\t"$8"\t"$11"\t"$2"\t"$10  }' rmsk.txt > hervS1.bed
awk '$13 ~ /ERV/ { print $6"\t"$7"\t"$8"\t"$11"\t"$2"\t"$10  }' rmsk.txt > hervS2.bed
awk '$11 ~ /HERV/  { print $6"\t"$7"\t"$8"\t"$11"\t"$2"\t"$10  }' rmsk.txt > hervS3.bed

awk '$11 ~ /LTR/ { print $6"\t"$7"\t"$8"\t"$11"\t"$2"\t"$10  }' rmsk.txt > ltr.bed