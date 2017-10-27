#in dir that contains rmsk.txt
awk '$11 ~ /ERV/ { for(i=6;i<=13;i++)printf "%s\t",$i;printf "\n"  }' rmsk.txt > hervS1.bed
awk '$13 ~ /ERV/ { for(i=6;i<=13;i++)printf "%s\t",$i;printf "\n"  }' rmsk.txt > hervS2.bed
awk '$11 ~ /HERV/ { for(i=6;i<=13;i++)printf "%s\t",$i;printf "\n"  }' rmsk.txt > hervS3.bed