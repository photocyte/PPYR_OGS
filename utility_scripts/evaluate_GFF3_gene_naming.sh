
cat $1 | grep "gene" | cut -d "=" -f 2 | cut -d ";" -f 1 | sort | uniq -c
