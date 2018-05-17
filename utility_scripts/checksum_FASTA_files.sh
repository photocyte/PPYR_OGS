

seqkit sort -s $1 | grep -v ">" | openssl md5
