echo "------------------------------------"
seqkit stat $1
echo "idchk":
seqkit seq -n -i $1 | sort | openssl md5
echo "seqchk":
seqkit seq -u $1 | seqkit sort -s | grep -v ">" | openssl md5

