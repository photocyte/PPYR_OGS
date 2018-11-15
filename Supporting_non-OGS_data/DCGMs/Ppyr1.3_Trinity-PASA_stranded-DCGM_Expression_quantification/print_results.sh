for LIB in ./*results/*/abundance.tsv
do
echo $LIB
cat $LIB | grep -P "asmbl_27914.p1|asmbl_27909.p1|asmbl_27913.p1" | cut -f 1,5 | sed 's/Ppyr1.3_Trinity-PASA_stranded-DCGM_transdecoder_//g'
done
