#kallisto index -i Ppyr1.3_Trinity-PASA_stranded-DCGM.mRNA ../Ppyr1.3_Trinity-PASA_stranded-DCGM.mRNA.fa.gz

REFNAME=Ppyr1.3_Trinity-PASA_stranded-DCGM.mRNA

BS=100

rm -rf non_strand_specific_results
mkdir non_strand_specific_results
for F_READ in ./non_strand_specific_reads/*_1.fq.gz
do
SAMPLE_NAME=$(basename $F_READ)
SAMPLE_NAME=${SAMPLE_NAME%_1.fq.gz}
R_READ=${F_READ%_1.fq.gz}_2.fq.gz
bsub -n 8 -R "span[hosts=1]" "kallisto quant -t 8 -b $BS -i $REFNAME -o non_strand_specific_results/$SAMPLE_NAME $F_READ $R_READ 1>$SAMPLE_NAME.kallisto.stdout.log 2>$SAMPLE_NAME.kallisto.stderr.log"
done

rm -rf strand_specific_results
mkdir strand_specific_results
for F_READ in ./strand_specific_reads/*_1.fq.gz
do
SAMPLE_NAME=$(basename $F_READ)
SAMPLE_NAME=${SAMPLE_NAME%_1.fq.gz}
R_READ=${F_READ%_1.fq.gz}_2.fq.gz
bsub -n 8 -R "span[hosts=1]" "kallisto quant --rf-stranded -t 8 -b $BS -i $REFNAME -o strand_specific_results/$SAMPLE_NAME $F_READ $R_READ 1>$SAMPLE_NAME.kallisto.stdout.log 2>$SAMPLE_NAME.kallisto.stderr.log"
done
