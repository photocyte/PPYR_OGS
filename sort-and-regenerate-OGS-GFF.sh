#! /bin/bash

if [ -z "$1" ] 
then
echo "First argument should be the path to the GFF3 file to process"
exit
fi

BASE=$(basename $1)
BASE=${BASE%.gff3}

echo "Sorting with gt..."
gt gff3 -tidy -sort -retainids $1 > tmp.gt.gff3

echo "Extracting CDS features..."
gt extractfeat -join -seqid -usedesc -retainids -coords -type CDS -seqfile Ppyr_Genome_release.fa tmp.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" > ${BASE}.CDS.fa
echo "Extracting peptide features..."
gt extractfeat -join -seqid -usedesc -retainids -coords -type CDS -translate -gcode 1 -seqfile Ppyr_Genome_release.fa tmp.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" > ${BASE}.pep.fa
echo "Extracting mRNA features..."
gt extractfeat -join -seqid -usedesc -retainids -coords -type mRNA -seqfile Ppyr_Genome_release.fa tmp.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" > ${BASE}.mRNA.fa
echo "Extracting gene features..."
gt extractfeat -join -seqid -usedesc -retainids -coords -type CDS -seqfile Ppyr_Genome_release.fa tmp.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" > ${BASE}.gene.fa

echo "Sorting with igvtools..."
/Applications/IGVTools/igvtools sort tmp.gt.gff3 tmp.gt.igv.gff3
echo "Overwriting original file with sorted version..."
mv -f tmp.gt.igv.gff3 $1
rm -f tmp.gt.gff3

