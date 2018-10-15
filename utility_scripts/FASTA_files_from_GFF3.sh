#! /bin/bash

echo "This is a script intended to take a Official Gene Set file which is very nearly in a standardized format (e.g., you edited the existing OGS GFF3 file with a simple text editor), and reformat that GFF3 plus produce the derivative FASTA files to ensure they are in a standardized format"

if [ -z "$1" ];
then
echo "First argument should be the path to the GFF3 file to process."
echo "The naming of the GFF3 should be in the style PPYR_OGS1.0.gff3, where PPYR = the Locus Tag Prefix for the genome in question, and 1.0 = the release # you want to build, e.g. 1.0, 1.1, 1.2, 1.3... If you edited the existing GFF3 file in the repository, conviniently it should already have this naming scheme"
exit
fi

if [ -z "$2" ]; then
    echo "Second argument should be path to the genome reference file"
    exit
fi

BASE=$(basename $1)
if [ ${BASE: -3} == ".gz" ]
then
  BASE=${BASE%.gff3.gz}
  echo "Sorting with gt..."
  zcat $1 | grep -v "#" | gt gff3 -tidy -sort -retainids > tmp.${BASE}.gt.gff3
  echo "Done sorting."
else
  BASE=${BASE%.gff3}
  echo "Sorting with gt..."
  cat $1 | grep -v "#" | gt gff3 -tidy -sort -retainids > tmp.${BASE}.gt.gff3
  echo "Done sorting."
fi

echo "Extracting CDS features..."
gt extractfeat -join -seqid -usedesc -retainids -coords -type CDS -seqfile $2 tmp.${BASE}.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" | gzip > ${BASE}.CDS.fa.gz
echo "Extracting peptide features..."
gt extractfeat -join -seqid -usedesc -retainids -coords -type CDS -translate -gcode 1 -seqfile $2 tmp.${BASE}.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" | gzip > ${BASE}.pep.fa.gz
echo "Extracting mRNA features..."
gt extractfeat -join -seqid -usedesc -retainids -coords -type exon -seqfile $2 tmp.${BASE}.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" | gzip > ${BASE}.mRNA.fa.gz
echo "Extracting gene features..."
gt extractfeat -seqid -usedesc -retainids -coords -type gene -seqfile $2.fa tmp.${BASE}.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" | gzip > ${BASE}.gene.fa.gz

echo "Sorting with igvtools..."
igvtools sort tmp.${BASE}.gt.gff3 tmp.${BASE}.gt.igv.gff3
echo "gzip compressing..."
gzip tmp.${BASE}.gt.igv.gff3

echo "Overwriting original file with sorted & compressed version..."
mv -f tmp.${BASE}.gt.igv.gff3.gz ${BASE}.gff3.gz
rm -f tmp.${BASE}.gt.igv.gff3
rm -f tmp.${BASE}.gt.gff3

