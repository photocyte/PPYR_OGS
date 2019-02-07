#! /bin/bash

echo "This is a script intended to take a Official Gene Set file which is very nearly in a standardized format (e.g., you edited the existing OGS GFF3 file with a simple text editor), and reformat that GFF3 plus produce the derivative FASTA files to ensure they are in a standardized format"
set -e
set -o pipefail

if [ -z "$1" ] 
then
echo "First argument should be the path to the GFF3 file to process."
echo "The naming of the GFF3 should be in the style PPYR_OGS1.0.gff3, where PPYR = the Locus Tag Prefix for the genome in question, and 1.0 = the release # you want to build, e.g. 1.0, 1.1, 1.2, 1.3... If you edited the existing GFF3 file in the repository, conviniently it should already have this naming scheme"
exit
fi

if [ ! -f "Genome_release.fa" ]; then
    echo "Genome_release.fa not found. Please symlink the most recent/relevant genome release to this file. This should be the FASTA file which the GFF3 file is assosciated with"
    exit
fi

BASE=$(basename $1)
BASE=${BASE%.gff3}

echo "Adding start/stop codons with GAG..."
rm -rf gag_output
mkfifo tmp0.gff3
cat $1 | grep -v "stop_codon" | grep -v "start_codon" > tmp0.gff3 &
python2 /lab/solexa_weng/testtube/GAG-2.0.1/gag.py --fix_start_stop -f ./Genome_release.fa -g tmp0.gff3
rm -f tmp0.gff3
echo "Done adding codons"

echo "Sorting with gt..."
cat gag_output/genome.gff gag_output/genome.ignored.gff | grep -v "#" | gt gff3 -tidy -sort -retainids > tmp.${BASE}.gt.gff3
echo "Done gt sorting."
echo "Sorting with igvtools..."
igvtools sort tmp.${BASE}.gt.gff3 tmp.${BASE}.gt.igv.gff3
echo "Done sorting."

echo "Extracting CDS features..."
gt extractfeat -join -seqid -usedesc -retainids -coords -type CDS -seqfile Genome_release.fa tmp.${BASE}.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" | gzip > ${BASE}.CDS.fa.gz
echo "Extracting peptide features..."
gt extractfeat -join -seqid -usedesc -retainids -coords -type CDS -translate -gcode 1 -seqfile Genome_release.fa tmp.${BASE}.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" | gzip > ${BASE}.pep.fa.gz
echo "Extracting mRNA features..."
gt extractfeat -join -seqid -usedesc -retainids -coords -type exon -seqfile Genome_release.fa tmp.${BASE}.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" | gzip > ${BASE}.mRNA.fa.gz
echo "Extracting gene features..."
gt extractfeat -seqid -usedesc -retainids -coords -type gene -seqfile Genome_release.fa tmp.${BASE}.gt.gff3 | seqkit replace -p "\(joined\)|\(translated\)" -r "" | gzip > ${BASE}.gene.fa.gz

echo "Recording FASTA file checksums..."
rm -f fasta-checksums.txt
echo "${BASE}.CDS.fa.gz" > fasta-checksums.txt
./utility_scripts/checksum_FASTA_files.sh ${BASE}.CDS.fa.gz >> fasta-checksums.txt
echo "${BASE}.pep.fa.gz" >> fasta-checksums.txt
./utility_scripts/checksum_FASTA_files.sh ${BASE}.pep.fa.gz >> fasta-checksums.txt
echo "${BASE}.mRNA.fa.gz" >> fasta-checksums.txt
./utility_scripts/checksum_FASTA_files.sh ${BASE}.mRNA.fa.gz >> fasta-checksums.txt
echo "${BASE}.gene.fa.gz" >> fasta-checksums.txt
./utility_scripts/checksum_FASTA_files.sh ${BASE}.gene.fa.gz >> fasta-checksums.txt

echo "Overwriting original file with igvtools sorted version..."
mv -f tmp.${BASE}.gt.igv.gff3 $1
rm -f tmp.${BASE}.gt.gff3

