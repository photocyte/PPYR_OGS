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

##Removed this, as GAG messed with my locus_tag attribute?
#########################
#rm -rf gag_output
#echo "Adding start/stop codons with GAG..."
#cat $1 | grep -v "stop_codon" | grep -v "start_codon" > tmp0.gff3
#python2 /lab/solexa_weng/testtube/GAG-2.0.1/gag.py --fix_start_stop -f ./Genome_release.fa -g tmp0.gff3
#rm -f tmp0.gff3
#echo "Done adding codons"
#########################

echo "gt and igv sorting"
nextflow run ./utility_scripts/doubleSort.nf --gff $1 --base ${BASE} -resume

echo "extracting features with nextflow extract_gene_features.nf"
nextflow run ./utility_scripts/extract_gff_features.nf --gff tmp.${BASE}.gt.igv.gff3 --fasta Genome_release.fa -resume -with-trace

echo "Renaming nextflow files..."
mv Genome_release.fa.CDS.fa.gz ${BASE}.CDS.fa.gz
mv Genome_release.fa.pep.fa.gz ${BASE}.pep.fa.gz
mv Genome_release.fa.mRNA.fa.gz ${BASE}.mRNA.fa.gz
mv Genome_release.fa.gene.fa.gz ${BASE}.gene.fa.gz
##echo "Deleting nextflow work directory..."
##mv work DELETEME
##mkdir tmp_empty_dir
##time rsync -a --delete tmp_empty_dir/ DELETEME/
##rm -rf DELETEME tmp_empty_dir

echo "Recording FASTA file checksums..."
rm -f fasta-checksums.txt
echo "#########" > fasta-checksums.txt
seqkit stat ${BASE}.CDS.fa.gz >> fasta-checksums.txt
./utility_scripts/checksum_FASTA_files.sh ${BASE}.CDS.fa.gz >> fasta-checksums.txt
md5sum ${BASE}.CDS.fa.gz >> fasta-checksums.txt
echo "#########" >> fasta-checksums.txt
seqkit stat ${BASE}.pep.fa.gz >> fasta-checksums.txt
./utility_scripts/checksum_FASTA_files.sh ${BASE}.pep.fa.gz >> fasta-checksums.txt
md5sum ${BASE}.pep.fa.gz >> fasta-checksums.txt
echo "#########" >> fasta-checksums.txt
seqkit stat ${BASE}.mRNA.fa.gz >> fasta-checksums.txt
./utility_scripts/checksum_FASTA_files.sh ${BASE}.mRNA.fa.gz >> fasta-checksums.txt
md5sum ${BASE}.mRNA.fa.gz >> fasta-checksums.txt
echo "#########" >> fasta-checksums.txt
seqkit stat ${BASE}.gene.fa.gz >> fasta-checksums.txt
./utility_scripts/checksum_FASTA_files.sh ${BASE}.gene.fa.gz >> fasta-checksums.txt
md5sum ${BASE}.gene.fa.gz >> fasta-checksums.txt
echo "#########" >> fasta-checksums.txt
echo "Genome_release.fa" >> fasta-checksums.txt
echo `basename $(readlink -f Genome_release.fa)` >> fasta-checksums.txt
seqkit stat Genome_release.fa >> fasta-checksums.txt
./utility_scripts/checksum_FASTA_files.sh Genome_release.fa >> fasta-checksums.txt
md5sum Genome_release.fa >> fasta-checksums.txt

echo "Overwriting original file with igvtools sorted version..."
mv -f tmp.${BASE}.gt.igv.gff3 $1
rm -f tmp.${BASE}.gt.gff3

