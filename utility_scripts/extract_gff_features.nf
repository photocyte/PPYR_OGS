params.splitBy = 5
fasta_ch = Channel.fromPath(params.fasta)
gff_ch = Channel.fromPath(params.gff)

fasta_ch.into{fasta_ch1 ; fasta_ch2 ; fasta_ch3}
gff_ch.into{gff_ch1 ; gff_ch2}

fasta_ch3.println()
gff_ch2.println()
//fasta_name = Channel.from("Genome_release.fa")
fasta_name = Channel.from(fasta_ch2.getVal().toFile().name)

fasta_ch1.splitFasta(by:params.splitBy,file:true).set{fastaChunks}

fastaChunks.combine(gff_ch1).set{combinedCmds}

process filterGFF {
conda "seqkit genometools"
input:
 set file(chunk),file(gff) from combinedCmds
output:
 set file("${chunk}"),file("filtered.gff3") optional true into filteredChunks
tag "${chunk}"
script:
"""
seqkit fx2tab -n -i ${chunk} | tr -s "\t" > fasta_records.txt ##A trailing \t is actually important for accurate grepping
cat ${gff} | grep -f fasta_records.txt | gt gff3 -tidy -sort -retainids > filtered.gff3

if [[ -s filtered.gff3 ]]
then
echo "pass"
else
rm -f filtered.gff3
fi
"""
}

process makeGtIndex {
conda "genometools"

input:
 set file(chunk),file(filteredGff) from filteredChunks
output:
 set file("${chunk}"),file("${chunk}.des"),file("${chunk}.sds"),file("${chunk}.md5"),file("${chunk}.esq"),file("${chunk}.ssp"),file("${chunk}.ois"),file("${filteredGff}") into indexedChunks

script:
"""
gt seq ${chunk}
##Below was used for testing. Not needed anymore
##chmod u-w ${chunk}.des ${chunk}.sds ${chunk}.md5 ${chunk}.esq ${chunk}.ssp ${chunk}.ois
"""
}

feature_types = Channel.from( "CDS", "pep", "mRNA", "gene" )
feature_types.combine(indexedChunks).set{extractCmds}

process extractFeatures {
conda "seqkit genometools"
input:
 set val(featureType),file(fastaChunk),file(des),file(sds),file(md5),file(esq),file(ssp),file(ois),file(filteredGff) from extractCmds

output:
 set val("${featureType}"),file("${fastaChunk}.${featureType}.fa.gz") into extractedFeatureFastas 
tag "${featureType} ${fastaChunk}"
script:
"""
##Note extra backslashes for nextflow
if [ "$featureType" == "CDS" ]
then
 echo "Extracting CDS features..."
 gt extractfeat -join -seqid -usedesc -retainids -coords -type CDS -seqfile ${fastaChunk} ${filteredGff} | seqkit replace -p "\\(joined\\)|\\(translated\\)" -r "" | gzip > ${fastaChunk}.CDS.fa.gz
elif [ "$featureType" == "pep" ] 
then
echo "Extracting peptide features..."
  gt extractfeat -join -seqid -usedesc -retainids -coords -type CDS -translate -gcode 1 -seqfile ${fastaChunk} ${filteredGff} | seqkit replace -p "\\(joined\\)|\\(translated\\)" -r "" | gzip > ${fastaChunk}.pep.fa.gz
elif [ "$featureType" == "mRNA" ]
then
echo "Extracting mRNA features..."
  gt extractfeat -join -seqid -usedesc -retainids -coords -type exon -seqfile ${fastaChunk} ${filteredGff} | seqkit replace -p "\\(joined\\)|\\(translated\\)" -r "" | gzip > ${fastaChunk}.mRNA.fa.gz
elif [ "$featureType" == "gene" ]
then
echo "Extracting gene features..."
  gt extractfeat -seqid -usedesc -retainids -coords -type gene -seqfile ${fastaChunk} ${filteredGff} | seqkit replace -p "\\(joined\\)|\\(translated\\)" -r "" | gzip > ${fastaChunk}.gene.fa.gz
fi 
"""
}

extractedFeatureFastas.groupTuple().set{groupedFeatureFastas}

fasta_name.combine(groupedFeatureFastas).set{totalGroup}

process mergeFastas {
//publishDir './' , mode:'copy'
conda "seqkit"
input:
 set val(fn),val(featureType),file(allFiles) from totalGroup
output:
 file "${fn}.${featureType}.fa.gz"
tag "${fn}.${featureType}.fa.gz"
script:
 """
 ls -f1 | grep ".fa.gz" > files.txt
 cat files.txt | xargs cat | seqkit sort -n | gzip > ${fn}.${featureType}.fa.gz
 """
}
