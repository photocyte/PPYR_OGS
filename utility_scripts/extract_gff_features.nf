nextflow.enable.dsl=2

process gffToScaffoldList {
input:
 path gff_ch1
output:
 path "gff_scaffolds.txt"
tag "${gff_ch1}"
script:
"""
cut -f 1 ${gff_ch1} | grep -Pv "^#" | uniq | sort | uniq > gff_scaffolds.txt
"""
}

process partitionGenome {
cache 'deep'
conda "ucsc-fasplit seqkit"
input:
 path fasta_ch1
 path gffScaffolds_ch
output:
 path "split/*.fa"
tag "${fasta_ch1}"
script:
"""
mkdir split
faSplit about <(cat ${fasta_ch1} | seqkit grep -f ${gffScaffolds_ch}) 10000000 split/
##If the scaffold isn't in the gff, just delete it.
##for f in ./split/*.fa
##do
##sort <(cat ${gffScaffolds_ch}) <(seqkit fx2tab -ni \${f}) | uniq -d > \${f}.dups.txt
##done
"""
}


process filterGFF {
conda "seqkit genometools-genometools grep"
cache 'deep'
input:
 tuple path(chunk),path(gff)
output:
 tuple path("${chunk}"),path("filtered.gff3") optional true
tag "${chunk}"
script:
"""
seqkit fx2tab -n -i ${chunk} | tr -s "\t" | sed 's/^/\\^/g' > fasta_records.txt ##A trailing \t is actually important for accurate grepping
grep -Gf fasta_records.txt ${gff} | gt gff3 -tidy -sort -retainids > filtered.gff3

if [[ -s filtered.gff3 ]]
then
echo "pass"
else
rm -f filtered.gff3
fi
"""
}

process makeGtIndex {
conda "genometools-genometools"

input:
 tuple path(chunk),path(filteredGff)
output:
 tuple path("${chunk}"),path("${chunk}.des"),path("${chunk}.sds"),path("${chunk}.md5"),path("${chunk}.esq"),path("${chunk}.ssp"),path("${chunk}.ois"),path("${filteredGff}")

script:
"""
gt seq ${chunk}
if [ ! -f "${chunk}.ssp" ]
then
    touch "${chunk}.ssp"
fi
##Below was used for testing. Not needed anymore
##chmod u-w ${chunk}.des ${chunk}.sds ${chunk}.md5 ${chunk}.esq ${chunk}.ssp ${chunk}.ois
"""
}

process extractFeatures {
conda "seqkit genometools-genometools"
input:
 tuple val(featureType),path(fastaChunk),path(des),path(sds),path(md5),path(esq),path(ssp),path(ois),path(filteredGff)

output:
 tuple val("${featureType}"),path("${fastaChunk}.${featureType}.fa.gz")
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
elif [ "$featureType" == "indExon" ]
then
echo "Extracting individual exon features..."
  gt extractfeat -seqid -usedesc -retainids -coords -type exon -seqfile ${fastaChunk} ${filteredGff} | seqkit replace -p "\\(joined\\)|\\(translated\\)" -r "" | gzip > ${fastaChunk}.indExon.fa.gz
elif [ "$featureType" == "indCDS" ]
then
echo "Extracting individual CDS features..."
  gt extractfeat -seqid -usedesc -retainids -coords -type CDS -translate -gcode 1 -seqfile ${fastaChunk} ${filteredGff} | seqkit replace -p "\\(joined\\)|\\(translated\\)" -r "" | gzip > ${fastaChunk}.indCDS.fa.gz
elif [ "$featureType" == "gene" ]
then
echo "Extracting gene features..."
  gt extractfeat -seqid -usedesc -retainids -coords -type gene -seqfile ${fastaChunk} ${filteredGff} | seqkit replace -p "\\(joined\\)|\\(translated\\)" -r "" | gzip > ${fastaChunk}.gene.fa.gz
fi 
"""
}


process mergeFastas {
publishDir './' , mode:'link' , overwrite: true
//conda "seqkit"
input:
 tuple val(fn),val(featureType),path(allFiles)
output:
 path "${fn}.${featureType}.fa.gz"
tag "${fn}.${featureType}.fa.gz"
script:
 """
 ls -f1 | grep ".fa.gz" > files.txt
 cat files.txt | xargs cat | seqkit sort -n | gzip > ${fn}.${featureType}.fa.gz
 """
}

workflow {
params.splitBy = 5
fasta_ch = Channel.fromPath(params.fasta)
gff_ch = Channel.fromPath(params.gff)



gffToScaffoldList(gff_ch)


fastaChunks = partitionGenome(fasta_ch,gffToScaffoldList.out)
fastaChunks.flatten().combine(gff_ch).set{combinedCmds}

filterGFF(combinedCmds)
indexedChunks = makeGtIndex(filterGFF.out)

feature_types = Channel.from( "CDS", "pep", "mRNA", "gene" , "indCDS" , "indExon" )
feature_types.combine(indexedChunks).set{extractCmds}
extractFeatures(extractCmds)
extractFeatures.out.groupTuple().set{groupedFeatureFastas}

mergeCmds = fasta_ch.map{ it.getFileName() }.combine(groupedFeatureFastas)
mergeFastas(mergeCmds)

}
