params.fasta = ""

fastaFile = Channel.fromPath(params.fasta)
fastaFile.into{fastaFile_ch1 ; fastaFile_ch2}
fasta_name = Channel.value(fastaFile_ch2.getVal().toFile().name.tokenize(".").get(0))

//process splitFasta {
//conda "ucsc-fasplit"
//input:
// file fasta from fastaFile_ch1
// val fasta_name
//output:
// file "split*" into splitFiles
//tag "${fasta}"
//script:
//"""
// faSplit about ${fasta} 1000000 split_${fasta_name}
//"""
//}


process calculateGaps {
cache 'deep'
cpus 1
conda "seqkit"
input:
 file fasta from fastaFile_ch1.splitFasta(size:'1.MB',file:true)
output:
 file "*_located_gaps.gtf" into gtfFiles
tag "${fasta}"
script:
"""
seqkit locate -j ${task.cpus} --gtf -Pi -p "N+" ${fasta} | sed 's/location/gapped_bases/g' | sed 's/gene_id /ID=/g' | sed 's/; //g' > ${fasta}_located_gaps.gtf
"""
}

process convertSeqkitGTFtoGFF {
conda "gawk"
input:
 file splitGtfFile from gtfFiles
output:
 file "${splitGtfFile}_gaps.gff" into gffFiles
tag "${splitGtfFile}"
 script:
"""
NEWNAME="${splitGtfFile}"
NEWNAME=\${NEWNAME%_located_gaps.gtf}
##Note the extra backslashes for Nextflow
cat ${splitGtfFile} | gawk 'match(\$0, /^(.+)\tSeqKit\tgapped_bases\t([0-9]+)\t([0-9]+)\t.+\$/, a) {print a[0] a[1] ":" a[2] "-" a[3]}' > ${splitGtfFile}_gaps.gff
"""
}

process combineAndPublish {
publishDir './Supporting_non-OGS_data/Gaps/',mode:'copy',overwrite:true
conda "genometools"
input:
 file gf from gffFiles.collect()
 val fasta_name
output:
 file "${fasta_name}_Ns2gff.gff3"
tag "${fasta_name}_Ns2gff.gff3"
script:
"""
ls -f1 | grep ".gff" > files.txt
cat files.txt | xargs cat | gt gff3 -tidy -sort -retainids -fixregionboundaries > ${fasta_name}_Ns2gff.gff3
"""
}
