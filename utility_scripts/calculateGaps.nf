params.fasta = ""

fastaFile = Channel.fromPath(params.fasta)

process calculateGaps {
cache 'deep'
publishDir './Supporting_non-OGS_data/Gaps/',mode:copy,overwrite:true
input:
 file fastaFile
output:
 file "*_Ns2gff.gff3"
script:
"""
NEWNAME="${fastaFile}"
NEWNAME=\${NEWNAME%.fasta}
NEWNAME=\${NEWNAME%.fa}
NEWNAME=\${NEWNAME%.fnt}
seqkit locate --gtf -Pi -p "[N]+" ${fastaFile} | sed 's/location/gapped_bases/g' | sed 's/gene_id /ID=/g' | sed 's/; //g' | gawk 'match($0, /^.+gapped_bases\t([0-9]+)\t([0-9]+)\t.+$/, a) {print a[0] a[1] "_" a[2]}' | gt gff3 -tidy -sort -retainids -fixregionboundaries > ${NEWNAME}_Ns2gff.gff3
"""
}
