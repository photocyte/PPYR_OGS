nextflow.preview.dsl=2
Channel.fromPath(params.gff).set{gffInput}

process gtSort {
conda "genometools-genometools"
cache 'deep'
input:
 file gffInput
output:
 file "tmp.${gffInput}.gt.gff3" into gtSorted

script:
"""
echo "Sorting with gt..."
cat ${gffInput} | grep -v "#" | gt gff3 -tidy -sort -retainids > tmp.${gffInput}.gt.gff3
echo "Done gt sorting."
"""
}

process igvtoolsSort {
conda "igvtools"
publishDir './' , mode:'copy'
input:
 file gtSorted
output:
 file "tmp.${params.base}.gt.igv.gff3" into doubleSorted
script:
"""
echo "Sorting with igvtools..."
igvtools sort ${gtSorted} tmp.${params.base}.gt.igv.gff3
echo "Done sorting."
"""
}

process EVM_gff3_validator {
conda "evidencemodeler"
input: 
 file doubleSorted
echo true
script:
"""
#echo \$CONDA_PREFIX
#echo \$CONDA_PREFIX_1
\$CONDA_PREFIX/opt/evidencemodeler-*/EvmUtils/gff3_gene_prediction_file_validator.pl ${doubleSorted}
"""
}

workflow doubleSort_wf {
take: gff
 gtSort(gff)
 igvtoolsSort(gtSort.out)
 EVM_gff3_validator(igvtoolsSort.out)

emit:
 igvtoolsSort.out
}

workflow {
 doubleSort_wf(params.gff)
}
