nextflow.preview.dsl=2

process gtSort {
conda "genometools-genometools"
cache 'deep'
input:
 path gffInput
output:
 path "tmp.${gffInput}.gt.gff3"

script:
"""
echo "Sorting with gt..."
cat ${gffInput} | grep -v "#" | gt gff3 -tidy -sort -retainids > tmp.${gffInput}.gt.gff3
echo "Done gt sorting."
"""
}

process igvtoolsSort {
conda "igvtools"
publishDir './results/doubleSort' , mode:'link'
input:
 path gtSorted
output:
 path "tmp.${gtSorted}.gt.igv.gff3"
script:
"""
echo "Sorting with igvtools..."
igvtools sort ${gtSorted} tmp.${gtSorted}.gt.igv.gff3
echo "Done sorting."
"""
}

process EVM_gff3_validator {
conda "evidencemodeler"
input: 
 path doubleSorted
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
main: 
 gtSort(gff)
 igvtoolsSort(gtSort.out)
 EVM_gff3_validator(igvtoolsSort.out)
emit:
 igvtoolsSort.out
}

workflow {
 Channel.fromPath(params.gff).set{gffInput}
 doubleSort_wf(gffInput)
}
