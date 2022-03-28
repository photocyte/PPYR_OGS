nextflow.enable.dsl=2

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

process EVM_gff3_simple_validator {
conda "evidencemodeler"
input: 
 path doubleSorted
echo true
script:
"""
#echo \$CONDA_PREFIX
#echo \$CONDA_PREFIX_1
\$CONDA_PREFIX/opt/evidencemodeler-*/EvmUtils/gff3_gene_prediction_file_validator.pl <(cat ${doubleSorted} | grep 
"""
}

process gfacs_validator {
conda "perl-gfacs coreutils"
input: 
 path doubleSorted
echo true
shell:
'''
GFACS_PATH="$(dirname $(readlink -f $(which gfacs.pl)))"
ln -s $GFACS_PATH/* .
mkdir output
gfacs.pl -f refseq_gff -O output !{doubleSorted}
'''
}

workflow doubleSort_wf {
take: gff
main: 
 gtSort(gff)
 igvtoolsSort(gtSort.out)
 //EVM_gff3_simple_validator(igvtoolsSort.out)
 gfacs_validator(igvtoolsSort.out)
emit:
 igvtoolsSort.out
}

workflow {
 Channel.fromPath(params.gff).set{gffInput}
 doubleSort_wf(gffInput)
}
