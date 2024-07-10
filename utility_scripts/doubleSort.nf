nextflow.enable.dsl=2

process gtSort {
conda "bioconda::genometools-genometools conda-forge:coreutils"
cache 'deep'
input:
 path gffInput
output:
 path "output/${gffInput}"

script:
"""
mkdir -p output
echo "Sorting with gt..."
cat ${gffInput} | gt gff3 -tidy -sort -retainids -checkids -fixregionboundaries > output/${gffInput}
## I don't have -addintrons, as is not aware of existing introns. It will duplicate them
## It also doesn't give an informative ID to the intron, or any ID at all. They only have Parent= attributes
echo "Done gt sorting."
"""
}

process agatAddIntrons {
conda "bioconda::agat conda-forge:coreutils" //Also available via Docker/Apptainer
cache 'deep'
input:
 path gffInput
output:
 path "output/${gffInput}"

script:
"""
mkdir -p output
cat ${gffInput} | awk '!/\tintron\t/' > nointron.${gffInput}
agat_sp_add_introns.pl --gff nointron.${gffInput} --out output/${gffInput} 
"""
}

process igvtoolsSort {
conda "bioconda::igvtools conda-forge:coreutils"
publishDir './results/doubleSort' , mode:'link'
input:
 path gtSorted
output:
 path "output/${gtSorted}"
script:
"""
mkdir -p output
echo "Sorting with igvtools..."
igvtools sort ${gtSorted} output/${gtSorted}
echo "Done sorting."
"""
}

process EVM_gff3_simple_validator {
conda "bioconda::evidencemodeler conda-forge:coreutils"
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
conda "bioconda::perl-gfacs conda-forge::coreutils"
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
 //agatAddIntrons(gff)
 gtSort(gff)
 igvtoolsSort(gtSort.out)
 //EVM_gff3_simple_validator(igvtoolsSort.out)
 //gfacs_validator(igvtoolsSort.out)
emit:
 igvtoolsSort.out
}

workflow {
 Channel.fromPath(params.gff).set{gffInput}
 doubleSort_wf(gffInput)
}
