params.contamination = "/Users/tim/Downloads/20190813_lateralis_contamination.txt"
params.gffs = ""
params.genome = ""

gffFiles = params.gffs.tokenize(",")

gffFiles_ch = Channel.fromPath(gffFiles).flatten()
gffFiles_ch.into{ gffFiles_ch1 ; gffFiles_ch2 }

genomeFile = Channel.fromPath(params.genome)
genomeFile.into{ genomeFile_ch1 ; genomeFile_ch2 ; genomeFile_ch3 ; genomeFile_ch4 }

process parse_contamination_file {
input:
  file contaminationFile from file(params.contamination)
  file genome from genomeFile_ch1
output:
 file "exclude.txt" optional true into excluded_ch
 file "trim_notsplit.txt" optional true into trim_reformat_ch
tag "${contaminationFile}"
script:
"""
##Note the extra backslashes for nextflow
cat ${contaminationFile} |  grep -Poz '(?s)Exclude:.*?\\n\\n' | sed 's/Exclude://g' | sed '/^[[:space:]]*\$/d' | tail -n +2 > exclude.txt
cat ${contaminationFile} |  grep -Poz '(?s)Trim:.*?\\n\\n' | sed 's/Trim://g' | sed '/^[[:space:]]*\$/d' | tail -n +2 > trim_notsplit.txt


echo "Checking for zero size files..."
if [[ -s exclude.txt ]]
then
echo "pass"
else
echo "exclude.txt is size zero"
rm -f exclude.txt
fi

if [[ -s trim_notsplit.txt ]]
then
echo "pass"
else
echo "trim_notsplit.txt is size zero"
rm -f trim_notsplit.txt
fi
"""
}

process ncbiTrimToGff {
input:
 file trim from trim_reformat_ch
output:
 file "trim.gff" into trim_ch
script:
"""
#!/usr/bin/env python
import os
read_handle = open("${trim}","rU")
write_handle = open("trim.gff","w")

##E.g. Alat1.3_scaffold_2172	21117	1..167,1578..1707,16427..16460	Wolbachia endosymbiont
for line in read_handle.readlines():
    splitline = line.split("\t")
    regions = splitline[2].split(",")
    ##This double loop makes extra lines for the multiple regions in a single line
    for r in regions:
        coords = r.split("..")
        start_int = int(coords[0])
        end_int = int(coords[1])
        span_int = end_int - start_int
        span_str = str(span_int)
        write_handle.write("\t".join([ splitline[0], "python", "trim_me", coords[0] , coords[1] , span_str , "+" , "." , "ID="+splitline[0]+":"+coords[0]+"-"+coords[1]+";description="+splitline[3] ]))
read_handle.close()
write_handle.close()
"""
}

process gffSort {
publishDir "trim_output/intersectingGffs/",mode:'copy',overwrite:true
conda "gawk genometools"
input:
 file trim_txt from trim_ch
output:
 file "trim.gff3" into trimGff3_ch1, trimGff3_ch2
script:
"""
gt gff3 -tidy -sort -retainids ${trim_txt} | grep -v "###" > trim.gff3
"""

}

gffFiles_ch2.combine(trimGff3_ch1).set{trimCmds_ch1}

process intersectingGffs {
errorStrategy "finish"
conda "bedtools grep"
input:
 set file(gff),file(trim) from trimCmds_ch1
output:
 file "${gff}_intersect.gff" into intersectingGffs_ch
script:
"""
zless ${gff} | grep -Pv "^#" > tmp.gff
bedtools intersect -a tmp.gff -b ${trim} > ${gff}_intersect.gff
rm -f tmp.gff
"""
}

process trimGenome {
publishDir "trim_output",overwrite:true
conda "seqkit bedtools"
input:
 file genome from genomeFile_ch3
 file gff from trimGff3_ch2
output:
 file "${genome}_trimmed.fa" into trimmedGenome
 file "${genome}_subset-original.fa" into originalSubset, originalSubset_ch2
 file "${genome}_subset-trimmed.fa" into trimmedSubset
tag "${genome}"
script:
"""
cat ${gff} | grep -v "#" | cut -f 1 > scaffolds.txt
seqkit grep --delete-matched -f scaffolds.txt ${genome} > ${genome}_subset-original.fa
seqkit grep -v -f scaffolds.txt ${genome} > ${genome}_remainder.fa
bedtools maskfasta -fi ${genome}_subset-original.fa -bed ${gff} -fo ${genome}_masktmp.fa

##Some of the gaps might show up at the beginning of the sequence, so have to trim them.  
seqkit replace --by-seq -p "^[Nn]+" -r "" ${genome}_masktmp.fa | seqkit replace --by-seq -p "[Nn]+\$" -r "" > ${genome}_subset-trimmed.fa
seqkit sort -n ${genome}_subset-trimmed.fa ${genome}_remainder.fa > ${genome}_trimmed.fa

rm -f ${genome}_masktmp.fa
"""
}

process publishIntersectingGffs {
publishDir "trim_output/intersectingGffs/",mode:'copy',overwrite:true
echo true
input:
 file iGffs from intersectingGffs_ch.collect()
output:
 file iGffs
script:
"""
wc -l ./*.gff
"""
}

originalSubset.splitFasta(by:1,file:true).set{originalSplit}
trimmedSubset.splitFasta(by:1,file:true).set{trimmedSplit}

originalSplit.merge(trimmedSplit).set{ groupedFastas }

process prefilterGff {
input:
 file gff from gffFiles_ch1
 file fasta from originalSubset_ch2
output:
 file "target.${gff}.gff3" optional true into targetGff_ch1
 file "ignored.${gff}.gff3" into ignoredGff
 file "prelift.${gff}.gff3" into preliftGff
script:
"""
seqkit fx2tab --only-id -n ${fasta} | tr -s "\t" > target_scaffolds.txt
echo "##gff-version 3" >> target_scaffolds.txt
gt gff3 -tidy -sort -retainids -fixregionboundaries ${gff} > prelift.${gff}.gff3
grep -f target_scaffolds.txt prelift.${gff}.gff3 > target.${gff}.gff3
grep -v -f target_scaffolds.txt prelift.${gff}.gff3 > ignored.${gff}.gff3

if [[ \$(wc -l <target.${gff}.gff3) -le 1 ]]
then
    echo "No targets were found. Deleting target gff3 file."
    rm -f target.${gff}.gff3
fi
"""
}

targetGff_ch1.combine(groupedFastas).set{ liftoverCmds  }

process liftoverGffsPerScaffold {
errorStrategy "finish"
//maxForks 1
input:
 set file(gff),file(oldFa),file(newFa) from liftoverCmds
output:
 file "liftover_output/srt*.gff3" optional true into liftedGffs
tag "${newFa}"
script:
"""
nextflow run /Users/tim/Source/git/doSameSpeciesLiftOver_nextflow/doSameSpeciesLiftOver.nf --old ${oldFa} --new ${newFa} --gff ${gff}
rm -rf work
"""
}

process mergeLiftedGffs {
publishDir "trim_output",mode:"copy",overwrite:true
conda "genometools grep"
input:
 file lgff from liftedGffs.collectFile(name: 'merged_lifted_gffs.gff',skip:2)
 file igff from ignoredGff
 file pgff from preliftGff
output:
 file "merged_lifted.gff3"
 file pgff
script:
"""
cat ${lgff} ${igff} | grep -Pv "^#" | gt gff3 -tidy -sort -retainids > merged_lifted.gff3
"""
}
