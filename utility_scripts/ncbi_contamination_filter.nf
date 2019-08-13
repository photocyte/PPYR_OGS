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
 file "exclude_confirmed.txt" optional true into excluded_ch
 file "trim_confirmed.txt" optional true into trim_ch
script:
"""
##Note the extra backslashes for nextflow
cat ${contaminationFile} |  grep -Poz '(?s)Exclude:.*?\\n\\n' | sed 's/Exclude://g' | sed '/^[[:space:]]*\$/d' | tail -n +2 > exclude.txt || true
cat ${contaminationFile} |  grep -Poz '(?s)Trim:.*?\\n\\n' | sed 's/Trim://g' | sed '/^[[:space:]]*\$/d' | tail -n +2 > trim.txt || true

seqkit fx2tab -n -i ${genome} | tr -s "\t" | sed 's/^/\\^/g' > genome_record_regex.txt

##Note: grep and pipefail can screw things up
##If grep doesn't find a match, it returns exitcode 1
##which makes the script quit at that point.  So the rest of the script won't be executed

echo "Checking NCBI contamination file against genome..."
cat exclude.txt | grep -f genome_record_regex.txt > exclude_confirmed.txt || true
cat trim.txt | grep -f genome_record_regex.txt > trim_confirmed.txt || true
cat exclude.txt | grep -v -f genome_record_regex.txt > exclude_nomatch.txt || true
cat trim.txt | grep -v -f genome_record_regex.txt > trim_nomatch.txt || true

echo "Checking for zero size files..."
if [[ -s exclude_confirmed.txt ]]
then
echo "pass"
else
echo "exclude_confirmed.txt is size zero"
rm -f exclude_confirmed.txt
fi

if [[ -s trim_confirmed.txt ]]
then
echo "pass"
else
echo "trim_confirmed.txt is size zero"
rm -f trim_confirmed.txt
fi
"""
}

process ncbiExcludeToGrepPatterns {
input:
 file excludedFile from excluded_ch
output:
 file "scaffolds_to_exclude-regex.txt" into excludeRegexFile_ch1, excludeRegexFile_ch2
script:
"""
cat ${excludedFile} | cut -f 1 | sed 's/^/^/g' | sed 's/\$/\t/g' > scaffolds_to_exclude-regex.txt ##Add ^ to start and literal tab to end of files, for accurate grepping
"""
}

process ncbiTrimtoGff3 {
conda "gawk genometools"
input:
 file trim_txt from trim_ch
output:
 file "trim.gff3" into trimGff3_ch1, trimGff3_ch2
script:
"""
##Note, extra backslahes for nextflow
cat ${trim_txt} | gawk 'match(\$0, /(.+)\t[0-9]+\t([0-9]+)\\.\\.([0-9]+)/,a) { print a[1] "\tgawk" "\t" "trim_me" "\t" a[2] "\t" a[3] "\t.\t+\t.\t" "ID=" a[1] ":" a[2] "-" a[3]}' | gt gff3 -tidy -sort -retainids > trim.gff3
"""

}

gffFiles_ch1.combine(excludeRegexFile_ch1).into{filterCmds_ch1; filterCmds_ch2}
gffFiles_ch2.combine(trimGff3_ch1).set{trimCmds_ch1}

process intersectingGffs {
conda "bedtools"
input:
 set file(gff),file(trim) from trimCmds_ch1
output:
 file "${gff}_intersect.gff" into intersectingGffs_ch
script:
"""
bedtools intersect -a ${gff} -b ${trim} > ${gff}_intersect.gff
"""
}

process filterGffsKeep {
conda "genometools"
input:
  set file(gff),file(excludeRegexFile) from filterCmds_ch1 
output:
  file "${gff}_filtered.gff" into filteredGffs
tag "${gff}"
script:
"""
zless ${gff} | grep -v -f ${excludeRegexFile} | gt gff3 -tidy -sort -retainids -fixregionboundaries > ${gff}_filtered.gff
"""
}

process filterGffsRemove {
input:
  set file(gff),file(excludeRegexFile) from filterCmds_ch2 
output:
  file "${gff}_removed.gff" into removedGffs
tag "${gff}"
script:
"""
zless ${gff} | grep -f ${excludeRegexFile} > ${gff}_removed.gff || true
"""
}

process filterGenome {
echo true
publishDir "filteredGenome/",overwrite:true
conda "seqkit"
input:
 file erf from excludeRegexFile_ch2
 file genome from genomeFile_ch2
output:
 file "${genome}_without-exclusions.fa" into filteredGenome_ch
 file "${genome}_excluded.fa"
tag "${genome}"
script:
"""
##Actually, seqkit grep doesn't need the regex patterns
cat ${erf} | sed 's/^^//g' | sed 's/[[:space:]]\$//g' > scaffolds.txt
seqkit grep -v -f scaffolds.txt ${genome} > ${genome}_without-exclusions.fa
seqkit grep -f scaffolds.txt ${genome} > ${genome}_excluded.fa
seqkit stat ${genome} ${genome}_filtered.fa ${genome}_excluded.fa
"""
}

process trimGenome {
publishDir "trimmedGenome/",overwrite:true
input:
 file genome from genomeFile_ch3
 file gff from trimGff3_ch2
output:
 file "${genome}_trimmed.fa" into trimmedGenome

script:
"""
bedtools maskfasta -fi ${genome} -bed ${gff} -fo ${genome}_tmp.fa
##Some of the gaps might show up at the beginning of the sequence, so have to trim them.  
seqkit replace --by-seq -p "^[Nn]+" -r "" ${genome}_tmp.fa | seqkit replace --by-seq -p "[Nn]+\$" -r "" > ${genome}_trimmed.fa
seqkit stat ./*.fa
"""
}

process compareOriginalTrimmedGenomes {
publishDir "trimmedGenome/",overwrite:true
input:
 file trim from trimmedGenome
 file genome from genomeFile_ch4
output:
 file "${genome}_only-diff.fa"
 file "${trim}_only-diff.fa"
script:
"""
seqkit rmdup -s -D dup_seqs.txt ${genome} ${trim} --ignore-case --md5 > /dev/null
cat dup_seqs.txt | cut -f 2 | cut -f 1 -d "," > scaffolds.txt
seqkit grep -v -f scaffolds.txt ${genome} > ${genome}_only-diff.fa
seqkit grep -v -f scaffolds.txt ${trim} > ${trim}_only-diff.fa
"""
}

process publishIntersectingGffs {
publishDir "intersectingGffs/",mode:'copy',overwrite:true
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

process publishRemovedGffs {
publishDir "removedGffs/",mode:'copy',overwrite:true
echo true
input:
 file rGffs from removedGffs.collect()
output:
 file rGffs
script:
"""
wc -l ./*.gff
"""
}

process publishFilteredGffs {
publishDir "filteredGffs/",overwrite:true
echo true
input:
 file fGffs from filteredGffs.collect()
output:
 file fGffs
script:
"""
wc -l ./*.gff
"""
}
