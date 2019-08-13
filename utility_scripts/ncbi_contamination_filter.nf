params.contamination = "/Users/tim/Downloads/20190813_lateralis_contamination.txt"
params.gffs = ""
params.genome = ""

gffFiles = params.gffs.tokenize(",")

gffFiles_ch = Channel.fromPath(gffFiles).flatten()

process parse_contamination_file {
input:
  file contaminationFile from file(params.contamination)
output:
 file "exclude.txt" optional true into excluded_ch
 file "trim.txt" optional true
script:
"""
##Note the extra backslashes for nextflow
cat ${contaminationFile} |  grep -Poz '(?s)Exclude:.*?\\n\\n' | sed 's/Exclude://g' | sed '/^[[:space:]]*\$/d' | tail -n +2 > exclude.txt

cat ${contaminationFile} |  grep -Poz '(?s)Trim:.*?\\n\\n' | sed 's/Trim://g' | sed '/^[[:space:]]*\$/d' | tail -n +2 > trim.txt

if [[ -s exclude.txt ]]
then
echo "pass"
else
rm -f exclude.txt
fi

if [[ -s trim.txt ]]
then
echo "pass"
else
rm -f trim.txt
fi
"""
}

process excludeToGrepPatterns {
input:
 file excludedFile from excluded_ch
output:
 file "scaffolds_to_exclude-regex.txt" into excludeRegexFile_ch1, excludeRegexFile_ch2
script:
"""
cat ${excludedFile} | cut -f 1 | sed 's/^/^/g' | sed 's/\$/\t/g' > scaffolds_to_exclude-regex.txt ##Add ^ to start and literal tab to end of files, for accurate grepping
"""
}

gffFiles_ch.combine(excludeRegexFile_ch1).into{filterCmds_ch1; filterCmds_ch2}

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
validExitStatus 0,1
input:
  set file(gff),file(excludeRegexFile) from filterCmds_ch2 
output:
  file "${gff}_removed.gff" into removedGffs
tag "${gff}"
script:
"""
zless ${gff} | grep -f ${excludeRegexFile} > ${gff}_removed.gff
"""
}

process filterGenome {
echo true
publishDir "filteredGenome/",overwrite:true
conda "seqkit"
input:
 file erf from excludeRegexFile_ch2
 file genome from file(params.genome)
output:
 file "${genome}_filtered.fa"
 file "${genome}_removed.fa"
tag "${genome}"
script:
"""
##Actually, seqkit grep doesn't need the regex patterns
cat ${erf} | sed 's/^^//g' | sed 's/[[:space:]]\$//g' > scaffolds.txt
seqkit grep -v -f scaffolds.txt ${genome} > ${genome}_filtered.fa
seqkit grep -f scaffolds.txt ${genome} > ${genome}_removed.fa
seqkit stat ${genome} ${genome}_filtered.fa ${genome}_removed.fa
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
