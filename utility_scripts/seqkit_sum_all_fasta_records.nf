nextflow.enable.dsl=2

process do_md5 {
input:
 path(file)
output:
 env FOO
shell:
'''
FOO=$(openssl md5 !{file} | grep -oE '.{32}$')
'''
}

process seqkit_sum {
//publishDir 'results'
input:
 tuple path(fasta),path(filename),val(md5)
output:
 path "${fasta}.seqkit.sum.txt"
shell:
'''
RECORD_NAME=$(seqkit fx2tab !{fasta} | cut -f 1)
seqkit sum --basename --all !{fasta} | sed "s^!{fasta}^${RECORD_NAME}^g" | tr -d '\n' > !{fasta}.seqkit.sum.txt
echo "\t!{filename}\t!{md5}" >> !{fasta}.seqkit.sum.txt
'''
}

process publish {
publishDir 'results', mode: 'link', overwrite: true
input:
 path(report)
 val(md5)
 path(filename)
//If doing path(filename) rather than val(filename) it gets converted to a basename
output:
 path "out/md5_${md5}__${filename}__${report}"
shell:
'''
mkdir -p out
sort -nr -k 4,4 !{report} > out/md5_!{md5}__!{filename}__!{report}
'''
}

workflow {
fasta_ch = Channel.fromPath(params.fasta)

do_md5(fasta_ch)

file_md5_tuple = fasta_ch.combine(do_md5.out)

fastaChunks = fasta_ch.splitFasta(by:1,file:true)

fastaChunks.combine(file_md5_tuple).set{fastaTuples}

//fastaTuples.view()

seqkit_sum(fastaTuples).collectFile(name:'collected_seqkit_sum.txt',sort: true).set{ collected_report }

publish(collected_report,do_md5.out,fasta_ch)

}
