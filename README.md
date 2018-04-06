## *Photinus pyralis* geneset
This is the Official Geneset for [*Photinus pyralis* (Big Dipper Firefly)](https://en.wikipedia.org/wiki/Photinus_pyralis)

Genome version: Ppyr1.3  
Geneset version: PPYR_OGS1.0

This GIT repository is meant to keep track of any changes, and hopefully have a useful commit log which describes what changes occured.

### Dependencies
(Executable from your local command line)

* [gt](http://genometools.org/index.html)
* [igvtools](https://software.broadinstitute.org/software/igv/download) (also available in [IGV](https://software.broadinstitute.org/software/igv/home) GUI)
* [seqkit](https://github.com/shenwei356/seqkit)
* [gffutils](http://daler.github.io/gffutils/installation.html)

### Reporting gene model problems

Report issues [here](issues)

### Making direct gene model or annotation changes

 1. Fork this repository
 2. Download and edit the [PPYR_OGS1.0.gff3](./PPYR_OGS1.0.gff3) file
 3. Sort the GFF3 file and regenerate dependent files (CDS, mRNA, peptide) files using the [sort-OGS-GFF.sh](sort-OGS-GFF.sh) script.
 4. Commit your changes back to your repository, with an informative message for the changes that were made
 5. Submit a pull request to this repository