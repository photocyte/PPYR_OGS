## *Photinus pyralis* geneset
This is the Official Geneset (OGS) for [*Photinus pyralis* (Big Dipper Firefly)](https://en.wikipedia.org/wiki/Photinus_pyralis)

As reported in this preprint: [https://www.biorxiv.org/content/early/2018/02/25/237586](https://www.biorxiv.org/content/early/2018/02/25/237586)

Genome version: Ppyr1.3  
Geneset version: PPYR_OGS1.0

This GIT repository is meant to keep track of any changes to the OGS, and hopefully have a useful audit log which describes what changes occured.

Also see: [fireflybase.org](http://www.fireflybase.org)

### Dependencies
(Executable from your local command line)

* [gt](http://genometools.org/index.html)
* [igvtools](https://software.broadinstitute.org/software/igv/download) (also available in [IGV](https://software.broadinstitute.org/software/igv/home) GUI)
* [seqkit](https://github.com/shenwei356/seqkit)
* [gffutils](http://daler.github.io/gffutils/installation.html)

### Reporting gene model problems

Report issues [here](https://github.com/photocyte/PPYR_OGS/issues)

### Making direct gene model or annotation changes

 1. Fork this repository
 2. Download and edit the [PPYR_OGS*.gff3](./PPYR_OGS1.0.gff3) file
 3. Sort the GFF3 file and regenerate dependent files (CDS, mRNA, peptide) files using the [sort-and-regenerate-OGS-GFF.sh](sort-and-regenerate-OGS-GFF.sh) script.
 4. Commit your changes back to your repository, with an informative message for the changes that were made
 5. Submit a pull request to this repository
 6. If all looks good I'll merge the pull request

### Future
 
Get a proper relational database ala [Chado](http://gmod.org/wiki/Chado_-_Getting_Started), or a collaborative GUI editing interface, ala [Apollo](http://genomearchitect.github.io).
