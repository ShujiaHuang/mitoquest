## Overview:

`databases`: Files and their sources include:
- dbSNP.chrM.vcf.gz (GRCh38, Build157) [from the dbSNP website](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz) and only extract all chrM dbSNP.
- gnomAD dataset [from the gnomAD browser](https://gnomad.broadinstitute.org/downloads). 
- HelixMTdb dataset [from the Helix website](https://www.helix.com/pages/mitochondrial-variant-database).
- MITOMAP dataset [from MITOMAP](https://www.mitomap.org/MITOMAP/resources).
- Disease-associated variants [from MITOMAP](https://www.mitomap.org/MITOMAP/resources).
- Disease-associated variants [from ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/).
- Haplogroup-associated variants extracted [from PhyloTree](https://www.phylotree.org/).
- List of loci in the mitochondrial genome [from MITOMAP](https://www.mitomap.org/foswiki/bin/view/MITOMAP/GenomeLoci).

`insilicos`: 
- PhyloP conservation scores from alignment of 100 vertebrates, [from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.100way.phyloP100way/).
- APOGEE prediction scores [from MitImpact](https://mitimpact.css-mendel.it/).
- MitoTip prediction scores [from MITOMAP](https://www.mitomap.org/MITOMAP/MitoTipScores).
- HmtVar disease scores retrieved [from HmtVar](https://www.hmtvar.uniba.it/).

`other_annotations`:
- UniProt Knowledgebase protein annotations obtained [from the UniProt FTP site](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/). Note annotations for mtDNA-encoded proteins were extracted from all available bed files for Homo sapiens (available in UniProt directory *UP000005640_9606_beds*).
- Residues involved in complex I proton transfer curated from [PMID:32972993](https://pubmed.ncbi.nlm.nih.gov/32972993/).
- tRNA position numbers, source file obtained [from MitoTIP](https://github.com/sonneysa/MitoTIP/) and manually converted to a list of tRNA positions for each mtDNA coordinate.
- RNA modifications and tRNA domains, obtained as previously described [by Lake et al](https://academic.oup.com/bioinformatics/article/38/10/2967/6567356).
- RNA base type extracted from manually curated data reported previously [by Lake et al](https://academic.oup.com/bioinformatics/article/38/10/2967/6567356), using file `all_RNA_bases.tsv`.
- Bases involved in rRNA:rRNA bridges curated from [PMID:25838379](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4501431/).
- Alignment of mitochondrial reference sequences for human and chimpanzee, obtained using reference sequences from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/). 
- Data on mitochondrial genome variation in the UK Biobank, obtained from [Hong, Battle et al](https://www.nature.com/articles/s41467-023-41785-7).

`synthetic_vcf`: Files in this directory obtained as described [previously](https://github.com/broadinstitute/gnomad-mitochondria/tree/main/gnomad_mitochondria/manuscript_analyses).

*The following are for nuclear genes:*
- Nuclear gene missense constraint metrics obtained from [Karczewski et al](https://www.nature.com/articles/s41586-020-2308-7).
- Essential genes obtained [from IMPC](https://www.ebi.ac.uk/mi/impc/essential-genes-search/).
- Developmental disorder genes obtained [from DECIPHER](https://www.deciphergenomics.org/ddd/ddgenes).

