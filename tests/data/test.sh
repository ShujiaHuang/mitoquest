../../bin/hmtk -f ce.fa.gz -o tt range.bam 
../../bin/hmtk -f ce.fa.gz -o tt -b bamfile.list
../../bin/hmtk -t 8 -r CHROMOSOME_III -f ce.fa.gz -o tt range.bam -b bamfile.list
../../bin/hmtk -f ce.fa.gz -o tt range.bam -r CHROMOSOME_MtDNA
../../bin/hmtk -f ce.fa.gz -o tt -r CHROMOSOME_I:1000-2400 range.bam 
../../bin/hmtk -t 4 -f ce.fa.gz -o tt -r CHROMOSOME_I:1000-2400 range.bam range.bam

../../bin/hmtk -t 4 -f Homo_sapiens.GRCh38.chrM_rCRS.decoy.fa -r chrM -o tt 21200715BFF2A.cram 15100105TLL4A.cram 00115121204M17BFF2.cram


