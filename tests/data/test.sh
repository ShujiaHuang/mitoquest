../../bin/mitoquest caller -f ce.fa.gz -o tt range.bam 
../../bin/mitoquest caller -f ce.fa.gz -o tt -b bamfile.list
../../bin/mitoquest caller -t 8 -r CHROMOSOME_III -f ce.fa.gz -o tt range.bam -b bamfile.list
../../bin/mitoquest caller -t 4 -f ce.fa.gz -o tt -r CHROMOSOME_I:1000-2400 range.bam range.bam

../../bin/mitoquest caller -f ref_mt.fasta -o tt test.bam --filename-has-samplename

../../bin/mitoquest caller -t 4 -f chrM_rCRS.decoy.fa -r chrM -o tt 21200715BFF2A.cram 15100105TLL4A.cram 00115121204M17BFF2.cram
