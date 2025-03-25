../../bin/mitoquest caller -f ce.fa.gz -o tt range.bam 
../../bin/mitoquest caller -f ce.fa.gz -o tt -b bamfile.list
../../bin/mitoquest caller -t 8 -r CHROMOSOME_III -f ce.fa.gz -o tt range.bam -b bamfile.list
../../bin/mitoquest caller -t 4 -f ce.fa.gz -o tt -r CHROMOSOME_I:1000-2400 range.bam range.cram

../../bin/mitoquest caller -f ref_mt.fasta -o tt test.bam --filename-has-samplename

../../bin/mitoquest caller -t 4 -f chrM_rCRS.decoy.fa -r chrM -o tt smp1.cram smp2.cram smp3.cram

../../bin/mitoquest caller -t 4 -f chrM_rCRS.decoy.fa -r chrM:5745-5747 -o t 19106769BFF2B.sorted.markdup.realigner.cram > log

