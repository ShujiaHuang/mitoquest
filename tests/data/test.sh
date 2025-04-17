../../bin/mitoquest caller -f ce.fa.gz -o tt range.bam 
../../bin/mitoquest caller -f ce.fa.gz -o tt -b bamfile.list
../../bin/mitoquest caller -t 8 -r CHROMOSOME_III -f ce.fa.gz -o tt range.bam -b bamfile.list
../../bin/mitoquest caller -t 4 -f ce.fa.gz -o tt -r CHROMOSOME_I:1000-2400 range.bam range.cram

../../bin/mitoquest caller -f ref_mt.fasta -o tt test.bam --filename-has-samplename

../../bin/mitoquest caller -t 4 -f chrM_rCRS.decoy.fa.gz -r chrM -o tt smp1.cram smp2.cram smp3.cram

../../bin/mitoquest caller -t 4 -f chrM_rCRS.decoy.fa.gz -r chrM:5745-5747 -o t smp4.cram

## Annotation
python ../../tools/mito_annotate.py -d ~/Projects/mitoquest/data -i tt_bak -o t.ann.txt -v t.ann.vcf

