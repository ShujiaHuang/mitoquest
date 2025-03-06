../../bin/hmtk -f ce.fa.gz -o tt range.bam 
../../bin/hmtk -f ce.fa.gz -o tt range.bam -r CHROMOSOME_MtDNA
../../bin/hmtk -f ce.fa.gz -o tt -r CHROMOSOME_I:1000-2400 range.bam 
../../bin/hmtk -t 4 -f ce.fa.gz -o tt -r CHROMOSOME_I:1000-2400 range.bam range.bam
