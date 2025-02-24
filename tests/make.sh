g++ -O3 -fPIC test_fasta.cpp ../src/io/fasta.cpp ../src/io/utils.cpp ../htslib/libhts.a -I ../src -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_fasta && ./test_fasta && rm -f test_fasta

g++ -O3 -fPIC test_bamheader.cpp ../src/io/bam_header.cpp ../src/io/utils.cpp ../htslib/libhts.a -I ../src -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_bamheader && ./test_bamheader && rm -f test_bamheader

g++ -O3 -fPIC test_bamrecord.cpp ../src/io/bam_record.cpp ../src/io/fasta.cpp ../src/io/bam_header.cpp ../src/io/utils.cpp ../htslib/libhts.a -I ../src -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_bamrecord && ./test_bamrecord && rm -f test_bamrecord

g++ -O3 -fPIC test_bam.cpp ../src/io/bam.cpp ../src/io/bam_header.cpp ../src/io/bam_record.cpp ../src/io/utils.cpp ../htslib/libhts.a -I ../src -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_bam && ./test_bam && rm -f test_bam

g++ -O3 -fPIC test_threadpool.cpp -I ../src -o test_threadpool && ./test_threadpool && rm -f test_threadpool

g++ -O3 -fPIC test_utils.cpp ../src/io/utils.cpp -I ../src -o test_utils && ./test_utils && rm -f test_utils

