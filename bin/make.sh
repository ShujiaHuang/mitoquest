g++ -O3 -fPIC ../scripts/count_align_fragments.cpp ../htslib/libhts.a -I ../htslib -lz -lm -lpthread -lcurl -llzma -lbz2 -o count_align_fragments
g++ -O3 -fPIC ../src/*.cpp ../src/io/*.cpp ../htslib/libhts.a -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o mtvariantcaller
