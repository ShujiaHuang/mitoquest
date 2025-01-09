#g++ -o count_align_fragments count_align_fragments.cpp -lhts -lz -lcurl -pthread -O3
g++ -O3 -fPIC ../scripts/count_align_fragments.cpp ../htslib/libhts.a -I ../htslib -lz -lm -lpthread -lcurl -llzma -lbz2 -o count_align_fragments
