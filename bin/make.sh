#g++ -O3 -fPIC ../scripts/count_align_fragments.cpp ../htslib/libhts.a -I ../htslib -lz -lm -lpthread -lcurl -llzma -lbz2 -lssl -lcrypto -o count_align_fragments
#g++ -o mtvariantcaller -O3 -fPIC ../src/*.cpp ../src/io/*.cpp ../htslib/libhts.a -I ../htslib -I ../src -lz -lbz2 -lm -llzma -lpthread -lcurl -lssl -lcrypto

g++ -o count_align_fragments ../scripts/count_align_fragments.cpp -lhts -lz -lcurl -pthread -O3 -fPIC
g++ -o mitoquest -O3 -fPIC ../src/*.cpp ../src/io/*.cpp -I ../src -lhts -lz -lbz2 -lm -llzma -lpthread -lcurl -lssl -lcrypto

