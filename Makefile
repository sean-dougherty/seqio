#CXXFLAGS=-g -std=c++11 -Wall -Werror
CXXFLAGS=-O2 -std=c++11 -Wall -Werror
shared_flags=-fPIC -shared
includes=-I src -I src/seqio

src_seqio=$(shell ls src/seqio/*.cpp src/util.cpp)
inc_seqio=$(shell ls src/seqio/*.h)
target_seqio=bld/lib/libseqio.so

src_pna=$(shell ls src/tools/pna/*.cpp src/util.cpp)
target_pna=bld/bin/pna

src_fasta=$(shell ls src/tools/fasta/*.cpp)
inc_fasta=src/tools/fasta/kseq.h
target_fasta=bld/bin/fasta

public_inc=$(shell ls src/seqio/*.h)

all: $(target_seqio) $(target_pna) $(target_fasta)

$(target_seqio): $(src_seqio) $(inc_seqio) Makefile
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(src_seqio) $(includes) $(shared_flags) -o $@ -lz -lrt

$(target_pna): $(src_pna) $(inc_seqio) $(target_seqio) Makefile
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(src_pna) $(includes) -o $@ -lseqio -L bld/lib -lrt

$(target_fasta): $(src_fasta) $(inc_seqio) $(target_seqio) Makefile
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(src_fasta) $(includes) -o $@ -lseqio -L bld/lib -lz -lrt -lboost_filesystem -lboost_system

install:
	cp $(target_seqio) /usr/local/lib
	cp $(target_pna) /usr/local/bin
	cp $(target_fasta) /usr/local/bin
	cp $(public_inc) /usr/local/include

clean:
	rm -rf bld