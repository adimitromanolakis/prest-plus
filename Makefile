
# Makefile for prest
#
# configure your command line options here
# specifically if you have a recent g++ compiler
# add to OPT: -mtune=native -march=native


OPT=-O2 
CC=g++

# OPT += -mtune=native -march=native -ffast-math


default: prest convert-geno

prest: prest4.cpp
	$(CC) $(OPT) -o prest prest4.cpp

convert-geno: convert-geno.cpp
	$(CC) $(OPT) -o convert-geno convert-geno.cpp

clean:
	rm -f prest4 convert-geno *.o

tar:
	rm -rf prest-plus
	mkdir prest-plus/
	cp *.cpp *.h *.pdf* Makefile reltypes.* prest-plus/
	tar cvzf prest-plus.tar.gz prest-plus/
	rm -rf prest-plus/

#
# Depreceated
#
prest-ted: prest_ted_5_v14.c
	$(CC) $(OPT) -ggdb -o prest-ted prest_ted_5_v14.c

