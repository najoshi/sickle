PROGRAM_NAME = sickle
VERSION = 0.93
CC = gcc
CFLAGS = -Wall -pedantic -DVERSION=$(VERSION)
DEBUG = -g
OPT = -O3
ARCHIVE = $(PROGRAM_NAME)_$(VERSION)

# Mac OS X - Linux may need different linking
LDFLAGS = -lz

default: build

sliding.o: src/sliding.c src/kseq.h src/sickle.h
	$(CC) $(CFLAGS) $(OPT) -c $?

trim_single.o: src/trim_single.c src/sickle.h src/kseq.h
	$(CC) $(CFLAGS) $(OPT) -c $?

trim_paired.o: src/trim_paired.c src/sickle.h src/kseq.h
	$(CC) $(CFLAGS) $(OPT) -c $?

sickle.o: src/sickle.c src/sickle.h
	$(CC) $(CFLAGS) $(OPT) -c $?

clean:
	rm -rf *.o src/*.gch ./sickle

distclean: clean
	rm -rf *.tar.gz

dist:
	tar -zcf $(ARCHIVE).tar.gz src Makefile

build: sliding.o trim_single.o trim_paired.o sickle.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(OPT) $? -o sickle

debug:
	$(MAKE) build "CFLAGS=-Wall -pedantic -g -DDEBUG"

