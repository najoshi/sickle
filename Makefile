PROGRAM_NAME = sickle
VERSION = 1.33
CC = gcc
CFLAGS = -Wall -pedantic -DVERSION=$(VERSION)
DEBUG = -g
OPT = -O3
ARCHIVE = $(PROGRAM_NAME)_$(VERSION)
LDFLAGS=
LIBS = -lz
SDIR = src

.PHONY: clean default build distclean dist debug

default: build

sliding.o: $(SDIR)/sliding.c $(SDIR)/kseq.h $(SDIR)/sickle.h
	$(CC) $(CFLAGS) $(OPT) -c $(SDIR)/$*.c

trim_single.o: $(SDIR)/trim_single.c $(SDIR)/sickle.h $(SDIR)/kseq.h
	$(CC) $(CFLAGS) $(OPT) -c $(SDIR)/$*.c

trim_paired.o: $(SDIR)/trim_paired.c $(SDIR)/sickle.h $(SDIR)/kseq.h
	$(CC) $(CFLAGS) $(OPT) -c $(SDIR)/$*.c

sickle.o: $(SDIR)/sickle.c $(SDIR)/sickle.h
	$(CC) $(CFLAGS) $(OPT) -c $(SDIR)/$*.c

print_record.o: $(SDIR)/print_record.c $(SDIR)/print_record.h
	$(CC) $(CFLAGS) $(OPT) -c $(SDIR)/$*.c

clean:
	rm -rf *.o $(SDIR)/*.gch ./sickle

distclean: clean
	rm -rf *.tar.gz

dist:
	tar -zcf $(ARCHIVE).tar.gz src Makefile README.md sickle.xml LICENSE

build: sliding.o trim_single.o trim_paired.o sickle.o print_record.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(OPT) $? -o sickle $(LIBS)

debug:
	$(MAKE) build "CFLAGS=-Wall -pedantic -g -DDEBUG"

