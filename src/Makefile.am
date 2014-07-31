AM_CFLAGS = -Wall -pedantic -O3

bin_PROGRAMS = sickle

sickle_SOURCES = kseq.h print_record.c print_record.h sickle.c sickle.h sliding.c trim_paired.c trim_single.c

sickle_CPPFLAGS = $(ZLIB_CPPFLAGS)
sickle_LDADD = $(ZLIB_LIBS)
