# sickle - A windowed adaptive trimming tool for FASTQ files using quality

## About

Most modern sequencing technologies produce reads that have
deteriorating quality towards the 3'-end. Incorrectly called bases
here negatively impact assembles, mapping, and downstream
bioinformatics analyses.

Sickle is a tool that uses sliding windows to determine when quality
is sufficiently low to trim the 3'-end of reads.

## Requirements 

Sickle requires a C compiler; GCC or clang are recommended. Sickle
relies on Heng Li's kseq.h, which is bundled with the source.

Sickle also requires Zlib, which can be obtained at
<http://www.zlib.net/>.

## Building and Installing Sickle

To build Sickle, enter:

    make

Then, copy or move "sickle" to a directory in your $PATH.

## Usage

Sickle has two modes to work with both paired-end and single-end
reads: `sickle pe` and `sickle se`.


