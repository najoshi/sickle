#!/bin/bash

FILE_O=$(mktemp)
FILE_P=$(mktemp)
FILE_S=$(mktemp)

trap "rm -f $FILE_O $FILE_P $FILE_S" EXIT

../src/sickle pe \
  -t illumina \
  -f "$(dirname "$0")"/segfault-on-empty-reads-1.fastq \
  -r "$(dirname "$0")"/segfault-on-empty-reads-2.fastq \
  -o $FILE_O \
  -p $FILE_P \
  -s $FILE_S
