#!/bin/bash

pwd

FILE_O=$(mktemp)
FILE_P=$(mktemp)
FILE_S=$(mktemp)

trap "rm -f $FILE_O $FILE_P $FILE_S" EXIT

sickle pe \
  -t illumina \
  -f "$(dirname "$0")"/test.f.fastq \
  -r "$(dirname "$0")"/test.r.fastq \
  -o $FILE_O \
  -p $FILE_P \
  -s $FILE_S
