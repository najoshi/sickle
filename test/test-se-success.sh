#!/bin/bash

FILE=$(mktemp)

trap "rm -f $FILE" EXIT

sickle se \
  -t illumina \
  -f "$(dirname "$0")"/test.fastq \
  -o $FILE
