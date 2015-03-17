#!/bin/bash

if which valgrind &> /dev/null ; then
  FILE_O=$(mktemp)
  FILE_P=$(mktemp)
  FILE_S=$(mktemp)

  MEMCHECK_OUT=$(mktemp)

  trap "rm -f $FILE_O $FILE_P $FILE_S $MEMCHECK_OUT" EXIT

  valgrind --tool=memcheck --leak-check=full \
    ../src/sickle pe \
      -t illumina \
      -f "$(dirname "$0")"/test.f.fastq \
      -r "$(dirname "$0")"/test.r.fastq \
      -o $FILE_O \
      -p $FILE_P \
      -s $FILE_S \
        2> $MEMCHECK_OUT

  cat $MEMCHECK_OUT
  grep -q 'ERROR SUMMARY: 0' $MEMCHECK_OUT
else
  echo "no valgrind, skipping test"
  exit 77
fi
