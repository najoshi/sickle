#!/bin/bash

if which valgrind &> /dev/null ; then
  FILE=$(mktemp)

  MEMCHECK_OUT=$(mktemp)

  trap "rm -f $FILE $MEMCHECK_OUT" EXIT

  valgrind --tool=memcheck --leak-check=full \
    ../src/sickle se \
      -t illumina \
      -f "$(dirname "$0")"/test.fastq \
      -o $FILE \
        2> $MEMCHECK_OUT

  cat $MEMCHECK_OUT
  grep -q 'ERROR SUMMARY: 0' $MEMCHECK_OUT
else
  echo "no valgrind, skipping test"
  exit 77
fi
