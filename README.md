# sickle - A windowed adaptive trimming tool for FASTQ files using quality

## About

Most modern sequencing technologies produce reads that have
deteriorating quality towards the 3'-end and some towards the 5'-end
as well. Incorrectly called bases in both regions negatively impact
assembles, mapping, and downstream bioinformatics analyses.

Sickle is a tool that uses sliding windows along with quality and
length thresholds to determine when quality is sufficiently low to
trim the 3'-end of reads and also determines when the quality is
sufficiently high enough to trim the 5'-end of reads.  It will also
discard reads based upon the length threshold.  It takes the quality
values and slides a window across them whose length is 0.1 times the
length of the read.  If this length is less than 1, then the window is
set to be equal to the length of the read.  Otherwise, the window
slides along the quality values until the average quality in the
window rises above the threshold, at which point the algorithm
determines where within the window the rise occurs and cuts the read
and quality there for the 5'-end cut.  Then when the average quality
in the window drops below the threshold, the algorithm determines
where in the window the drop occurs and cuts both the read and quality
strings there for the 3'-end cut.  However, if the length of the
remaining sequence is less than the minimum length threshold, then the
read is discarded entirely (or replaced with an "N" record). 5'-end 
trimming can be disabled.

Sickle supports three types of quality values: Illumina, Solexa, and
Sanger. Note that the Solexa quality setting is an approximation (the
actual conversion is a non-linear transformation). The end
approximation is close. Illumina quality refers to qualities encoded
with the CASAVA pipeline between versions 1.3 and 1.7.  Illumina
quality using CASAVA >= 1.8 is Sanger encoded.

Note that Sickle will remove the 2nd fastq record header (on the "+"
line) and replace it with simply a "+". This is the default format for
CASAVA >= 1.8.

Sickle also supports gzipped file inputs and optional gzipped outputs. By default,
Sickle will produce regular (i.e. not gzipped) output, regardless of the input.
Sickle also has an option to truncate reads with Ns at the first N position.

There is also a sickle.xml file included in the package that can be used to add sickle to your
local [Galaxy](http://galaxy.psu.edu/) server.

## Citation
Sickle doesn't have a paper, but you can cite it like this:

    Joshi NA, Fass JN. (2011). Sickle: A sliding-window, adaptive, quality-based trimming tool for FastQ files 
    (Version 1.33) [Software].  Available at https://github.com/najoshi/sickle.

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
reads: `sickle se` and `sickle pe`.

Running sickle by itself will print the help:

    sickle

Running sickle with either the "se" or "pe" commands will give help
specific to those commands:

    sickle se
    sickle pe

### Sickle Single End (`sickle se`)

`sickle se` takes an input fastq file and outputs a trimmed version of
that file.  It also has options to change the length and quality
thresholds for trimming, as well as disabling 5'-trimming and enabling
truncation of sequences with Ns.

#### Examples

    sickle se -f input_file.fastq -t illumina -o trimmed_output_file.fastq
    sickle se -f input_file.fastq -t illumina -o trimmed_output_file.fastq -q 33 -l 40
    sickle se -f input_file.fastq -t illumina -o trimmed_output_file.fastq -x -n
    sickle se -t sanger -g -f input_file.fastq -o trimmed_output_file.fastq.gz
    sickle se --fastq-file input_file.fastq --qual-type sanger --output-file trimmed_output_file.fastq

### Sickle Paired End (`sickle pe`)

`sickle pe` can operate with two types of input.  First, it can take
two paired-end files as input and outputs two trimmed paired-end files
as well as a "singles" file.  The second form starts with a single
combined input file of reads where you have already interleaved the
reads from the sequencer.  In this form, you also supply a single
output file name as well as a "singles" file.  The "singles" file
contains reads that passed filter in either the forward or reverse
direction, but not the other.  Finally, there is an option (-M) to only 
produce one interleaved output file where any reads that did not pass 
filter will be output as a FastQ record with a single "N" (whose quality 
value is the lowest possible based upon the quality type), thus 
preserving the paired nature of the data.  You can also change the length 
and quality thresholds for trimming, as well as disable 5'-trimming and 
enable truncation of sequences with Ns.

#### Examples

    sickle pe -f input_file1.fastq -r input_file2.fastq -t sanger \
    -o trimmed_output_file1.fastq -p trimmed_output_file2.fastq \
    -s trimmed_singles_file.fastq

    sickle pe -f input_file1.fastq -r input_file2.fastq -t sanger \
    -o trimmed_output_file1.fastq -p trimmed_output_file2.fastq \
    -s trimmed_singles_file.fastq -q 12 -l 15

    sickle pe -f input_file1.fastq -r input_file2.fastq -t sanger \
    -o trimmed_output_file1.fastq -p trimmed_output_file2.fastq \
    -s trimmed_singles_file.fastq -n

    sickle pe -c combo.fastq -t sanger -m combo_trimmed.fastq \
    -s trimmed_singles_file.fastq -n

    sickle pe -t sanger -g -f input_file1.fastq -r input_file2.fastq \
    -o trimmed_output_file1.fastq.gz -p trimmed_output_file2.fastq.gz \
    -s trimmed_singles_file.fastq.gz

    sickle pe -c combo.fastq -t sanger -M combo_trimmed_all.fastq

    sickle pe --pe-file1 input_file1.fastq --pe-file2 input_file2.fastq --qual-type sanger \
    --output-pe1 trimmed_output_file1.fastq --output-pe2 trimmed_output_file2.fastq \
    --output-single trimmed_singles_file.fastq
