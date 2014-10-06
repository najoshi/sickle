use strict;
use warnings;
use 5.014;
use Test::More;
use Sickle::Test;

my $fastq = $0 =~ s/\.t$/.fastq/r;

my @se = (qw(se -t sanger -l 1 -w 1 --quiet -f), $fastq, qw(-o /dev/stdout));

my ($output, $err);

($output, $err) = sickle_ok(@se);
is $output, <<'', "default; no -n or -N";
@middle
ABCDEFGHIJKLMNOPQRSTUVWXYZ
+
~~~~~~~~~~~~~~~~~~~~~~~~~~
@start_no_middle
ANBCDEFGHIJKLMOPQRSTUVWXYZ
+
~~~~~~~~~~~~~~~~~~~~~~~~~~
@start_trimmed_and_middle
CDEFGHIJKLMNOPQRSTUVW
+
~~~~~~~~~~~~~~~~~~~~~
@start_trimmed_no_middle_end
CDEFGHIJKLMOPQRSTUVWX
+
~~~~~~~~~~~~~~~~~~~~~

($output, $err) = sickle_ok(@se, "-n");
is $output, <<'', "with -n";
@middle
ABCDEFGHIJKLM
+
~~~~~~~~~~~~~
@start_no_middle
A
+
~
@start_trimmed_and_middle
CDEFGHIJKLM
+
~~~~~~~~~~~
@start_trimmed_no_middle_end
CDEFGHIJKLMOPQRSTUVWX
+
~~~~~~~~~~~~~~~~~~~~~

($output, $err) = sickle_ok(@se, "-N");
is $output, <<'', "with -N";
@start_trimmed_no_middle_end
CDEFGHIJKLMOPQRSTUVWX
+
~~~~~~~~~~~~~~~~~~~~~

done_testing;
