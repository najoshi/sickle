use strict;
use warnings;

package Sickle::Test;
use base 'Exporter';

use File::Basename qw< dirname >;
use Capture::Tiny qw< capture >;
use Test::More;

our @EXPORT = qw( sickle sickle_ok );
our $SICKLE = dirname(__FILE__) . "/../../../sickle";

# Linux libc and glibc tunable malloc(); on detected heap corruption, this will
# cause a message to stderr and immediate abort.  It's a poor-man's valgrind
# that you can enable always.  See malloc(3) for more details.
$ENV{MALLOC_CHECK_} = 2;

sub sickle {
    my ($command, @opts) = @_;
    my ($stdout, $stderr, $exit) = capture {
        system $SICKLE, $command, @opts
    };
    return ($stdout, $stderr, $exit);
}

sub sickle_ok {
    local $Test::Builder::Level = $Test::Builder::Level + 1;
    my ($stdout, $stderr, $exit) = sickle(@_);
    ok $exit >> 8 == 0, "exit code"
        or diag "error running command: " . join(" ", "sickle", @_) . "\n"
               ."exit code = " . ($exit >> 8) . " (system() return value = $exit)\n"
               ."stderr = \n$stderr";
    return ($stdout, $stderr);
}

1;
