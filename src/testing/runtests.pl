#!/usr/bin/env perl
use strict;
use warnings;
use TAP::Harness;
use File::Path qw(make_path remove_tree);

my $harness = TAP::Harness->new();

my $resultsDir = "results";
remove_tree($resultsDir) if -d $resultsDir;
make_path($resultsDir);

my $tests = q/find / . "$ENV{srcdir}/src/testing/tests" . q/ -name "*.pl"/;
$tests = qx/$tests/;
print "Found tests: $tests\n";
chomp $tests;
my @tests = split( /\s+/, $tests );

my $ret = $harness->runtests( sort @tests );

exit !$ret->all_passed;
