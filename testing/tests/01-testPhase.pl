#!/bin/perl
use Test::More tests => 1;
use Test::Files;
use File::Spec;

my $expected_file =
  File::Spec->catfile(qw/ .. samples test01 test.1Col.hapfuse.expected.vcf /);
my $results_file = File::Spec->catfile(qw| results test01.1Col.hapfuse.vcf |);

system "../bin/hapfuse -o $results_file ../samples/test01/test.1ColPart*";

compare_ok( $expected_file, $results_file, "hapfuse phases" );
