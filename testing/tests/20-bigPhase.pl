#!/bin/perl
use Test::More tests => 1;
use Test::Files;
use File::Spec;

SKIP: {
      skip "vcftools does not exist on system", 1 unless qx/which vcftools/ =~ m/.*vcftools$/;

      my $expected_file =
        File::Spec->catfile(qw/ .. samples test01 test.1Col.hapfuse.expected.vcf /);
      my $resultsDir = qw/results/;
      my $results_file =         File::Spec->catfile($resultsDir, q|test20.hapfuse.vcf|);
      system "../bin/hapfuse -o $results_file ../samples/test20/*.vcf";

      # pull haplotypes out of vcf
      system "vcftools --vcf $results_file --IMPUTE --out basename -s .vcf $results_file";
#      compare_ok( $expected_file, $results_file, "hapfuse phases" );
  };



