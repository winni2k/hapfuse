#!/bin/perl
use Test::More tests => 1;
use Test::Files;
use File::Spec;
use File::Basename;
use File::Path qw(make_path remove_tree);

# this test tests APP field parsing
SKIP: {
    skip "vcftools does not exist on system", 1
      unless qx/which vcftools/ =~ m/.*vcftools$/;

    my $resultsDir = qw|results/test20|;
    make_path($resultsDir);
    my $expected_file =
      File::Spec->catfile(qw/ .. samples test20 expected.impute.hap /);
    my $resultsName = "test20.hapfuse";
    my $results_file =
      File::Spec->catfile( $resultsDir, $resultsName . q/.impute.hap/ );

    my $results_vcf =
      File::Spec->catfile( $resultsDir, $resultsName . q/.vcf/ );
    system "../bin/hapfuse -o $results_vcf ../samples/test20/*.bin.vcf";

    # pull haplotypes out of vcf
    system
      "vcftools --vcf $results_vcf --IMPUTE --out ${resultsDir}/$resultsName";
    compare_ok( $expected_file, $results_file, "hapfuse phases" );
}

