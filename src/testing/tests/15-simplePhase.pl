#!/bin/perl
use Test::More;
use Test::Files;
use File::Spec;
use File::Basename;
use File::Path qw(make_path remove_tree);

my $sd = $ENV{srcdir};
unless ( qx/which vcftools/ =~ m/.*vcftools$/ ) {
    plan skip_all => "vcftools does not exist on system";
}
else {
    plan tests => 1;
}

# this test tests the GP field parsing
my $testTag = 'test15';

my $resultsDir = "results/$testTag";
make_path($resultsDir);
my $expectedHap = File::Spec->catfile( $sd, q/samples/, $testTag,
"$testTag.expected.chr20_20000120_20153108.STv1.2.17.C500.B0.i0.m20.K100.t2.f.q0.05.h.s.bin.vcf.hap"
);

my $resultsName = "$testTag.hapfuse";
my $resultsHap =
  File::Spec->catfile( $resultsDir, $resultsName . q/.impute.hap/ );

my $results_vcf = File::Spec->catfile( $resultsDir, $resultsName . q/.vcf/ );
system
  "./hapfuse -o $results_vcf $sd/samples/$testTag/$testTag.chr20*.bin.vcf";

# pull haplotypes out of vcf
system "vcftools --vcf $results_vcf --IMPUTE --out ${resultsDir}/$resultsName";

compare_ok( $expectedHap, $resultsHap, "hapfuse phases" );

