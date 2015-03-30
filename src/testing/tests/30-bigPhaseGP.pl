#!/usr/bin/env perl
use strict;
use warnings;

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
    plan tests => 2;
}

# this test tests the GP field parsing
my $testTag = 'test30';

my $resultsDir = "results/$testTag";
make_path($resultsDir);
my $expectedHap = File::Spec->catfile( $sd, q/samples/, $testTag,
"$testTag.expected.chr20.consensus.STv1.2.13.C100.K100.20_20059716_20399169.first2Samp.vcf.impute.hap"
);
my $expectedLegend = File::Spec->catfile($sd, q/samples/, $testTag,
"$testTag.expected.chr20.consensus.STv1.2.13.C100.K100.20_20059716_20399169.first2Samp.vcf.impute.legend"
);

my $resultsName = "$testTag.hapfuse";
my $resultsHap =
  File::Spec->catfile( $resultsDir, $resultsName . q/.impute.hap/ );
my $resultsLegend =
  File::Spec->catfile( $resultsDir, $resultsName . q/.impute.legend/ );

my $results_vcf = File::Spec->catfile( $resultsDir, $resultsName . q/.vcf/ );
my $cmd =   "./hapfuse -Ov -waverage -tGT,GP -TGT,GP -o $results_vcf $sd/samples/$testTag/$testTag.chr20*.bin.vcf";
print "Call: $cmd\n";
system $cmd;

# pull haplotypes out of vcf
system
"vcftools --vcf $results_vcf --from-bp 20059716 --to-bp 20399169 --chr 20 --IMPUTE --out ${resultsDir}/$resultsName";

compare_ok( $expectedLegend, $resultsLegend, "extracted correct sites" );
compare_ok( $expectedHap,    $resultsHap,    "hapfuse phases" );


