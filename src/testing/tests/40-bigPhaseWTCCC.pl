#!/bin/perl
use strict;
use warnings;
use Test::More;
use Test::Files;
use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Glob ':bsd_glob';
use File::Slurp;

my $sd = $ENV{srcdir};

unless ( qx/which vcftools/ =~ m/.*vcftools$/ ) {
    plan skip_all => "vcftools does not exist on system";
}
else {
    plan tests => 4;
}

# this test tests the GP field parsing
my $testTag = 'test40';

my $resultsDir = "results/$testTag";
make_path($resultsDir);
my $expectedHap =
  File::Spec->catfile( $sd, q/samples/, $testTag, "$testTag.ligated.i2.flip.hap" );
my $expectedLegend =
  File::Spec->catfile( $sd, q/samples/, $testTag,
    "$testTag.ligated.i2.legend" );

my $resultsName = "$testTag.hapfuse";
my $resultsHap =
  File::Spec->catfile( $resultsDir, $resultsName . q/.impute.hap/ );
my $resultsLegend =
  File::Spec->catfile( $resultsDir, $resultsName . q/.impute.modified.legend/ );

# define input hap files file
my $inputHaps =
  File::Spec->catfile( $resultsDir,
    $testTag . ".$resultsName.WTCCC.inputHaps" );
my $inputSamps = $inputHaps;
$inputSamps =~ s/Haps$/Samps/;

my @chunkHaps = bsd_glob(
"../samples/$testTag/$testTag.chr20_20*.STv1.2.13.C100.K100.first2Samp.bin.vcf.hap"
);
my @chunkSamps = map { my $s = $_; $s =~ s/\.hap$/.sample/; $s } @chunkHaps;

write_file( $inputHaps,  join( "\n", @chunkHaps ) );
write_file( $inputSamps, join( "\n", @chunkSamps ) );

my $results_vcf = File::Spec->catfile( $resultsDir, $resultsName . q/.vcf/ );
system "./hapfuse -w step -o $results_vcf -h $inputHaps -s $inputSamps";

# pull haplotypes out of vcf
system "vcftools --vcf $results_vcf --IMPUTE --out ${resultsDir}/$resultsName"
  . q( && perl -pne 's/^20-\S+/20/; s/allele/a/g; s/pos/position/; s/ID/id/' )
  . " < ${resultsDir}/$resultsName.impute.legend > $resultsLegend";

compare_ok( $expectedLegend, $resultsLegend, "extracted correct sites" );
compare_ok( $expectedHap,    $resultsHap,    "hapfuse phases" );


# and run on reversed chunks
my $inputHaps_rev = "$inputHaps.rev";
write_file( $inputHaps_rev,  join( "\n",reverse @chunkHaps ) );

my $results_vcf_rev = File::Spec->catfile( $resultsDir, $resultsName . q/.rev.vcf/ );
system "./hapfuse -w step -o $results_vcf_rev -h $inputHaps_rev -s $inputSamps";

# pull haplotypes out of vcf
system "vcftools --vcf $results_vcf_rev --IMPUTE --out ${resultsDir}/$resultsName.rev"
  . q( && perl -pne 's/^20-\S+/20/; s/allele/a/g; s/pos/position/; s/ID/id/' )
  . " < ${resultsDir}/$resultsName.rev.impute.legend > $resultsLegend.rev";

my $resultsHap_rev = $resultsHap;
$resultsHap_rev =~ s/\.impute\.hap$/.rev.impute.hap/;
compare_ok( $expectedLegend, "$resultsLegend.rev", "extracted correct sites" );
compare_ok( $expectedHap,    $resultsHap_rev,    "hapfuse phases" );

