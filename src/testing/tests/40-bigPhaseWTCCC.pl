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
    plan tests => 9;
}

# this test tests the GP field parsing
my $testTag = 'test40';

my $resultsDir = "results/$testTag";
make_path($resultsDir);
my $expectedHap =
  File::Spec->catfile( $sd, q/samples/, $testTag,
    "$testTag.ligated.i2.flip.hap" );
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
system "./hapfuse -Ov -w step -o $results_vcf -h $inputHaps -s $inputSamps";

# pull haplotypes out of vcf
system "vcftools --vcf $results_vcf --IMPUTE --out ${resultsDir}/$resultsName"
  . q( && perl -pne 's/^20-\S+/20/; s/allele/a/g; s/pos/position/; s/ID/id/' )
  . " < ${resultsDir}/$resultsName.impute.legend > $resultsLegend";

compare_ok( $expectedLegend, $resultsLegend, "extracted correct sites" );
compare_ok( $expectedHap,    $resultsHap,    "hapfuse phases" );

# and run on reversed chunks
my $inputHaps_rev = "$inputHaps.rev";
write_file( $inputHaps_rev, join( "\n", reverse @chunkHaps ) );

my $results_vcf_rev =
  File::Spec->catfile( $resultsDir, $resultsName . q/.rev.vcf/ );
system
  "./hapfuse -Ov -w step -o $results_vcf_rev -h $inputHaps_rev -s $inputSamps";

# pull haplotypes out of vcf
system
"vcftools --vcf $results_vcf_rev --IMPUTE --out ${resultsDir}/$resultsName.rev"
  . q( && perl -pne 's/^20-\S+/20/; s/allele/a/g; s/pos/position/; s/ID/id/' )
  . " < ${resultsDir}/$resultsName.rev.impute.legend > $resultsLegend.rev";

my $resultsHap_rev = $resultsHap;
$resultsHap_rev =~ s/\.impute\.hap$/.rev.impute.hap/;
compare_ok( $expectedLegend, "$resultsLegend.rev", "extracted correct sites" );
compare_ok( $expectedHap,    $resultsHap_rev,      "hapfuse phases" );

### now try the same, but output as WTCCC haps/sample files
###
my $tag = "wtccc";
my $resultsHap_wtccc =
  File::Spec->catfile( $resultsDir, $resultsName . ".$tag.haps" );
my $resultsSample_wtccc =
  File::Spec->catfile( $resultsDir, $resultsName . ".$tag.sample" );

my $cmd =
"./hapfuse -Ow -w step -o $resultsHap_wtccc,$resultsSample_wtccc -h $inputHaps -s $inputSamps";
print "Call: " . $cmd . "\n";
system $cmd;

# convert haps to haps/sample
my $resLeg_wtccc =
  File::Spec->catfile( $resultsDir, $resultsName . ".$tag.i2.legend" );
my $resHap_wtccc =
  File::Spec->catfile( $resultsDir, $resultsName . ".$tag.i2.hap" );
$cmd = "(echo 'id position a0 a1'; cut -f1,3-5 -d' ' $resultsHap_wtccc)"
  . " > $resLeg_wtccc ";
print $cmd;
system $cmd;
$cmd = "cut -f6- -d' ' $resultsHap_wtccc > $resHap_wtccc";
print $cmd;
system $cmd;
compare_ok( $expectedLegend, $resLeg_wtccc,
    "output as wtccc haps -- i2 legend ok" );
compare_ok( $expectedHap, $resHap_wtccc, "output as wtccc haps -- i2 hap ok" );
compare_ok( $chunkSamps[0], $resultsSample_wtccc,
    "output as wtccc sample -- ok" );

### now try and output to bcf, but only GT field
# use test30 bcfs
my @chunkVCFs = reverse bsd_glob(
    "../samples/test30/test30.chr20_20*.STv1.2.13.C100.K100.first2Samp.bin.vcf"
);

$tag = q/from_VCF_only_GT/;
$results_vcf = File::Spec->catfile( $resultsDir, $resultsName . ".$tag.vcf" );
$cmd = "./hapfuse -Ov -w step -tGT -o $results_vcf " . join( " ", @chunkVCFs );
print $cmd. "\n";
system $cmd;

# pull haplotypes out of vcf
system
  "vcftools --vcf $results_vcf --IMPUTE --out ${resultsDir}/$resultsName.$tag"
  . q( && perl -pne 's/^20-\S+/20/; s/allele/a/g; s/pos/position/; s/ID/id/' )
  . " < ${resultsDir}/$resultsName.$tag.impute.legend > $resultsLegend.$tag";
compare_ok( $expectedLegend, "$resultsLegend.$tag",
    "using only GT gives ok haps" );

my $resultsHap_gt = $resultsHap;
$resultsHap_gt =~ s/\.impute\.hap/.$tag.impute.hap/;
compare_ok( $expectedHap, $resultsHap_gt, "hapfuse with output GT phases" );

