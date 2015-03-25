#!/bin/perl
use strict;
use warnings;
use Test::More tests => 8;
use Test::Files;
use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Glob ':bsd_glob';
use File::Slurp;
use FindBin qw($Bin);
use lib "$Bin/..";
use TestHelpers qw(:ALL);

my $tag           = "test10";
my $sd            = $ENV{srcdir};
my $expected_file = File::Spec->catfile( $sd, q/samples/,
    $tag, $tag . '.expected.madeUpData.fused.vcf' );

my $resDir = File::Spec->catfile( qw| results |, $tag );
make_path($resDir) unless -d $resDir;

my $resultsName = "madeUpData";
my $results_vcf =
  File::Spec->catfile( $resDir, $tag . ".$resultsName.fused.vcf" );

my @chunkVcfs = bsd_glob("../samples/$tag/$tag.madeUpData*.vcf");
my $cmd = "./hapfuse -Ov -w average -tGT,APP -TGT,GP,APP -Ov -o $results_vcf "
  . join( ' ', @chunkVcfs );
print "Call: $cmd\n";
system $cmd;

SKIP: {
    skip "vcftools does not exist on system", 1
      unless qx/which vcftools/ =~ m/.*vcftools$/;

    # pull haplotypes out of vcf
    my $expectedHap =
      File::Spec->catfile( qw/.. samples/, $tag,
        $tag . ".expected.$resultsName.hap" );
    my $resultsHap =
      File::Spec->catfile( $resDir, $tag . ".$resultsName.impute.hap" );
    system
      "vcftools --vcf $results_vcf --IMPUTE --out ${resDir}/$tag.$resultsName";
    compare_ok( $resultsHap, $expectedHap, "hapfuse phases fuzzy APPs ok" );

}

compare_filter_ok( $results_vcf, $expected_file, \&vcfComp,
    "hapfuse phases and gets APPs right" );

# let's do the same by converting to WTCCC haps first
my @chunkHaps =
  map { my $s = $_; $s =~ s/\.vcf/.hap/; $s } @chunkVcfs;
my @chunkSamps = map { my $s = $_; $s =~ s/\.vcf/.sample/; $s } @chunkVcfs;
my $inputHaps =
  File::Spec->catfile( $resDir, $tag . ".$resultsName.WTCCC.inputHaps" );
my $inputSamps = $inputHaps;
$inputSamps =~ s/Haps$/Samps/;

write_file( $inputHaps,  join( "\n", @chunkHaps ) );
write_file( $inputSamps, join( "\n", @chunkSamps ) );

my $results_vcf_wtccc = $results_vcf;
$results_vcf_wtccc =~ s/\.vcf$/.wtccc.vcf/;

my $expected_wtccc_file = $expected_file;
$expected_wtccc_file =~ s/\.vcf$/.wtccc.vcf/;

$cmd =
  "./hapfuse -Ov -w step -o $results_vcf_wtccc -h $inputHaps -s $inputSamps";
print "Call: $cmd\n";
system $cmd;

compare_filter_ok( $results_vcf_wtccc, $expected_wtccc_file, \&vcfComp,
    "hapfuse phases from haplotypes" );

# Let's also test out of order chunk phasing works
write_file( $inputHaps, join( "\n", reverse @chunkHaps ) );

my $results_vcf_wtccc_rev = $results_vcf;
$results_vcf_wtccc_rev =~ s/\.vcf$/.wtccc.rev.vcf/;

$cmd =
"./hapfuse -Ov -w step -o $results_vcf_wtccc_rev -h $inputHaps -s $inputSamps";
print "Call: $cmd\n";
system $cmd;

compare_filter_ok( $results_vcf_wtccc_rev, $expected_wtccc_file, \&vcfComp,
    "hapfuse phases from haplotypes, chunks in reverse order" );

### Now redo, step analysis, but output as wtccc_haps file
my $expected_wtccc_haps_file = $expected_wtccc_file;
$expected_wtccc_haps_file =~ s/\.vcf/.wtccc.haps/;
my $expected_wtccc_sample_file = $expected_wtccc_haps_file;
$expected_wtccc_sample_file =~ s/\.haps/.sample/;
my $base = $results_vcf_wtccc_rev;
$base =~ s/\.vcf//;
my $results_haps_wtccc_rev   = "$base.haps";
my $results_sample_wtccc_rev = "$base.sample";
$cmd =
"./hapfuse -Ow -w step -o $results_haps_wtccc_rev,$results_sample_wtccc_rev -h $inputHaps -s $inputSamps";
print "Call: $cmd\n";
system $cmd;
compare_ok( $results_haps_wtccc_rev, $expected_wtccc_haps_file,
    "hapfuse phases from haplotypes to wtccc haps, chunks in reverse order" );
compare_ok( $results_sample_wtccc_rev, $expected_wtccc_sample_file,
    "hapfuse phases from haplotypes to wtccc haps, chunks in reverse order" );

### now redo using prefix argument instead
my $results_haps_wtccc_rev2   = "$base.2.hap";
my $results_sample_wtccc_rev2 = "$base.sample";
$cmd = "./hapfuse -Ow -w step -o $base.2 -h $inputHaps -s $inputSamps";
print "Call: $cmd\n";
system $cmd;
system "gunzip -f $base.2.hap.gz";
compare_ok( $results_haps_wtccc_rev2, $expected_wtccc_haps_file,
"hapfuse phases from haplotypes to wtccc haps, chunks in reverse order, prefix supplied"
);
compare_ok( $results_sample_wtccc_rev2, $expected_wtccc_sample_file,
"hapfuse phases from haplotypes to wtccc haps, chunks in reverse order, prefix supplied"
);
