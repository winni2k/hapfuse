#!/bin/perl
use strict;
use warnings;
use Test::More tests => 4;
use Test::Files;
use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Glob ':bsd_glob';
use File::Slurp;

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
my $cmd = "./hapfuse -Ov -w average -tGT,APP -TGT,GP,APP -Ov -o $results_vcf " . join( ' ', @chunkVcfs );
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

$cmd = "./hapfuse -Ov -w step -o $results_vcf_wtccc -h $inputHaps -s $inputSamps";
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

sub vcfComp {
    my $line = shift;
    $line = removeHeaderLines($line);
    $line = roundFloats($line);
    return $line;
}

sub removeHeaderLines {
    my $line = shift;
    if ( $line =~ m/^##/ ) {
        return q//;
    }
    return $line;
}

sub roundFloats {
    my $line = shift;
    $line =~ s/(\d+[\.\d]{0,3})\d*/$1/g;

    # wow never thought I'd need a look-ahead assertion...
    $line =~ s/([,:])(\d+)(?=[,:\t])/$1$2.00/g;
    return $line;
}
