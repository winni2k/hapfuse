#!/bin/perl
use strict;
use warnings;
use Test::More tests => 2;
use Test::Files;
use File::Spec;
use File::Path qw(make_path remove_tree);

my $tag           = "test10";
my $expected_file = File::Spec->catfile( qw/ .. samples /,
    $tag, $tag . '.expected.madeUpData.fused.vcf' );

my $resDir = File::Spec->catfile( qw| results |, $tag );
make_path($resDir) unless -d $resDir;

my $resultsName = "madeUpData";
my $results_vcf =
  File::Spec->catfile( $resDir, $tag . ".$resultsName.fused.vcf" );

system "../bin/hapfuse -o $results_vcf ../samples/$tag/$tag.madeUpData*.vcf";

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

compare_filter_ok( $results_vcf, $expected_file, \&removeHeaderLines,
    "hapfuse phases and gets APPs right" );

sub removeHeaderLines {
    my $line = shift;
    if ( $line =~ m/^##/ ) {
        return q//;
    }
    return $line;
}
