#!/usr/bin/env perl
use Test::More tests => 1;
use Test::Files;
use File::Spec;
use File::Basename;
use File::Path qw(make_path remove_tree);

my $sd = $ENV{srcdir};

# this test tests APP field parsing and fusing
# APP is quite a bit more accurate than using GP!
SKIP: {
    skip "vcftools does not exist on system", 1
      unless qx/which vcftools/ =~ m/.*vcftools$/;

    my $resultsDir = File::Spec->catdir(qw|results test20|);
    make_path($resultsDir);
    my $expected_file = File::Spec->catfile(
        $sd, qw/samples test20 test20.expected.chr20.consensus.STv1.2.13.C100.K100.20_20059716_20399169.first2Samp.secondIndHapSwapped.impute.editedForApp.hap /
    );
    my $resultsName = "test20.hapfuse";
    my $results_file =
      File::Spec->catfile( $resultsDir, $resultsName . q/.impute.hap/ );

    my $results_vcf =
      File::Spec->catfile( $resultsDir, $resultsName . q/.vcf/ );
    system "./hapfuse -Ov -w average -o $results_vcf $sd/samples/test20/*.bin.vcf";

    # pull haplotypes out of vcf
    system
      "vcftools --vcf $results_vcf --IMPUTE --out ${resultsDir}/$resultsName";

    compare_ok( $expected_file, $results_file, "hapfuse phases" );
}

