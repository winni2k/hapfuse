#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Test::Files;
use File::Spec;
use File::Basename;
use File::Path qw(make_path remove_tree);

my $sd = $ENV{srcdir};


plan tests => 4;

# this tests erroneous input to hapfuse
my $testTag = 'test05';

my $resultsDir = "results/$testTag";
make_path($resultsDir);

my $resultsName = "$testTag.error";

my $results_vcf = File::Spec->catfile( $resultsDir, $resultsName . q/.vcf/ );
my $cmd =   "./hapfuse -Ov  -tGT -TGT,GP -o $results_vcf $sd/samples/test20/test20.chr20*.bin.vcf 2>&1";
print "Call: $cmd\n";
my $res = qx/$cmd/;
like($res, qr/Need to specify APP or GP as input tag if specifying GP as output tag/, "Missing for GP out");


$cmd =   "./hapfuse -Ov  -tGT -TGT,APP -o $results_vcf $sd/samples/test20/test20.chr20*.bin.vcf 2>&1";
print "Call: $cmd\n";
$res = qx/$cmd/;
like($res, qr/Need to specify APP as input tag if specifying APP as output tag/, "Missing for APP out");


$cmd =   "./hapfuse -Ov  -tGT,GP -TGT,APP -o $results_vcf $sd/samples/test20/test20.chr20*.bin.vcf 2>&1";
print "Call: $cmd\n";
$res = qx/$cmd/;
like($res, qr/Need to specify APP as input tag if specifying APP as output tag/, "Missing for APP out 2");


$cmd =   "./hapfuse -Ov  -tAPP -TGT,APP -o $results_vcf $sd/samples/test20/test20.chr20*.bin.vcf 2>&1";
print "Call: $cmd\n";
$res = qx/$cmd/;
like($res, qr(Need to specify input GT tag if fusing BCF/VCF files), "Missing input GT");

