package TestHelpers;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(vcfComp removeHeaderLines roundFloats);
%EXPORT_TAGS = ( ALL => [ map { '&'.$_ } @EXPORT_OK ] );

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
    $line =~ s/([,:])(\d+)(?=[,:\t\n])/$1$2.00/g;
    return $line;
}

1;
