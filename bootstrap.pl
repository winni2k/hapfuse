#!/usr/bin/perl
# bootstrap.pl                   wkretzsch@gmail.com
#                                19 Feb 2015
use warnings;
use strict;
use File::Path qw(make_path);
use File::Basename;
use File::Slurp;
use autodie;
use FindBin qw/$Bin/;

# run autoconf, etc.
print STDERR "autoreconf --install --symlink...\n";
system('autoreconf --install --symlink');

# compile gtest
my $gtestDir = "gtest-1.7.0";
my $cmd = "cd $gtestDir/make && make";
print STDERR $cmd ."\n";
system($cmd);

# find package version
my @configFile = read_file("$Bin/configure.ac");
my @lines = grep { m/^\s*AC_INIT/ } @configFile;
die "Error in configure.ac file: too many AC_INIT lines" if @lines > 1;
die "could not find AC_INIT" if @lines < 1;

my $line = shift @lines;
chomp $line;
$line =~ s/ //g;
my @line = split( /[\(\),\[\]]+/, $line );
my ( undef, $binName, $version, $email ) = @line;

my $buildDir = "$Bin/$binName.$version";

print STDERR "Creating Build dir at " . basename($buildDir) . "\n";
make_path($buildDir);
chdir "$buildDir";

print STDERR "Configuring...\n";
system "../configure";

print STDERR "make all";
system "make all";

print STDERR "Returning to base dir\n";
chdir "$Bin";

print STDERR "\nAll files built in $buildDir\n";

__END__

=head1 NAME

bootstrap.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION
Stub documentation for bootstrap.pl, 
created by template.el.

It looks like the author of this script was negligent
enough to leave the stub unedited.


=head1 AUTHOR

Warren Winfried Kretzschmar, E<lt>wkretzsch@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Warren Winfried Kretzschmar

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
