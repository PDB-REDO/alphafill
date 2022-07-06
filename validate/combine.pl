#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

open(my $h1, "<validation-report-unoptimized.tsv") or die "kon 1 niet openen: $!";
open(my $h2, "<validation-report-optimized.tsv") or die "kon 2 niet openen: $!";

my (%d2);

while (my $line = <$h2>)
{
	chomp($line);
	my ($afid, $pdbid, $compid, $asymid, $seqid, @rest) = split(m/\t/, $line);

	$d2{"$afid-$pdbid-$asymid"} = \@rest;
}

while (my $line = <$h1>)
{
	chomp($line);
	my ($afid, $pdbid, $compid, $asymid, $seqid, @rest) = split(m/\t/, $line);

	my @r2 = @{$d2{"$afid-$pdbid-$asymid"}};

	print "$line", join("\t", @r2), "\n";
}
