#!/usr/bin/env perl

use warnings;
use strict;

my $usage = "\n\tperl convert_ltrdetector_to_finder.pl genome.fasta LtrDetector.Detector.bed > out.finder.scn\n\n";
my ($genome, $bed) = @ARGV;
die $usage unless defined $bed;

print "## LTR_FINDER\n";
print "## predictions are reported in the following way\n";
print "## s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chr\n";
print "## where:\n";
print "## s = starting position\n";
print "## e = ending position\n";
print "## l = length\n";
print "## ret = LTR-retrotransposon\n";
print "## lLTR = left LTR\n";
print "## rLTR = right LTR\n";
print "## sim = similarity\n";
print "## seq-nr = sequence number\n";

open my $genome_fh, '<', $genome or die "Cannot open $genome: $!\n";
my @list;
while (<$genome_fh>) {
	push @list, $1 if /^>(\S+)/;
}
close $genome_fh;

open my $bed_fh, '<', $bed or die "Cannot open $bed: $!\n";
my $seq_id = -1;
my $chr;
while (<$bed_fh>) {
	if (/^SeqID/i) {
		$seq_id++;
		die "More BED sequence sections than FASTA sequences\n" if $seq_id >= scalar @list;
		$chr = $list[$seq_id];
		next;
	}
	next if /^\tStart/i || /^\s*$/;

	my @fields = split;
	die "Expected 18 columns in LtrDetector output\n" unless scalar @fields == 18;
	my ($from, $to, $lLTR_end, $rLTR_str, $identity) = @fields[1,2,4,5,7];
	my ($lTSD_s, $lTSD_e, $rTSD_s, $rTSD_e) = @fields[8,9,10,11];
	$from = $lTSD_e + 1 if $from eq $lTSD_s;
	$to = $rTSD_s - 1 if $to eq $rTSD_e;

	my $LTR_len = $to - $from + 1;
	my $lLTR_len = $lLTR_end - $from + 1;
	my $rLTR_len = $to - $rLTR_str + 1;
	my $similarity = $identity / 100;
	print "$from $to $LTR_len $from $lLTR_end $lLTR_len $rLTR_str $to $rLTR_len $similarity $seq_id $chr\n";
}
close $bed_fh;
