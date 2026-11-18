#!/usr/bin/env perl

#convert LtrDetector default output to LTR_retriever input format
#Usage: perl convert_ltrdetector.pl genome.fasta LtrDetector.Detector.bed > out.scn
#Author:Shujun OU (shujun.ou.1@gmail.com) 04/14/2019


use warnings;
use strict;

my $usage = "\n\tperl convert.pl genome.fasta LtrDetector.bed.scn > out.scn\n\n";

print "#this summary file is created from LtrDetector tabular output by convert_ltrdetector.pl (Shujun Ou, shujun.ou.1\@gmail.com, 04/14/2019)
#start end len lLTR_str lLTR_end lLTR_len rLTR_str rLTR_end rLTR_len identity seqid chr direction TSD lTSD rTSD\n";
#start end len lLTR_str lLTR_end lLTR_len rLTR_str rLTR_end rLTR_len identity seqid chr direction TSD lTSD rTSD motif superfamily family age(ya)

my $genome = $ARGV[0];
my $scn = $ARGV[1];
open List, "grep \\> $genome | " or die $usage;
open File, "<$scn" or die "$usage";

my @list;
while (<List>){
	s/>//g;
	my $id = (split)[0];
	push @list, $id;
	}
close List;

my ($from, $to, $chr, $LTR_len, $lLTR_len, $rLTR_len, $lLTR_end, $rLTR_str, $identity, $TSD, $lTSD, $rTSD, $direction, $seq_id);
$from=$to=$LTR_len=$from=$lLTR_end=$lLTR_len=$rLTR_str=$to=$rLTR_len=$identity=$chr=$direction=$TSD=$lTSD=$rTSD=' ';
$seq_id=-1;
while (<File>){
	if (/^SeqID/i){
		$seq_id++;
		$chr=$list[$seq_id];
		next;
		}
	next if /^\tStart/i;

	my @fields = split;
	($from, $to, $lLTR_end, $rLTR_str, $identity, $direction) = @fields[1,2,4,5,7,14];

	# correct start end coordinate of candidates
	my ($lTSD_s, $lTSD_e, $rTSD_s, $rTSD_e) = @fields[8,9,10,11];
	$lTSD = $lTSD_s eq '---' ? $lTSD_s : "$lTSD_s-$lTSD_e";
	$rTSD = $rTSD_s eq '---' ? $rTSD_s : "$rTSD_s-$rTSD_e";
	$from = $lTSD_e + 1 if $from eq $lTSD_s;
	$to = $rTSD_s - 1 if $to eq $rTSD_e;

	$LTR_len = $to - $from + 1;
	$lLTR_len = $lLTR_end - $from + 1;
	$rLTR_len = $to - $rLTR_str + 1;

	$direction = "NA" if $direction eq 0;
	print "$from $to $LTR_len $from $lLTR_end $lLTR_len $rLTR_str $to $rLTR_len $identity $seq_id $chr $direction NA $lTSD $rTSD \n";
	$from=$to=$LTR_len=$lLTR_len=$rLTR_len=$lLTR_end=$rLTR_str=$identity=$direction=$TSD=$lTSD=$rTSD="NA";
	}
close File;
