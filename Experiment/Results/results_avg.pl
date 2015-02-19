#! /usr/bin/perl
# Author: Dominic Bett
# Bioinformatics Project: Experiment on Diffusion Kernel

use strict;
use warnings;
use extract_kernel_values;

my $dk_NetworkRef = readFile2Array("../Diffusion/dk_avg_Network_950.txt");
my $dk_avgOriginalRef = readFile2Array("../../Propable PPI/Raw/avg_Network_950.txt");
my $correspondingPPIRef = getCorrespondingPPIValues($dk_NetworkRef);


my %dk_expt = ();
foreach (@$correspondingPPIRef) {
	my @temp = split("\t", $_);
	$temp[0] =~ s/^\s+|\s+$//g;
	$temp[1] =~ s/^\s+|\s+$//g;
	$temp[2] =~ s/^\s+|\s+$//g;
	my $key = $temp[0]."\t".$temp[1];
	#$key =~ s/^\s+|\s+$//g;
	my $value = $temp[2];
	#$value =~ s/^\s+|\s+$//g;
	$dk_expt{$key} = $value;
}

my %dk_original = ();
foreach (@$dk_avgOriginalRef) {
	my @temp = split("\t", $_);
	my $key = $temp[0]."\t".$temp[1];
	#$key =~ s/^\s+|\s+$//g;
	my $value = $temp[2];
	#$value =~ s/^\s+|\s+$//g;
	$dk_original{$key} = $value;
}

foreach (sort keys %dk_expt) {
#	print $_."\t".$dk_expt{$_},"\n";
}

#print $dk_expt{"ARG1\tCHIA"},"\n";
#print "+++++++++++++++++\n";
foreach (sort keys %dk_original) {

#	print $_,"\n";
}

my @dk_new = ();
foreach (sort keys %dk_original) {
	#print $_,"\n";
	#print $dk_expt{$_};
	push (@dk_new, $_."\t".$dk_expt{$_});
}

writeArraySort(\@dk_new, "res_Network_950.txt");

#writeArraySort($correspondingPPIRef, "../Kernel Values/Network_950.txt");

exit;