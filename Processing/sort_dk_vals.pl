#! usr/bin/perl
#author: Dominic Bett

use strict;
use warnings;
use laplacian;
use extract_kernel_values;

my $array1 = readFile2Array("../Raw_Data/Allergy_and_Asthma.txt");
my $array2 = readFile2Array("../Actual DK Values/ppi_kernel_Asthma_and_Allergy.txt");
my @sorted = ();
my %hash2 = ();

print "Sorting PPIs....\n";

foreach (@$array2) {
	my @temp = split("\t", $_);
	$hash2{$temp[0]."\t".$temp[1]} = $temp[2];
}

foreach (@$array1) {
	my @temp = split("\t", $_);
	my $key = $temp[0]."\t".$temp[1];
	push(@sorted, join("\t", $key, $hash2{$key}));
}

print "Sorted!!\n";

my $aFile = "../Sorted_Kernel Output/Actual/sorted_actual_dk_Allergy_and_Asthma.txt";

unless(open(WF, '>'.$aFile)){
	print "Could not write to file: $aFile.\n";
	exit;
}

print "Writing to file: $aFile...\n";
foreach my $element (@sorted){	
	print WF "$element\n";
}

close WF;

print "Done!\n";

exit;