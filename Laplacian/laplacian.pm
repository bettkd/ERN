#! /usr/bin/perl
# Author: Dominic Bett
# Bioinformatics Project: Laplacian Matrix Module

use strict;
use warnings;

####################################################
# Subroutine to read file to array
# Parameters: The input file name
# Do: Read the content of the file into an array
# Return: A reference to the array with PPI from the file
####################################################	
sub readFileToArray{

	my($fName, $uniquePsRef, $header) = @_;

	unless(open(RF, $fName)){

		print "Could not read the file: $fName...\n";

		exit;
	}
	
	print "Reading file: $fName to array...\n";

	my @input = <RF>;
	chomp @input;

	close RF;
	
	#if (scalar split("\t", $input[0]) eq 1) {
	#	shift @input;
	#}
	my @rawPPI = ();
	foreach (@input){
		my @temp = split("\t", $_);
		push(@rawPPI, $temp[0]."\t".$temp[1]);
	}
	
	return \@rawPPI;
}

####################################################
# Subroutine to remove duplicates from data
# Parameters: A reference to array to the raw PPIs data
# Do: Put the data into a single sorted array
#		Get a set of unique pair of PPIs
# Return: A reference to an array with the unique pairs of proteins (PPI)
####################################################	
sub removeDuplicates{

	my ($protRef) = @_;

	my %uniquePairs = ();
	print "Removing duplicate entries and sorting pairs...\n";

	foreach (@$protRef) {

		$_ = join ("\t", sort split ("\t", $_));

		$uniquePairs{$_} ++;
	}

	my @uniquePairs = sort keys %uniquePairs;

	return \@uniquePairs;
}

####################################################
# Subroutine to write the PPIs to file
# Parameters: A reference to a unique set of PPIs,
#				A reference to a unique set of proteins
# Do: Write the cleaned data to file
# Return: N/A
####################################################
sub writeCleanPPIsToFile{

	my ($cleanPPIsFile, $uniquePairsRef) = @_;

	unless(open(WF, '>'.$cleanPPIsFile)){

		print "Could not write to file: $cleanPPIsFile...\n";

		exit;
	}

	print "Writing data to file: $cleanPPIsFile...\n";

	foreach (@$uniquePairsRef){

		print WF "$_\n";
	}

}

####################################################
# Subroutine to find unique proteins from array
# Parameters: Reference to an array of PPIs
# Do: Puts all protein pairs into a single array,
#		Gets a unique set of the proteins
# Return: A reference to the unique set of proteins
####################################################
sub findUniqueProt {

	my ($protRef) = @_;

	my $prots = join ("\t", @$protRef);
	my @prots = split ("\t", $prots);

	print "Finding unique proteins...\n";

	my %uniquePs = ();

	foreach (@prots){

		$uniquePs{$_}++;
	}

	my @uniquePs = sort keys %uniquePs;
	return \@uniquePs;
}

####################################################
# Subroutine to determine the Laplacian Matrix
# Parameters: A reference to the unique set of proteins
# Do: Initialize the laplacian matrix with zeros(0)
# Return: An initialized laplacian matrix
####################################################
sub initializeMatrix {

	my ($matrixSize) = @_;

	print "Creating a null matrix and initialializing to zero...\n";
	my @matrix = ();

	for (my $i = 0; $i < $matrixSize; $i++) {

			for (my $j = 0; $j < $matrixSize; $j++) {

				$matrix[$i][$j] = 0;
			}
		}

	return \@matrix;
}

####################################################
# Subroutine to determine the Laplacian using Hash
# Parameters: A reference to the unique set of proteins,
#				A reference to to an array of PPIs
# Do: Call a subroutine to nitialize an empty matrix,
#		Calculate the laplacian matrix
# Return: A reference to the array containing laplacian matrix
####################################################
sub findLaplacian {

	my ($uniquePsRef, $protRef) = @_;

	my $matrixRef = initializeMatrix (scalar @$uniquePsRef);
	my %protIndex = ();
	print "Calculating the laplacian matrix using hash...\n";

	for (my $i = 0; $i < scalar @$uniquePsRef; $i++) {

		$protIndex{$$uniquePsRef[$i]} = $i;
	}
	
	foreach my $p (@$protRef) {

		my ($px, $py) = split ("\t", $p);

		$$matrixRef[$protIndex{$px}][$protIndex{$py}] ++;
		$$matrixRef[$protIndex{$py}][$protIndex{$px}] ++;
		$$matrixRef[$protIndex{$px}][$protIndex{$px}] --;
		$$matrixRef[$protIndex{$py}][$protIndex{$py}] --;
	}

	return \@$matrixRef;
}

####################################################
# Subroutine to determine the Laplacian Matrix Using Brute-Force
# Parameters: A reference to the unique set of proteins,
#				A reference to to an array of PPIs
# Do: Call a subroutine to nitialize an empty matrix,
#		Calculate the laplacian matrix with brute-force
# Return: A reference to the array containing laplacian matrix
####################################################
sub calculateLaplacian {

	my ($uniquePsRef, $protRef) = @_;

	my $matrixRef = initializeMatrix (scalar @$uniquePsRef);
	print "Calculating the laplacian matrix using brute force...\n";

	foreach my $p (@$protRef) {

		my ($px, $py) = split ("\t", $p);

		for (my $i = 0; $i < scalar @$uniquePsRef; $i++) {

			for (my $j = 0; $j < scalar @$uniquePsRef; $j++) {

				if (@$uniquePsRef[$j] eq $px and @$uniquePsRef[$i] eq $py) {

					$$matrixRef[$j][$i] = 1;
					$$matrixRef[$i][$j] = 1;
					$$matrixRef[$j][$j] -= 1;
					$$matrixRef[$i][$i] -= 1;
				}
			}
		}
	}

	return \@$matrixRef;
}

####################################################
# Subroutine to write the matrix to file
# Parameters: The output file name,
#				A reference to the laplacian matrix,
#				A reference to the unique set of proteins,
#				The header name of the file
# Do: Write the matrix into a file in tab delimited format
# Return: N/A
####################################################
sub writeMatrixToFile{
	
	my($fName, $laplacianRef, $uniquePsRef) = @_;

	unless(open(WF, '>'.$fName)){

		print "Could not write to file: $fName...\n";

		exit;
	}

	print "Writing matrix to file: $fName...\n";

	print WF "\t";

	foreach (@$uniquePsRef){

		print WF "$_\t";
	}

	for(my $k = 0; $k < scalar @$laplacianRef; $k++) {

		print WF "\n$$uniquePsRef[$k]";

		for(my $l = 0; $l < scalar @$laplacianRef; $l++) {

			print WF "\t$$laplacianRef[$k][$l]";
		}
	}

	close WF;
}

1;