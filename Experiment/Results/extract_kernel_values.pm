#! /usr/bin/perl
# Author: Dominic Bett
# Bioinformatics Project: Module to extract corresponding difussion kernel values for PPIs

use strict;
use warnings;

####################################################
# Subroutine to read file to array
# Parameters: The input file name
# Do: Read the content of the file into an array
# Return: A reference to the array with PPI from the file
####################################################	
sub readFile2Array{

	my($fName) = @_;

	unless(open(RF, $fName)){

		print "Could not read the file: $fName...\n";

		exit;
	}
	
	print "Reading file: $fName to array...\n";

	my @input = <RF>;
	chomp @input;

	close RF;
	return \@input;
}

##################################################
sub cleanData{
	my ($rawDataRef) = @_;
	
	my @cleanedData = ();

	print "Cleaning data...\n";
	foreach (@$rawDataRef) {
		$_ =~ s/\h+/\t/g;
		push (@cleanedData, $_);
	}

	return \@cleanedData;
}

##############################################################
sub nonAsthmaPPIs{
	my ($correspondingPPIRef, $asthmaPPIRef) = @_;

	print "Finding non-asthma PPIs...\n";
	my @asthmaPPI = ();
	foreach (@$asthmaPPIRef) {
		push (@asthmaPPI, join ("\t", sort split("\t", $_)));
	}

	my %corrPPI = ();
	foreach (@$correspondingPPIRef) {
		my @temp = split("\t", $_);
		my $value = pop (@temp);
		my $key = join ("\t", @temp);
		$corrPPI{$key} = $value;
	}

	my @corrPPKeys = keys %corrPPI;

	#print @corrPPKeys;

	my @union = my @isect = my @diff = ();
	my %union = my %isect = ();
	my %count = ();

	foreach my $e (@asthmaPPI, @corrPPKeys) { $count{$e}++ }

	foreach my $e (keys %count) {
	    push(@union, $e);
	    if ($count{$e} == 2) {
	        push @isect, $e;
	    } else {
	        push @diff, $e;
	    }
	}

	my @nonAsthmaPPIs = ();
	foreach (@diff) {
		push (@nonAsthmaPPIs, ($_."\t".$corrPPI{$_}));
	}

	my @asthmaKernelPPIs = ();
	foreach (@asthmaPPI) {
		push (@asthmaKernelPPIs, ($_."\t".$corrPPI{$_}));
	}

	return \(@nonAsthmaPPIs, @asthmaKernelPPIs);
}

####################################################
sub nonPPIWithHigherVals {
	my ($nonAsthmaPPI_KernelRef, $asthmaPPI_KernelRef) = @_;

	my %nonAsthmaPPI = ();
	foreach (@$nonAsthmaPPI_KernelRef) {
		my @temp = split("\t", $_);
		my $value = pop (@temp);
		my $key = join ("\t", @temp);
		$nonAsthmaPPI{$key} = $value;
	}

	my %asthmaPPI = ();
	foreach (@$asthmaPPI_KernelRef) {
		my @temp = split("\t", $_);
		my $value = pop (@temp);
		my $key = join ("\t", @temp);
		$asthmaPPI{$key} = $value;
	}

	my @temp = sort (values %asthmaPPI);
	my $lowestTrueVal = shift @temp;
	
	my @nonPPI_with_HigherVals = ();
	foreach (keys %nonAsthmaPPI) {
		if($nonAsthmaPPI{$_} > $lowestTrueVal) {
			push (@nonPPI_with_HigherVals, $_."\t".$nonAsthmaPPI{$_});
		}
	}
	return \@nonPPI_with_HigherVals;
}

####################################################
# Subroutine to get corresponding diffusion kernel values for PPI
# Parameters: Reference to the diffusion kernel matrix
# Do: Get all combinations of unique PPI and their kernel values
# Return: A reference to the PPI and their corresponding values
####################################################
sub getCorrespondingPPIValues {

	my ($diffKernelRef) = @_;
	
	my $matrixRef = readTableToMatrix($diffKernelRef);
	my %correspondingPPI = ();

	print "Getting unique PPI with their corresponding kernel values...\n";
	for ( my $i = 2; $i < scalar @$matrixRef; $i++) {
		for (my $j = 1; $j < $i; $j++) {

			#print "$$matrixRef[$i][$j]\n";
			#$$matrixRef[$i][0] =~ s/^\s+|\s+$//g;
			#$$matrixRef[0][$j] =~ s/^\s+|\s+$//g;
			#$$matrixRef[$i][$j] =~ s/^\s+|\s+$//g;
			$correspondingPPI{join ("\t", sort ($$matrixRef[$i][0], $$matrixRef[0][$j]))} = $$matrixRef[$i][$j];		
		}
	}

	my @corrPPI = ();
	foreach (keys %correspondingPPI){
		my @temp = split("\t", $_);
		if ($temp[0] ne $temp[1]){
			push(@corrPPI, $_."\t".$correspondingPPI{$_}); #change to \t
		}
	}
	return \@corrPPI;
}

####################################################
# Subroutine to read the diffusion kernel into a matrix with rows and columns
# Parameters: Reference to the array of kernel table
# Do: Separate elements in each row into columns
# Return: A reference to the matrix created
####################################################
sub readTableToMatrix {

	my ($diffKernelRef) = @_;

	my @matrix = ();

	#print "@$diffKernelRef\n";

	print "Reading the diffusion kernel into a matrix...\n";
	my $size = scalar @$diffKernelRef;
	for( my $i = 0; $i < $size; $i++) {

		my @columns = split ("\t", $$diffKernelRef[$i]);

		for (my $j = 0; $j < $size; $j++) {
			
			$matrix [$i][$j] = $columns[$j];
		}
	}
	return \@matrix;
}
####################################################
# Subroutine to extract the actual PPIs with their values
#	and to also extract the possible PPIs with their values
# Parameters: A reference to the PPI and their correlating kernel values,
#	a reference to the actual PPIs
# Do: Check if a PPI belongs to the actual PPI or not,
#	putting each category into array
# Return: A reference to the arrays with actual PPs and possible PPIs
#	with their corresponding values
####################################################
sub filterActual_PossiblePPI {

	my ($correspondingPPIRef, $actualPPIRef) = @_;

	my @actualPPIExtract = ();
	my @possiblePPIExtract = ();

	my $counter = 0;

	print "Obtaining the actual and the possible PPIs...\n";
	foreach (keys %$correspondingPPIRef) {

		if ($counter lt scalar @$actualPPIRef) {

			push (@actualPPIExtract, join ("\t", $$actualPPIRef[$counter],
				$$correspondingPPIRef{$$actualPPIRef[$counter]}));
			$counter++;
		} else {

			push (@possiblePPIExtract, join ("\t", $_,
				$$correspondingPPIRef{$_}));
		}
	}
	return (\@actualPPIExtract, \@possiblePPIExtract);
}

####################################################
# Subroutine to identify the true possible PPIs
# Parameters: References to the arrays with the actual
#	and the possible PPIs
# Do: Use the minimum value from the actual PPIs
#	to find those possible PPIs with higher value
# Return: A reference to the PPIs that tested true
####################################################
sub findTruePossiblePPIs {

	my ($actualPPIExtractRef, $possiblePPIExtractRef) = @_;

	my $minActualPPIVal = findMinPPIValue($actualPPIExtractRef);

	my @truePossiblePPIs = ();

	print "Finding the true possible PPIs...\n";
	foreach (@$possiblePPIExtractRef) {

		my @PPIExtract = (split "\t", $_ );
		my $ExtractVal = pop @PPIExtract;
		if ( $ExtractVal ge $minActualPPIVal) {

			push (@truePossiblePPIs, join ("\t", @PPIExtract));
		}

	}

	return \@truePossiblePPIs;
}

####################################################
# Subroutine to find the smallest value from the actual PPIs
# Parameters: A reference to the actual PPI with their values
# Do: Compare to find the smallest value from the actual PPIs
# Return: A reference to the smallest value established
####################################################
sub findMinPPIValue {

	my ($actualPPIExtractRef) = @_;

	my $minActualPPIValue = 1;

	print "Finding the minimum actual PPI value...\n";
	foreach (@$actualPPIExtractRef) {

		my @PPIExtract = (split "\t", $_ );
		if ( $PPIExtract[$#PPIExtract] lt $minActualPPIValue) {

			$minActualPPIValue = $PPIExtract[$#PPIExtract];
		}
	}
	return $minActualPPIValue;
}

####################################################
# Subroutine to write an array to file
# Parameters: A reference to the array and the file name
# Do: Create a file and write the array to it in a sorted order
# Return: N/A
####################################################
sub writeArraySort{
	my($input, $aFile) = @_;

	unless(open(WF, '>'.$aFile)){
		print "Could not write to file: $aFile.\n";
		exit;
	}

	print "Writing to file: $aFile...\n";
	foreach my $element (sort @$input){	
		print WF "$element\n";
	}

	close WF;
}
##############################################
#Subroutine to evaluate Propable PPI using the Minimum DK value of actual PPI
##############################################
sub evalPropablePPIs {

	my ($diffKernelRef, $laplMatrixRef, $actualDKRef) = @_;

	my $kernelRef = readTableToMatrix($diffKernelRef);
	my $laplacianRef = readTableToMatrix($laplMatrixRef);

	my $min_val = 1;
	foreach (@$actualDKRef) {
		my @temp = split("\t", $_);

		$min_val = $temp[2] if $temp[2] < $min_val;
	}

	my @propablePPIs = ();
	
	print "Evaluating propable PPIs...\n";
	for ( my $i = 1; $i < scalar @$kernelRef; $i++) {
		my %row_nonPPI = ();
		my %row_PPI = ();
		for (my $j = 1; $j < scalar @$kernelRef; $j++) {
			if ($$laplacianRef[$i][$j] eq 0) {
				$row_nonPPI{join ("\t", sort ($$kernelRef[$i][0], $$kernelRef[0][$j]))} = $$kernelRef[$i][$j];
			} elsif ($$laplacianRef[$i][$j] eq 1) {
				$row_PPI{join ("\t", sort ($$kernelRef[$i][0], $$kernelRef[0][$j]))} = $$kernelRef[$i][$j];
			}
		}
		
		# my $min_Actual = 1;
		# foreach (values %row_PPI) {
		# 	$min_Actual = $_ if ($_ < $min_Actual);
		# }
		# #print $min_Actual, "\n";

		foreach (keys %row_nonPPI) {
			push(@propablePPIs, $_."\t".$row_nonPPI{$_}) if($row_nonPPI{$_} > $min_val);
		}

	}
	return \@propablePPIs;
}

##############################################
#Subroutine to evaluate Propable PPI using the Average DK value of actual PPI
##############################################
sub evalPropablePPIs_Average {

	my ($diffKernelRef, $laplMatrixRef, $actualDKRef) = @_;

	my $kernelRef = readTableToMatrix($diffKernelRef);
	my $laplacianRef = readTableToMatrix($laplMatrixRef);

	my $sum_values = 0;
	foreach (@$actualDKRef) {
		my @temp = split("\t", $_);
		$sum_values += $temp[2];
	}
	my $avg_val = $sum_values/ scalar @$actualDKRef;


	my @propablePPIs = ();

	print "Evaluating propable PPIs...\n";
	for ( my $i = 1; $i < scalar @$kernelRef; $i++) {
		my %row_nonPPI = ();
		my %row_PPI = ();
		for (my $j = 1; $j < scalar @$kernelRef; $j++) {
			if ($$laplacianRef[$i][$j] eq 0) {
				$row_nonPPI{join ("\t", sort ($$kernelRef[$i][0], $$kernelRef[0][$j]))} = $$kernelRef[$i][$j];
			} elsif ($$laplacianRef[$i][$j] eq 1) {
				$row_PPI{join ("\t", sort ($$kernelRef[$i][0], $$kernelRef[0][$j]))} = $$kernelRef[$i][$j];
			}
		}
		
		
		# my $sum_of_vals = 0;
		# foreach (values %row_PPI) {
		# 	$sum_of_vals += $_;
		# }
		# my $average_of_actual = $sum_of_vals/(scalar keys %row_PPI);
		

		foreach (keys %row_nonPPI) {
			push(@propablePPIs, $_."\t".$row_nonPPI{$_}) if($row_nonPPI{$_} > $avg_val);
		}

	}
	return \@propablePPIs;
}

##############################################
sub removeDuplicates {
	my ($propablePPIsRef) = @_;

	my %uniquePPIs = ();

	foreach (@$propablePPIsRef) {
		$uniquePPIs{$_}++;
	}
	my @unique = keys %uniquePPIs;
	return \@unique;
}

##############################################
#Eval Propable PPI in with DK Values greater than diagonal value
################################################
sub evalUsingDiag {

	my ($diffKernelRef, $laplMatrixRef) = @_;

	my $kernelRef = readTableToMatrix($diffKernelRef);
	my $laplacianRef = readTableToMatrix($laplMatrixRef);

	my @propablePPIs = ();
	my @diagonals = ();
	print "Evaluating propable PPIs using diagonal...\n";
	for ( my $i = 1; $i < scalar @$kernelRef; $i++) {

		my %prob_row_PPI = ();
		for (my $j = 1; $j < scalar @$kernelRef; $j++) {

			if ($$laplacianRef[$i][$j] eq 0){
				if ($$kernelRef[$i][$j] gt $$kernelRef[$i][$i]) {
					$prob_row_PPI{join ("\t", sort ($$kernelRef[$i][0], $$kernelRef[0][$j]))} = $$kernelRef[$i][$j];
				}
			}

		}

		foreach (keys %prob_row_PPI) {
			push(@propablePPIs, $_."\t".$prob_row_PPI{$_});
		}

		push(@diagonals, $$laplacianRef[$i][$i]."\t".$$kernelRef[$i][$i])
	}

	
	return (\@propablePPIs, \@diagonals);
}

##############################################
#Eval Propable PPI in with DK Values greater than the minimum of the actual horizontal diagonal values
##############################################
sub evalPropablePPIsMethod3 {

	my ($diffKernelRef, $laplMatrixRef) = @_;

	my $kernelRef = readTableToMatrix($diffKernelRef);
	my $laplacianRef = readTableToMatrix($laplMatrixRef);

	my @propablePPIs = ();

	print "Evaluating propable PPIs using method 3...\n";
	for ( my $i = 2; $i < scalar @$kernelRef; $i++) {
		my %row_nonPPI = ();
		my %row_PPI = ();
		for (my $j = 1; $j < $i; $j++) {
			if ($$laplacianRef[$i][$j] eq 0) {
				$row_nonPPI{join ("\t", sort ($$kernelRef[$i][0], $$kernelRef[0][$j]))} = $$kernelRef[$i][$j];
			} elsif ($$laplacianRef[$i][$j] eq 1) {
				$row_PPI{join ("\t", sort ($$kernelRef[$i][0], $$kernelRef[0][$j]))} = $$kernelRef[$i][$j];
			}
		}

		my $min_Actual = 1;
		foreach (values %row_PPI) {
			$min_Actual = $_ if ($_ < $min_Actual);
		}
		#print $min_Actual, "\n";

		foreach (keys %row_nonPPI) {
			push(@propablePPIs, $_."\t".$row_nonPPI{$_}) if($row_nonPPI{$_} > $min_Actual);
		}

	}
	
	return \@propablePPIs;
}

#######################################
# Separate Missing PPI from Actual PPI from non-Existing PPI List
#######################################
sub separate_actual_missing {
	my ($uniquePsRef, $resultNetworkRefs) = @_;

	my @actualPPI = ();
	my @missingPPI = ();

	#my %prob_network = probableNet($resultNetworkRefs);
	my %prob_network = ();
	foreach (@$resultNetworkRefs) {
		my @temp = split("\t", $_);
		$prob_network{$temp[0]."\t".$temp[1]} = $temp[2];
	}

	my %prots = ();
	foreach (@$uniquePsRef) {
		$prots{$_}++;
	}

	print "Computing the actual and the missing PPIs...\n";

	foreach (keys %prob_network) {
		if (exists $prots{$_}) {
			push (@actualPPI, $_."\t".$prob_network{$_});
		} else {
			push (@missingPPI, $_."\t".$prob_network{$_});
		}
	}
	return (\@actualPPI, \@missingPPI);
}
=com
sub probableNet {
	my ($resultNetworkRefs) = @_;

	my %prob_network = ();
	foreach (@$resultNetworkRefs) {
		my @temp = split("\t", $_);
		$prob_network{$temp[0]."\t".$temp[1]} = $temp[2];
	}
	return %prob_network;
}
=cut
1;