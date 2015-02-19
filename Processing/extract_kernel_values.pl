#! /usr/bin/perl
# Author: Dominic Bett
# Bioinformatics Project: Extract corresponding difussion kernel values for PPIs

use strict;
use warnings;
use laplacian;
use extract_kernel_values;

#######################
=Already done... no need to run again
my $diffKernelFile = "../Diffusion/dk_Allergy_Asthma.matrix";
my $rawdataRef = readFile2Array($diffKernelFile);

my $diffKernelRef = cleanData($rawdataRef);

writeCleanPPIsToFile("CleanDK/clean_dk_Allergy_Asthma.txt", $diffKernelRef);

my $correspondingPPIRef = getCorrespondingPPIValues($diffKernelRef);

#my $fName = "PPI_val_extract.txt";
#writeArraySort($correspondingPPIRef, $fName);


my $asthmaPPIFile = "../Data/Allergy_and_Asthma.dat";
my $asthmaPPIRef = readFileToArray($asthmaPPIFile);

my ($nonAsthmaPPI_KernelRef, $asthmaPPI_KernelRef) = nonAsthmaPPIs($correspondingPPIRef, $asthmaPPIRef);
writeArraySort($nonAsthmaPPI_KernelRef , "../DK Non-Existing PPI/nonPPI_kernel_Asthma_and_Alergy.txt");
#writeArraySort($asthmaPPI_KernelRef , "../Actual DK Values/ppi_kernel_Network_950.txt");

#my $nonPPI_with_higherVals = nonPPIWithHigherVals($nonAsthmaPPI_KernelRef, $asthmaPPI_KernelRef);
#writeArraySort($nonPPI_with_higherVals, "non_PPI_with_Higher_Vals.txt");

=cut
###################
=comm method 3

my $propablePPIsRef = evalPropablePPIsMethod3($diffKernelRef, $laplMatrixRef);
writeArraySort($propablePPIsRef, "m3PropablePPIs.txt");
==cut
=comm alg 3
my ($propablePPIsRef, $diagonalsRef) = evalUsingDiag($diffKernelRef, $laplMatrixRef);
writeArraySort($propablePPIsRef, "diags_rawPropablePPIs.txt");
writeArraySort($diagonalsRef, "lapl_kern_diagonals.txt");

my $uniquePropablePPIs = removeDuplicates($propablePPIsRef);
writeArraySort($uniquePropablePPIs, "diags_UniquePropablePPIs.txt");
=cut
##############
#=comment alg2
my $diffKernelRef = readFile2Array("CleanDK/clean_dk_Network_950.txt");
my $laplMatrixRef = readFile2Array("../Laplacian/Network_950.matrix");
my $actualDKRef = readFile2Array("../Actual DK Values/ppi_kernel_Network_950.txt");

my $propablePPIsRef = evalPropablePPIs($diffKernelRef, $laplMatrixRef, $actualDKRef);
#writeArraySort($propablePPIsRef, "rawPropablePPIs.txt");
my $uniquePropablePPIs = removeDuplicates($propablePPIsRef);
writeArraySort($uniquePropablePPIs, "../Propable PPI/Raw/minim_Network_950.txt");


$propablePPIsRef = evalPropablePPIs_Average($diffKernelRef, $laplMatrixRef, $actualDKRef);
$uniquePropablePPIs = removeDuplicates($propablePPIsRef); #comment use Laplacian
writeArraySort($uniquePropablePPIs, "../Propable PPI/Raw/avg_Network_950.txt");
#=cut

######################
=comment out
my $diffKernelFile= "diff_kernel.matrix";
my $diffKernelRef = readFileToArray($diffKernelFile);

my $correspondingPPIRef = getCorrespondingPPIValues($diffKernelRef);

my $fName = "PPI_val_extract.txt";
writeHashSort($correspondingPPIRef, $fName);

my $cleanPPIFile = "clean_PPIs.dat";
my $actualPPIRef = readFileToArray($cleanPPIFile);
my ($actualPPIExtractRef, $possiblePPIExtractRef) = filterActual_PossiblePPI($correspondingPPIRef, $actualPPIRef);

my $actualPPIExtractFile = "actual_PPI_Extract.ppi";
my $possiblePPIExtractFile = "possible_PPI_Extract.ppi";
writeArraySort($actualPPIExtractRef, $actualPPIExtractFile);
writeArraySort($possiblePPIExtractRef, $possiblePPIExtractFile);

my $truePossiblePPIsRef = findTruePossiblePPIs($actualPPIExtractRef, $possiblePPIExtractRef);
my $truePossiblePPIsFile = "true_Possible_PPIs.ppi";
writeArraySort ($truePossiblePPIsRef, $truePossiblePPIsFile);

print "Done!\n";
=cut

=comm missing vs actual
my $uniquePsRef = readFileToArray("../Raw_Data/Allergy_and_Asthma.txt");
my $resultNetworkRefs = readFile2Array("../Propable PPI/Raw/minim_Network_950.txt");
my ($actualPPIRef, $missingPPIRef) = separate_actual_missing($uniquePsRef, $resultNetworkRefs);
writeArraySort($actualPPIRef, "../Propable PPI/Actual/minim_Network_950.txt");
writeArraySort($missingPPIRef, "../Propable PPI/Missing/minim_Network_950.txt");

$resultNetworkRefs = readFile2Array("../Propable PPI/Raw/avg_Network_950.txt");
($actualPPIRef, $missingPPIRef) = separate_actual_missing($uniquePsRef, $resultNetworkRefs);
writeArraySort($actualPPIRef, "../Propable PPI/Actual/avg_Network_950.txt");
writeArraySort($missingPPIRef, "../Propable PPI/Missing/avg_Network_950.txt");
=cut
exit;