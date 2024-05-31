# Evolutionary-Potential-Hypothesis
Intra- and interspecific response to testosterone

morphological_data.csv
ID - column containing the individual ID for each individual
Species - column containing the species to which each individual belongs
	SCME - Sceloporus merriami
	SCUN - Sceloporus undulatus
	SCVI - Sceloporus virgatus
Sex - column containing the sex of each individual
	F - Female
	M - Male
Treatment - column containing information about the implant recieved by an individual
	CONT - Control implant
	TEST - Testosterone implant
SVL - column containing the snout-to-vent length (SVL) of each individual, measured in millimeters
Mass - column containing the mass of each individual, measured in grams
Test_conc - column containing the circulating level of testosterone at the time of tissue collection for each individual, measured in nanograms per milliliter

ReadCounts.csv
Geneid - column containing the Gene ID
Chr - column containing information about which chromosome a gene resides on
Start - column containing information about the position on the chromosome a gene begins
End - column containing information about the position on the chromosome a gene ends
Strand - column containing information about which strand each exon resides on
Length - column containing information about the number of basepairs in a gene
Every remaining column is an individual ID, with each datum representing the number of transcripts counted for each gene (row)
