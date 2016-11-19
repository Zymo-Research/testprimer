# TestPrimer

TestPrimer is a command-line utility that runs *in silico* PCR with primers on microbial dataset and create performance report.

## Installation

TestPrimer is currently installed in AWS AMI `TestPrimer_v1.0.0` with **SILVA 128 SSU full-align** datasets in FASTA format.

## QuickStart

**Create files for forward and reverse primer pool separately:**

	#515f 11895-13861
	GTGCCAGCAGTCGCGGTAA
	GTGCCAGCAGGAGCGGTAA
	GTGCCACCAGCCGCGGTAA
	GTGCCAGAAGTCTCGGTAA
	GTGCCAGAAGCGTCGGTAA
	GTGCCAGAAGCCTCGGTAA

**Perform *in silico* PCR:**

	testprimer pcr --fasta /data/SILVA_128_SSURef_Nr99_tax_silva_full_align_trunc.fasta --forward <forward primer pool> --reverse <reverse primer pool>

**Generate coverage report:**

	testprimer report --sql <PCR result DB file> --coverage

**Perform *in silico* PCR and generate report jointly:**

	testprimer pcr -a /data/SILVA_128_SSURef_Nr99_tax_silva_full_align_trunc.fasta -f <forward primer pool> -r <reverse primer pool> --coverage

**Edit disease-related pathogen list in `~/.testprimer/config`:** 

	[TaxaCoverage]
	pathogens = Acinetobacter,Actinobacillus,Actinomyces,Aeromonas,Alcaligenes,Amycolata
