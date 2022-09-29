# Preprocess_aptamers_for_AptaSUITE
Remove aptamers with primer sequences within the random region and output **fastq files** for AptaSUITE

## *Workflow*

1. Merge FASTQ Files (R1 and R2 files; Single R1 file can also be applied)
2. Separate forward and reverse sequences from the merge fastq file
3. Reverse Complement reverse sequences
4. Merge all sequences (All sequences are forward now)
5. Demultiplex different barcodes and output the corresponding fasta files (such as (+)_barcode_name.fasta in the folder of raw)
6. Output the corresponding fasta files without primers (such as (+)_removeprimers_barcode_name.fasta in the folder of raw)
7. Remove aptamers with primer sequences within the random region and output **fastq files** for each barcode in the Clean_data (Input for AptaSUITE)

## 1. *Prepare your primer and bacodes files*

#### *Refer to the supplemetary part of "Preprocess_aptamers_for_AptaSUITE.R"*

* Create your ouput folder
* Primer and bacodes files should be in the path of ouput folder

## 2. *Install required packages*

<bar> Make sure you have installed the following packages:
* Biostrings
* fuzzywuzzyR
* QuasR
* tidyr
* dplyr
* optparse

## 3. *Run the following script*

``` Rscript Preprocess_aptamers_for_AptaSUITE.R --output {the path of your output folder} --data {the path of the data folder} ```
