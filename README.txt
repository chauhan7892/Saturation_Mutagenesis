*****************************


Date: 2023-01-12 12:04:38.431818
Saturation_Mutagenesis
Personnel: Pankaj Chauhan
*****************************


 
 
Question: 
Objective: 
Outline: 



## set the git master to github
git remote add -u sm https://github.com/chauhan7892/Saturation_Mutagenesis.git

## Read Fastq and write fasta file with avg Phred-score >= 20
python3 Program/Script/Fastq_read_with_Phred.py -c 20 -i Data/Raw_Data/NNS_library_1_test.fastq -o Data/Processed_Data/NNS_library_1_processed_test.fa Data/Processed_Data/NNS_library_1_missing_IDs_test.txt

## Run the alignment
python3 Program/Script/Run_BW_parallel.py -c 2 2 -i Data/Raw_Data/WT_sequence.fa Data/Processed_Data/NNS_library_1_processed_test.fa Data/Processed_Data/NNS_library_1_pickle -o Data/Result/NNS_library_1_test.txt

## Read Fastq and write fasta file to sql with avg Phred-score >= 20
## test
python3 Program/Script/Fastq_read_with_Phred_sql.py -c 20 -i test.sqlite fasta_table NNS_library_1 Data/Raw_Data/NNS_library_1_test.fastq -o Data/Processed_Data/NNS_library_1_missing_IDs_sql_test.txt

python3 Program/Script/Run_BW_parallel_sql.py -c 2 2 -i test.sqlite fasta_table NNS_library_1 Data/Raw_Data/WT_sequence.fa -o Data/Result/NNS_library_1_sql_test.txt


## nextflow
# nextflow run SaturMut.nf -resume -with-report -with-trace -with-timeline -with-dag dag.png
