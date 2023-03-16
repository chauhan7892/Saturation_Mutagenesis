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

python3 Program/Script/Run_BW_parallel_sql.py -c 2 2 -i test.sqlite fasta_table NNS_library_1 Data/Raw_Data/WT_sequence.fa -o Data/Result/NNS_library_1_sql_aligned_ids_test.txt

python3 Program/Script/Codon_Search_sql.py -i test.sqlite fasta_table NNS_library_1 Data/Raw_Data/WT_sequence.fa Data/Result/NNS_library_1_sql_aligned_ids_test.txt -o Data/Result/NNS_library_1_sql_codon.txt

python3 Program/Script/Codon_Freq.py -i Data/Result/NNS_library_1_sql_codon.txt -o Data/Result/NNS_library_1_sql_codon_freq.txt
## nextflow
nextflow run sm_sql.nf -c sm.config -resume -with-report -with-trace -with-timeline -with-dag dag.png

grep '0:151\|111:262\|211:362\|311:462\|411:562\|511:662'


python3 Program/Script/Fastq_read_with_Phred_sql.py -c 30 -i lib1.sqlite fasta_table NNS_library_1 Data/Raw_Data/NNS_library_1.fastq -o Data/Processed_Data/NNS_library_1_missing_IDs_sql.txt

python3 Program/Script/Run_BW_parallel_sql.py -c 2 2 -i lib1.sqlite fasta_table NNS_library_1 Data/Raw_Data/WT_sequence.fa -o Data/Result/NNS_library_1_sql_aligned_ids.txt

python3 Program/Script/Run_BW_parallel_sql_v2.py -c 2 2 -i lib1.sqlite fasta_table NNS_library_1 Data/Raw_Data/WT_sequence.fa -o Data/Result/NNS_library_1_sql_aligned_ids.txt

python3 Program/Script/Codon_Search_sql.py -i lib1.sqlite fasta_table NNS_library_1 Data/Raw_Data/WT_sequence.fa Data/Result/NNS_library_1_sql_aligned_ids.txt -o Data/Result/NNS_library_1_sql_aligned_codon.txt


python3 Program/Script/Fastq_read_with_Phred_sql_v2.py -c 30 -i test.sqlite NNS_library_1 Data/Raw_Data/NNS_library_1_test.fastq -o Data/Processed_Data/NNS_library_1_missing_IDs_sql_test.txt

python3 Program/Script/Run_BW_parallel_sql_v3.py -c 2 2 -i test.sqlite NNS_library_1 1000 Data/Raw_Data/WT_sequence.fa -o Data/Result/NNS_library_1_sql_aligned_ids_test.txt

python3 Program/Script/Codon_Search_sql_v2.py -i test.sqlite NNS_library_1 Data/Raw_Data/WT_sequence.fa Data/Result/NNS_library_1_sql_aligned_ids_test.txt -o Data/Result/NNS_library_1_sql_codon_test.txt

### longer set
python3 Program/Script/Fastq_read_with_Phred_sql_v2.py -c 30 -i lib1.sqlite NNS_library_1 Data/Raw_Data/NNS_library_1.fastq -o Data/Processed_Data/NNS_library_1_missing_IDs_sql.txt

python3 Program/Script/Run_BW_parallel_sql_v3.py -c 2 2 -i lib1.sqlite NNS_library_1 1000 Data/Raw_Data/WT_sequence.fa -o Data/Result/NNS_library_1_sql_aligned_ids_final.txt


