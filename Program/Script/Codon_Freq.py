from datetime import datetime
import os
import pandas as pd
# import pdb
import argparse
from tqdm import tqdm

## define a function to process each chunk
def process_chunk(chunk):
    # group by Coord, Primer, and Seq_Index, and count the number of occurrences of Mut_AA
    aa_wt_freq = chunk.groupby(['Coord', 'Seq_Index', 'Mut_AA'])[['WT_codon', 'Mut_codon']].apply(lambda x: (x['WT_codon'] == x['Mut_codon']).sum()).reset_index()
    aa_wt_freq = aa_wt_freq.rename(columns={0: 'WT_AA_Freq'})
    chunk = pd.merge(chunk, aa_wt_freq, on=['Coord', 'Seq_Index', 'Mut_AA'])
    chunk['AA_Freq'] = chunk.groupby(['Coord', 'Seq_Index', 'Mut_AA'])['Mut_AA'].transform('count')
    # drop duplicates
    chunk.drop_duplicates(['Coord', 'Seq_Index', 'Mut_AA'], inplace=True)
    chunk.reset_index(drop=True, inplace=True)
    # return the processed chunk
    return chunk

def main( ):
    start_time = datetime.now()

    parameters = args_.config_argument
    input_file = args_.input_argument
    output_file = args_.output_argument

    # define the chunksize
    BATCH_SIZE = int(parameters[0])# 100000

    # define the columns to select
    columns = ['Coord', 'Seq_Index', 'WT_codon', 'WT_AA', 'Mut_codon', 'Mut_AA']

    # create an iterator over the input file
    reader = pd.read_csv(input_file[0], sep='\t', usecols=columns, chunksize=BATCH_SIZE)
    # # process each chunk and concatenate the results
    chunks_result = [process_chunk(chunk) for chunk in tqdm(reader)]
    result_mut = pd.concat(chunks_result, ignore_index=True)

    # group and aggregate the Count column by summing
    result_mut_agg = result_mut.groupby(['Coord', 'Seq_Index', 'WT_AA', 'Mut_AA']).agg({'AA_Freq': 'sum', 'WT_AA_Freq': 'sum'}).reset_index()
 
    # drop duplicates
    result_mut_agg.drop_duplicates(['Coord', 'Seq_Index', 'Mut_AA'], inplace=True)
    # result_mut_agg.drop(['WT_codon',  'Mut_codon'], axis=1)

    # write the result to the output file
    result_mut_agg.to_csv(output_file[0], sep='\t', index=False)

    end_time = datetime.now()
    print(f"PROGRAM IS COMPLETED||Duration||H:M:S||{end_time - start_time}|")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script to read codons")
    parser.add_argument('-c', nargs='+', dest="config_argument",default="something", help="input parameters")
    parser.add_argument('-i', nargs='+', dest="input_argument",default="check inputs", help="input files")
    parser.add_argument('-o', nargs='+', dest="output_argument",default="check outputs", help="write to a file")
    args_ = parser.parse_args()
    main(  )
    print('***** well done *****')
