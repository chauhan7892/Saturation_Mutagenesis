from datetime import datetime
import os
import pandas as pd
# import pdb
import argparse

## define a function to process each chunk
def process_chunk(chunk):
    # group by Coord, Primer, and Seq_Index, and count the number of occurrences of Mut_AA
    chunk['AA_Freq'] = chunk.groupby(['Coord', 'Primer', 'Seq_Index'])['Mut_AA'].transform('count')
    # drop duplicates  
    chunk.drop_duplicates(inplace=True)
    # return the processed chunk
    return chunk


def main( ):
    start_time = datetime.now()
    input_file = args_.input_argument
    output_file = args_.output_argument

    # define the chunksize
    chunksize = 100000

    # define the columns to select
    columns = ['Coord', 'Primer', 'Seq_Index', 'Mut_AA', 'WT_AA']

    # create an iterator over the input file
    reader = pd.read_csv(input_file[0], sep='\t', usecols=columns, chunksize=chunksize)
    # # process each chunk and concatenate the results
    result_mut = pd.concat([process_chunk(chunk) for chunk in reader])

    # group and aggregate the Count column by summing
    result_mut_agg = result_mut.groupby(['Coord', 'Primer', 'Seq_Index', 'WT_AA', 'Mut_AA']).agg({'AA_Freq': 'sum'}).reset_index()
 
    # drop duplicates
    result_mut_agg.drop_duplicates(inplace=True)
    # write the result to the output file
    result_mut_agg.to_csv(output_file[0], sep='\t', index=False)

    end_time = datetime.now()
    print(f"PROGRAM IS COMPLETED||Duration||H:M:S||{end_time - start_time}|")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script to read codons")
    parser.add_argument('-i', nargs='+', dest="input_argument",default="check inputs", help="input files")
    parser.add_argument('-o', nargs='+', dest="output_argument",default="check outputs", help="write to a file")
    args_ = parser.parse_args()
    main(  )
    print('***** well done *****')
