from Burrow_Wheeler_sql_v2 import *
from datetime import datetime 
import os
from multiprocessing import Pool, Queue
import pandas as pd
# import pdb
import argparse
import sqlite3
from tqdm import tqdm
import cProfile
import pstats
import psutil
import resource

def getSeqsFasta(filepath):
    fasta_seq_dict = {}
    with open(filepath, 'r') as f:
        seqs = []
        header = ''
        for line in f:
            if line[0] == '\n':continue
            line = line.strip('\n')
            if line[0] == '>':
                if header == '':
                    header = line
                    continue
                else:
                    fasta_seq_dict[header] = ''.join(seqs)
                    seqs = []
                    header = line
            else:
                seqs.append(line)
        fasta_seq_dict[header] = ''.join(seqs)

    return fasta_seq_dict

# Set maximum memory usage to 80% of available RAM
MAX_MEMORY = int(psutil.virtual_memory().total * 0.8)

def main( ):

    start_time = datetime.now()
    parameters = args_.config_argument
    input_file = args_.input_argument
    output_file = args_.output_argument

    # Set memory limit for this process
    resource.setrlimit(resource.RLIMIT_AS, (MAX_MEMORY, MAX_MEMORY))

    ## Suffix array index, aligner mismatch and CPU parameters
    MULTIPLIER = 1
    MISMATCH_LEN = int(parameters[0])
    PROCESSORS = int(parameters[1])

    ## sqlite3 database
    db_path = input_file[0] # sql database 'fasta_database.sqlite'
    TABLE_NAME = input_file[1] ## table name "${sample name}"
    BATCH_SIZE = int(input_file[2]) ## set batch size


    ## Files
    wt_full_seq_file = input_file[3]
    result_file = output_file[0]

    ## Read files
    wt_dict = getSeqsFasta(wt_full_seq_file) # full gene
    text = wt_dict[list(wt_dict.keys())[0]]

    # a function used later for merging IDs at same corrdinates
    array_agg = lambda x: x.astype(str).str.cat(sep=',')

    ## open the mysql database connection and create a cursor
    try:
        conn = sqlite3.connect(db_path)
        conn.row_factory = lambda cursor, row: row[0]
        cursor = conn.cursor()
        
        ## set list to store dataframe for each batch
        grp_df_list = []
        count = 0
        # execute the query with batch size and offset
        offset = 0
        while True:
            cursor.execute(f"SELECT seq_id FROM {TABLE_NAME} LIMIT ? OFFSET ?", (BATCH_SIZE, offset))
            seq_ids_list = list(cursor.fetchall())
            if not seq_ids_list:
                # all rows have been fetched
                break      
            # # # process the batch of results
            # size of sequences Ids to be given to a processor
            each_processor_size = len(seq_ids_list)//(PROCESSORS) 
            ## divide the sequences Ids accordingly
            each_processor_chunk_ids = [seq_ids_list[i*each_processor_size:(i+1)*each_processor_size] for i in range(PROCESSORS)]  # get seq IDs
            # pdb.set_trace()  ### useful in tracing errors
            # Open mutliprocessing pool 
            pool = Pool(PROCESSORS)
            q = Queue()
            for chunk_ids in each_processor_chunk_ids:
                pool.apply_async(approxPatternMatchFreqWithClass, args=(db_path, TABLE_NAME, text, chunk_ids, MULTIPLIER, MISMATCH_LEN), callback=q.put)
                
            remaining_size = len(seq_ids_list)%(PROCESSORS)
            if remaining_size > 0: # if True
                remaining_chunk_ids = seq_ids_list[PROCESSORS*each_processor_size:len(seq_ids_list)] # Get the leftover sequence Ids
                pool.apply_async(approxPatternMatchFreqWithClass, args=(db_path, TABLE_NAME, text, remaining_chunk_ids, MULTIPLIER, MISMATCH_LEN), callback=q.put)
                
            pool.close() ## close the multi-processing 
            pool.join() ## join the outcomes from multi-processing
            q.put(None) ## important: add a None flag to Queue

            item_list = [q_item for sublist in iter(q.get, None) for q_item in sublist] ## get result class data

            pd_table = [[item.coord_position_first, item.seq_len, item.seq_id] for item in item_list]  ## make table (coordinate, sequence length, sequence ID)
            df = pd.DataFrame(pd_table, columns=['coord','seq_len','seq_id']) ## store in dataframe
            grp_df_local = df.groupby(['coord','seq_len']).agg({'seq_id': array_agg}) ## group sequences based on common coordinate and sequence length and concatenate
            grp_df_local.drop_duplicates(inplace=True)
            grp_df_list.append(grp_df_local)

            tqdm.write(f"Completed iteration {count}")
            # increment the offset for the next batch
            offset += BATCH_SIZE
            count += 1
            

    except sqlite3.Error as error:
        print("Error occurred:", error)

    finally:
        ## close the sql    
        cursor.close()
        conn.close()

    grp_df = pd.concat(grp_df_list)
    grp_df_full = grp_df.groupby(['coord','seq_len']).agg({'seq_id': array_agg})
    grp_df_full.drop_duplicates(inplace=True)
    # print( grp_df_full)

    with open(result_file,'w') as f_out:
        for index, row in grp_df_full.iterrows():
            f_out.write(f'#{index[0]}:{index[0]+index[1]}\n') #Note indexing is zero based and save only coordinates 
            seq_ids = ('\n').join([seq_id for seq_id in row['seq_id'].split(',')])
            f_out.write(f'{str(seq_ids)}\n')

    end_time = datetime.now()
    print(f"PROGRAM IS COMPLETED||Duration||H:M:S||{end_time - start_time}|")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script to run BW alignment in parallel")
    parser.add_argument('-c', nargs='+', dest="config_argument",default="something", help="input parameters")
    parser.add_argument('-i', nargs='+', dest="input_argument",default="check inputs", help="input files")
    parser.add_argument('-o', nargs='+', dest="output_argument",default="check outputs", help="write to a file")
    args_ = parser.parse_args()
    # cProfile.run('main()', sort='tottime')
    # p = pstats.Stats('stats')
    # p.sort_stats('cumulative').print_stats(20)
    main(  )

    print('***** well done *****')