from Burrow_Wheeler_sql import *
from datetime import datetime 
import os
from multiprocessing import Pool, Queue
import pandas as pd
# import pdb
import argparse
import sqlite3

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


### Requirement of sql table
### TABLE {table_name} (sample TEXT, seq_id TEXT, seq TEXT, PRIMARY KEY (sample, seq_id))

def main( ):

	start_time = datetime.now()
	parameters = args_.config_argument
	input_file = args_.input_argument
	output_file = args_.output_argument

	## Suffix array index, aligner mismatch and CPU parameters
	MULTIPLIER = 1
	MISMATCH_LEN = int(parameters[0])
	PROCESSORS = int(parameters[1])

	## sqlite3 database
	db_path = input_file[0] # sql database 'fasta_database.sqlite'
	TABLE_NAME = input_file[1] ## table name "fasta_table"
	SAMPLE = input_file[2]

	## Files
	wt_full_seq_file = input_file[3]
	result_file = output_file[0]

	## Read files
	wt_dict = getSeqsFasta(wt_full_seq_file) # full gene
	text = wt_dict[list(wt_dict.keys())[0]]

	## open the mysql database connection and create a cursor
	conn = sqlite3.connect(db_path)
	cursor = conn.cursor()
	## select seq_id from the table where sample matches the input sample
	# cursor.execute("SELECT seq_id FROM fasta_table WHERE sample = ?", (SAMPLE,)) ## hardcoded
	cursor.execute(f"SELECT seq_id FROM {TABLE_NAME} WHERE sample = ?", (SAMPLE,)) ## variable table
	seq_ids_tuple = cursor.fetchall()
	seq_ids_list = [x for x, in seq_ids_tuple]
	cursor.close()
	conn.close()

	# size of sequences Ids to be given to a processor
	each_processor_size = len(seq_ids_list)//(PROCESSORS) 

	## divide the sequences Ids accordingly
	each_processor_chunk_ids = [seq_ids_list[i*each_processor_size:(i+1)*each_processor_size] for i in range(PROCESSORS)]  # get seq IDs

	# pdb.set_trace()  ### useful in tracing errors

	# Open mutliprocessing pool 
	pool = Pool(PROCESSORS)
	q = Queue()

	for chunk_ids in each_processor_chunk_ids:
		pool.apply_async(approxPatternMatchFreqWithClass, args=(db_path, TABLE_NAME, SAMPLE, text, chunk_ids, MULTIPLIER, MISMATCH_LEN), callback=q.put)
		
	remaining_size = len(seq_ids_list)%(PROCESSORS)
	if remaining_size > 0: # if True
		remaining_chunk_ids = seq_ids_list[PROCESSORS*each_processor_size:len(seq_ids_list)] # Get the leftover sequence Ids
		pool.apply_async(approxPatternMatchFreqWithClass, args=(db_path, TABLE_NAME, SAMPLE, text, remaining_chunk_ids, MULTIPLIER, MISMATCH_LEN), callback=q.put)
		
	pool.close() ## close the multi-processing 
	pool.join() ## join the outcomes from multi-processing
	q.put(None) ## important: add a None flag to Queue

	item_list = [q_item for sublist in iter(q.get, None) for q_item in sublist] ## get result class data
	pd_table = [[item.coord_position_first, item.seq_len, item.seq_id] for item in item_list]  ## make table (coordinate, sequence length, sequence ID)
	df = pd.DataFrame(pd_table, columns=['coord','seq_len','seq_id']) ## store in dataframe
	array_agg = lambda x: ','.join(x.astype(str)) ## function to concatenate
	grp_df = df.groupby(['coord','seq_len']).agg({'seq_id': array_agg}) ## group sequences based on common coordinate and sequence length and concatenate

	# print(grp_df)

	# open the mysql database connection and create a cursor
	conn = sqlite3.connect(db_path)
	cursor = conn.cursor()
	with open(result_file,'w') as f_out:
		for index, row in grp_df.iterrows():
			# f_out.write(f'#WT {index[0]}-{index[0]+index[1]}\n') #Note indexing is zero based and save only coordinates 
			# seq_ids = ('\t').join([seq_id for seq_id in row['seq_id'].split(',')])
			# f_out.write(f'{seq_ids}\n\n')

			f_out.write(f'#WT {index[0]}:{index[0]+index[1]}\n{text[index[0]:index[0]+index[1]]}\n') ## save coordinates and sequences too
			for seq_id in row['seq_id'].split(','):
				cursor.execute(f"SELECT seq FROM {TABLE_NAME} WHERE (seq_id = ? AND sample = ?)", (seq_id, SAMPLE,))
				fasta_seq = cursor.fetchone()[0]
				f_out.write(f'{seq_id}\n{fasta_seq}\n')

		cursor.close() ##  close cursor
		conn.close() ## close database

	end_time = datetime.now()
	print(f"PROGRAM IS COMPLETED||Duration||H:M:S||{end_time - start_time}|")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="A script to run BW alignment in parallel")
	parser.add_argument('-c', nargs='+', dest="config_argument",default="something", help="input parameters")
	parser.add_argument('-i', nargs='+', dest="input_argument",default="check inputs", help="input files")
	parser.add_argument('-o', nargs='+', dest="output_argument",default="check outputs", help="write to a file")
	args_ = parser.parse_args()
	main(  )
	print('***** well done *****')
