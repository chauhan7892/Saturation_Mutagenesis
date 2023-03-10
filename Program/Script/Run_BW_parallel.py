from Burrow_Wheeler import *
from datetime import datetime 
import os
from multiprocessing import Pool, Queue
import pandas as pd
# import pdb
import argparse


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

# def deep_map(func, some_list):
#     if not isinstance(some_list, list):
#         return func(some_list)
#     return [deep_map(func, x) for x in some_list]


def main( ):

	start_time = datetime.now()
	parameters = args_.config_argument
	input_file = args_.input_argument
	output_file = args_.output_argument

	## Files
	wt_full_seq_file = input_file[0]
	query_seqs_file = input_file[1]
	pickle_file = input_file[2]
	result_file = output_file[0]


	## Suffix array index, aligner mismatch and CPU parameters
	MULTIPLIER = 1
	MISMATCH_LEN = int(parameters[0])
	PROCESSORS = int(parameters[1])


	## Read files
	wt_dict = getSeqsFasta(wt_full_seq_file) # full gene
	text = wt_dict[list(wt_dict.keys())[0]]
	text_substring = ''

	seq_dict = getSeqsFasta(query_seqs_file) # fasta reads 
	seq_ids_list = list(seq_dict.keys())

	# size of sequences Ids to be given to a processor
	each_processor_size = len(seq_ids_list)//(PROCESSORS) 

	## divide the sequences Ids accordingly
	each_processor_chunk_ids = [seq_ids_list[i*each_processor_size:(i+1)*each_processor_size] for i in range(PROCESSORS)]  # get seq IDs
	# each_processor_chunk = deep_map(lambda x: seq_dict[x], each_processor_chunk_ids) # get sequence values

	# pdb.set_trace()  ### useful in tracing errors

	## Open mutliprocessing pool 
	pool = Pool(PROCESSORS)
	q = Queue()

	## open a file, where you want to store the pickle data
	with open(pickle_file, 'wb') as f_pickle:
	## dump information to that file
		pickle.dump(seq_dict, f_pickle)

	for chunk_ids in each_processor_chunk_ids:
		x = pool.apply(approxPatternMatchFreqWithClass, args=(text, pickle_file, chunk_ids, MULTIPLIER, MISMATCH_LEN))
		q.put(x) ## put the result into Queue

	remaining_size = len(seq_ids_list)%(PROCESSORS)
	if remaining_size > 0: # if True
		remaining_chunk_ids = seq_ids_list[PROCESSORS*each_processor_size:len(seq_ids_list)] # Get the leftover sequence Ids
		# remaining_chunk = deep_map(lambda x: seq_dict[x], remaining_chunk_ids) 
		x = approxPatternMatchFreqWithClass(text, pickle_file, remaining_chunk_ids, MULTIPLIER, MISMATCH_LEN)
		q.put(x)


	pool.close() ## close the multi-processing 
	pool.join() ## join the outcomes from multi-processing
	q.put(None) ## important: add a None flag to Queue

	item_list = [q_item for sublist in iter(q.get, None) for q_item in sublist] ## get result class data
	table = [[item.coord_position_first, item.seq_len, item.seq_id] for item in item_list]  ## make table (coordinate, sequence length, sequence ID)

	# print(table)

	df = pd.DataFrame(table, columns=['coord','seq_len','seq_id']) ## store in dataframe
	array_agg = lambda x: ','.join(x.astype(str)) ## function to concatenate
	grp_df = df.groupby(['coord','seq_len']).agg({'seq_id': array_agg}) ## group sequences based on common coordinate and sequence length and concatenate

	# print(grp_df)

	with open(result_file,'w') as f_out:
		for index, row in grp_df.iterrows():
			f_out.write('#WT %s\n%s\n' %(index[0], text[index[0]:index[0]+index[1]]))
			for seq_id in row['seq_id'].split(','):
				f_out.write('%s\n%s\n' %(seq_id, seq_dict[seq_id]))

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
