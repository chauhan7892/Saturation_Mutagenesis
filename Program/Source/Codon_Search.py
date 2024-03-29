from datetime import datetime
import os
import pandas as pd
# import pdb
import argparse
from Burrow_Wheeler_sql import *


codon_dict = {
'NAA': 'X', 'NAC': 'X', 'NAG': 'X', 'NAT': 'X', 
'NCA': 'X', 'NCC': 'X', 'NCG': 'X', 'NCT': 'X', 
'NGA': 'X', 'NGC': 'X', 'NGG': 'X', 'NGT': 'X', 
'NTA': 'X', 'NTC': 'X', 'NTG': 'X', 'NTT': 'X', 
'NNA': 'X', 'NNC': 'X', 'NNG': 'X', 'NNT': 'X',
'NNN': 'X',
'ATG':'M',
'TGG':'W',
'TTT':'F', 'TTC':'F',
'TGT':'C', 'TGC':'C',
'TAT':'Y', 'TAC':'Y',
'CAA':'Q', 'CAG':'Q',
'AAT':'N', 'AAC':'N',
'CAT':'H', 'CAC':'H',
'GAA':'E', 'GAG':'E',
'GAT':'D', 'GAC':'D',
'AAA':'K', 'AAG':'K',
'ATT':'I', 'ATC':'I', 'ATA':'I',
'TAA':'*', 'TAG':'*', 'TGA':'*',
'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'TTA':'L', 'TTG':'L',
'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R'}


primers_list = ['CTACGGCTTCCAGCCCACAAAC','TGAAGCCTTTTGAGAGGGACATC','GACTTCACCGGCTGCGTGATC',\
    'CGTGTACGCCGACAGCTTTGTG','GCGTGGCCGACTATTCTGTGC','GCGGAGGGTCGGCTAGCCATATG']

def getSeqsFromFasta(filepath):
	# fasta_seq_dict = OrderedDict()
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


def processAlignedFile(file):
    seq_dict = {}
    alignment_reads = []
    ref_seq_flag = 0
    with open(file, 'r') as f_seq:
        for line in f_seq:
            if line[0] == '\n':continue
            if line[0] == '#':
                if alignment_reads:
                    seq_dict[ref_seq_coord] = alignment_reads
                    alignment_reads = []
                ref_seq_coord = line.strip('\n').split()[1]
                ref_seq_flag = 1
            elif ref_seq_flag == 2:
                ref_seq = line.strip('\n')
            else:
                if line[0] == '>':
                    seq_id = line.strip('\n')
                else:
                    seq = line.strip('\n')
                    alignment_reads.append((seq_id,seq))

            ref_seq_flag += 1
        seq_dict[ref_seq_coord] = alignment_reads

    return seq_dict


def translateCodon(seq,seq_pos):
    codon_triplet = seq[seq_pos:seq_pos+3]
    codon_amino = codon_dict[codon_triplet]
    return codon_amino, codon_triplet

def main( ):

    start_time = datetime.now()
    input_file = args_.input_argument
    output_file = args_.output_argument
    wt_dict = getSeqsFromFasta(input_file[0])
    ref_seq = wt_dict[list(wt_dict.keys())[0]]
    aligned_seqs_dict = processAlignedFile(input_file[1])

    with open(output_file[0],'w') as f_out:
        f_out.write(f'COORD\tWT_codon\tWT_AA\tMut_codon\tMut_AA\tSeq_Index\tPrimer\n')
        for coord in aligned_seqs_dict:
            start_ref = int(coord.split(":")[0])
            end_ref = int(coord.split(":")[1])
            text = ref_seq[start_ref:end_ref]
            # print(coord)
            primer_indices = approxPatternMatchFreq(text, primers_list, 1, 0)
            start_primer = -1
            end_primer = -1
            for i in range(len(primers_list)):
                x = primer_indices[i]
                if len(x)>0:
                    if int(x[0]) < 2:
                        start_primer, start_coord, start_len = i, x[0], len(primers_list[i])
                    elif int(x[0]) > 10:
                        end_primer, end_coord, end_len = i, x[0], len(primers_list[i])

            # print(start_primer, start_coord, start_len)
            # print(end_primer, end_coord, end_len)
            # print(text)
            for aligned_seq_tuple in aligned_seqs_dict[coord]:
                aligned_seq = aligned_seq_tuple[1]
                # print(aligned_seq)
                if start_primer != -1 and end_primer == -1:
                    query_seq_used = aligned_seq[start_coord+start_len:]
                    ref_seq_used = text[start_coord+start_len:]
                    coord_used = str(start_coord+start_len) + ":" + str(len(ref_seq_used))
                    if start_primer == 4: 
                        query_seq_used = query_seq_used[2:]
                        ref_seq_used = ref_seq_used[2:]
                        coord_used = str(start_coord+start_len+2) + ":" + str(len(ref_seq_used))
                elif start_primer  != -1 and end_primer != -1:
                    query_seq_used = aligned_seq[start_coord+start_len:end_coord]
                    ref_seq_used = text[start_coord+start_len:end_coord]
                    coord_used = str(start_coord+start_len) + ":" + str(end_coord)
                    if start_primer == 4: 
                        query_seq_used = query_seq_used[2:]
                        ref_seq_used = ref_seq_used[2:]
                        coord_used = str(start_coord+start_len+2) + ":" + str(end_coord)
                # print(query_seq_used)
                # print(ref_seq_used)

                ### read the codon
                codons_len = len(ref_seq_used)//3
                for i in range(codons_len):
                    current_index = i*3
                    try:
                        query_codon_aa,query_codon_triplet = translateCodon(query_seq_used,current_index)
                    except KeyError:
                        query_codon_aa = 'X'
                        query_codon_triplet = 'XXX'
                    try:
                        ref_codon_aa,ref_codon_triplet = translateCodon(ref_seq_used,current_index)
                    except KeyError:
                        ref_codon_aa = 'X'
                        ref_codon_triplet = 'XXX'
                    f_out.write(f"{coord_used}\t{ref_codon_triplet}\t{ref_codon_aa}\t{query_codon_triplet}\t{query_codon_aa}\t{int(current_index/3)}\t{start_primer}\n")


    end_time = datetime.now()
    print(f"PROGRAM IS COMPLETED||Duration||H:M:S||{end_time - start_time}|")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script to read codons")
    parser.add_argument('-i', nargs='+', dest="input_argument",default="check inputs", help="input files")
    parser.add_argument('-o', nargs='+', dest="output_argument",default="check outputs", help="write to a file")
    args_ = parser.parse_args()
    main(  )
    print('***** well done *****')
