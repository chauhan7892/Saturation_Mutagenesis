import numpy as np
from result_class import SequenceMappingInfo
# import pickle
import sqlite3

# this array also acts as partial suffix array
def suffixArray(text, MULTIPLIER = 1):

    if text[-1] != '$':
        text = text + '$'

    # index = sorted(range(len(text)), key = lambda i: (text + text)[i:]) # old
    index = sorted(range(len(text)), key = lambda i: (text)[i-1:])

    # return [(index.index(x), x) for x in index if not x % MULTIPLIER] # old
    return [x for x in index if not x % MULTIPLIER]


def bwtFromSuffixArray(text, sa = 'None', MULTIPLIER = 1):

    if sa == 'None':
        sa = suffixArray(text, MULTIPLIER)
        # print(sa)
    # return ''.join([text[(sa[i]+len(text)-1)%len(text)] for i in range(len(text))]) # old

    return ''.join([text[(sa[i]+len(text)-2)%len(text)] for i in range(len(text))])


def mapLast2First(bwt):
    first_col = sorted(bwt)
    map_index = []
    for char in bwt:
        index = first_col.index(char)
        map_index.append(index)
        first_col[index] = "#"

    return map_index


def mapFirst2Last(bwt):
    map_l2f = mapLast2First(bwt)
    map_f2l = [map_l2f.index(ind) for ind in range(len(map_l2f))]

    return map_f2l


def inverseBWT(bwt):
    first_col = sorted(bwt)
    map_f2l = mapFirst2Last(bwt)
    start = bwt.index("$")
    nextchar = first_col[start]
    wholestr = ""
    while nextchar != "$":
        wholestr += nextchar
        start = map_f2l[start]
        nextchar = first_col[start]

    return wholestr + "$"


def createCountDict(last_col, alphabet):
    count_dict = dict()

    for a in alphabet:
        count_dict[a] = [0] * (len(last_col) + 1)

    for i in range(len(last_col)):
        a_via_id = last_col[i]
        for a in alphabet:
            if a == a_via_id:
                count_dict[a][i + 1] = count_dict[a][i] + 1
            else:
                count_dict[a][i + 1] = count_dict[a][i]

    return count_dict


def createFirstOccur(count_dict, alphabet, bwt_len):
    first_occur = dict()
    curr_index = 0
    for a in sorted(alphabet):
        first_occur[a] = curr_index
        curr_index += count_dict[a][bwt_len]

    return first_occur


def bwMatching(pattern, last_col, first_occur, count_dict):
    top = 0
    bottom = len(last_col) - 1
    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            last_short = last_col[top : (bottom + 1)]
            if symbol in last_short:
                top = first_occur[symbol] + count_dict[symbol][top]
                bottom = first_occur[symbol] + count_dict[symbol][bottom + 1] - 1
            else:
                return False, False
        else:
            return top, bottom


def patternMatchFreq(text, patterns, MULTIPLIER):

    if text[-1] != '$':
        text = text + '$'

    # alphabet = ['$', 'A', 'C', 'G', 'T']
    alphabet = list(set(text))
    sa = suffixArray(text, MULTIPLIER = MULTIPLIER)
    # print(sa)
    bwt = bwtFromSuffixArray(text, sa = sa, MULTIPLIER = MULTIPLIER)
    # print(bwt)
    count_dict = createCountDict(bwt, alphabet)
    first_occur = createFirstOccur(count_dict, alphabet, len(bwt))

    num_matchs = []
    indices = []
    for pattern in patterns:
        top, bottom = bwMatching(pattern, bwt, first_occur, count_dict)
        if top:
            match_num = bottom - top + 1
            x = list((np.array(sa[top:bottom+1]) + len(bwt) -1)%len(bwt))
        else:
            match_num = 'nan'
            x = ['nan']
        num_matchs.append(match_num)
        indices.append(x)

    print (" ".join(map(str, num_matchs)))
    print (" ".join(map(str, indices)))


def approxMatch(pattern, text_pattern, MISMATCH_LEN):
    mismatch = 0
    for i in range(len(pattern)):
        if pattern[i] != text_pattern[i]:
            mismatch += 1
            if mismatch > MISMATCH_LEN:
                return False

    return True


def approxPatternMatchFreq(text, patterns, MULTIPLIER, MISMATCH_LEN):

    if text[-1] != '$':
        text = text + '$'

    # Alphabet = ['$', 'A', 'C', 'G', 'T']
    alphabet = list(set(text))
    sa = suffixArray(text, MULTIPLIER = MULTIPLIER)
    bwt = bwtFromSuffixArray(text, sa = sa, MULTIPLIER = MULTIPLIER)
    count_dict = createCountDict(bwt, alphabet)
    first_occur = createFirstOccur(count_dict, alphabet, len(bwt))
    
    pattern_indices_list = []
    for pattern in patterns:
        pattern_indices = set()
        n = len(pattern)
        k = n // (MISMATCH_LEN+1)
        seed_and_index_tuple = [(pattern[i*k:(i+1)*k], i*k) for i in range(MISMATCH_LEN)] + [(pattern[MISMATCH_LEN*k:n], MISMATCH_LEN*k)]
        for seed, seed_index in seed_and_index_tuple:
            top, bottom = bwMatching(seed, bwt, first_occur, count_dict)
            if top:
                seed_pos = list((np.array(sa[top:bottom+1]) + len(bwt) -1)%len(bwt))
                for pos in seed_pos:
                    pattern_index = pos - seed_index
                    if pattern_index < 0:
                        continue
                    elif pattern_index + len(pattern) > len(text):
                        continue
                    elif approxMatch(pattern, text[pattern_index: pattern_index + len(pattern)], MISMATCH_LEN):
                        pattern_indices.add(pattern_index)
            else:
                continue

        # return(list(pattern_indices))
        pattern_indices_list.append(list(pattern_indices))

    return pattern_indices_list

# ***********************************************************

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


def fasta_sql_data(cursor, TABLE_NAME, q_seq_id):
    cursor.execute(f"SELECT seq FROM {TABLE_NAME} WHERE (seq_id = ?)", (q_seq_id,))
    fasta_seq = cursor.fetchone()[0]
    return fasta_seq

def approxPatternMatchFreqWithClass(db_path, TABLE_NAME, text, patterns_ids, MULTIPLIER, MISMATCH_LEN):

    if text[-1] != '$':
        text = text + '$'

    # alphabet = ['$', 'A', 'C', 'G', 'T']
    # with open(pickle_file, 'rb') as f_pickle:
    #     seq_dict = pickle.load(f_pickle)

    alphabet = list(set(text))
    sa = suffixArray(text, MULTIPLIER)
    bwt = bwtFromSuffixArray(text, sa = sa, MULTIPLIER = MULTIPLIER)
    count_dict = createCountDict(bwt, alphabet)
    first_occur = createFirstOccur(count_dict, alphabet, len(bwt))
    mapped_result_list = []

    # open the mysql database connection and create a cursor
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    for seq_id in patterns_ids:
        pattern = fasta_sql_data(cursor, str(TABLE_NAME), str(seq_id))
        pattern_indices = set()
        n = len(pattern)
        k = n // (MISMATCH_LEN+1)
        seed_and_index_tuple = [(pattern[i*k:(i+1)*k], i*k) for i in range(MISMATCH_LEN)] + [(pattern[MISMATCH_LEN*k:n], MISMATCH_LEN*k)]
        for seed, seed_index in seed_and_index_tuple:
            top, bottom = bwMatching(seed, bwt, first_occur, count_dict)
            if top:
                seed_pos = list((np.array(sa[top:bottom+1]) + len(bwt) -1)%len(bwt))
                for pos in seed_pos:
                    pattern_index = pos - seed_index
                    if pattern_index < 0:
                        continue
                    elif pattern_index + len(pattern) > len(text):
                        continue
                    elif approxMatch(pattern, text[pattern_index: pattern_index + len(pattern)], MISMATCH_LEN):
                        pattern_indices.add(pattern_index)

            else:
                continue

        if pattern_indices:
            p = SequenceMappingInfo(seq_id,list(pattern_indices)[0],len(pattern)) # stored into class object
            mapped_result_list.append(p)

    # close the cursor and database connection
    cursor.close()
    conn.close()

    return mapped_result_list