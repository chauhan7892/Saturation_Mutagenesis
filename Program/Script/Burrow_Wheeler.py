import numpy as np
from result_class import MappedResult
import pickle


# def bwtFun(text):
#     index = sorted(range(len(text)), key = lambda i: (text + text)[(i + 1) :])
#     # print(index)
#     sortlist = ('').join([text[i] for i in index])
#     return sortlist


# this array also acts as partial suffix array
def suffixArray(text, multiplier = 1):

    if text[-1] != '$':
        text = text + '$'

    # index = sorted(range(len(text)), key = lambda i: (text + text)[i:]) # old
    index = sorted(range(len(text)), key = lambda i: (text)[i-1:])

    # return [(index.index(x), x) for x in index if not x % multiplier] # old
    return [x for x in index if not x % multiplier]


def bwtFromSuffixArray(text, sa = 'None', multiplier = 1):

    if sa == 'None':
        sa = suffixArray(text, multiplier)
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


def createFirstOccur(count_dict, alphabet, wht_len):
    first_occur = dict()
    curr_index = 0
    for a in sorted(alphabet):
        first_occur[a] = curr_index
        curr_index += count_dict[a][wht_len]

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


def patternMatchFreq(text, patterns):

    if text[-1] != '$':
        text = text + '$'

    # alphabet = ['$', 'A', 'C', 'G', 'T']
    alphabet = list(set(text))
    sa = suffixArray(text, multiplier = 1)
    # print(sa)
    bwt = bwtFromSuffixArray(text, sa = sa, multiplier = 1)
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


def approxMatch(pattern, text_pattern, d):
    error = 0
    for i in range(len(pattern)):
        if pattern[i] != text_pattern[i]:
            error += 1
            if error > d:
                return False

    return True


def approxPatternMatchFreq(text, patterns, d):

    if text[-1] != '$':
        text = text + '$'

    # Alphabet = ['$', 'A', 'C', 'G', 'T']
    alphabet = list(set(text))
    sa = suffixArray(text, multiplier = 1)
    bwt = bwtFromSuffixArray(text, sa = sa, multiplier = 1)
    count_dict = createCountDict(bwt, alphabet)
    first_occur = createFirstOccur(count_dict, alphabet, len(bwt))
    
    pattern_indices_list = []
    for pattern in patterns:
        pattern_indices = set()
        n = len(pattern)
        k = n // (d+1)
        seed_and_index_tuple = [(pattern[i*k:(i+1)*k], i*k) for i in range(d)] + [(pattern[d*k:n], d*k)]
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
                    elif approxMatch(pattern, text[pattern_index: pattern_index + len(pattern)], d):
                        pattern_indices.add(pattern_index)
            else:
                continue

        # return(list(pattern_indices))
        pattern_indices_list.append(list(pattern_indices))

    return pattern_indices_list

# ***********************************************************

def approxPatternMatchFreqWithClass(text, pickle_file, patterns_ids, d):

    if text[-1] != '$':
        text = text + '$'

    # alphabet = ['$', 'A', 'C', 'G', 'T']
    with open(pickle_file, 'rb') as f_pickle:
        seq_dict = pickle.load(f_pickle)

    alphabet = list(set(text))
    sa = suffixArray(text, multiplier = 1)
    bwt = bwtFromSuffixArray(text, sa = sa, multiplier = 1)
    count_dict = createCountDict(bwt, alphabet)
    first_occur = createFirstOccur(count_dict, alphabet, len(bwt))

    mapped_result_list = []
    for seq_id in patterns_ids:
        pattern = seq_dict[seq_id]
        pattern_indices = set()
        n = len(pattern)
        k = n // (d+1)
        seed_and_index_tuple = [(pattern[i*k:(i+1)*k], i*k) for i in range(d)] + [(pattern[d*k:n], d*k)]
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
                    elif approxMatch(pattern, text[pattern_index: pattern_index + len(pattern)], d):
                        pattern_indices.add(pattern_index)

            else:
                continue

        if pattern_indices:
            p = MappedResult(seq_id,list(pattern_indices)[0],len(pattern)) # stored into class object
            mapped_result_list.append(p)

    return mapped_result_list