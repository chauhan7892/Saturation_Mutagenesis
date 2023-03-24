import argparse 
import pandas as pd
import numpy as np
import os
from collections import Counter

def main( ):
    input_file = args_.input_argument
    output_file = args_.output_argument
    truncated_mids = []
    with open(input_file[0], 'r') as f_in:
        for line in f_in:
            mid = line.strip('\n')[1:]
            truncated_mids.append(mid)

    print(Counter(truncated_mids))
            
    tmp = output_file[0]
    
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="A script to make plot")
	parser.add_argument('-i', nargs='+', dest="input_argument",default="something.txt", help="something file")
	parser.add_argument('-o', nargs='+', dest="output_argument",default="something_.txt", help="something file")
	args_ = parser.parse_args()
	main(  )
	print('done')