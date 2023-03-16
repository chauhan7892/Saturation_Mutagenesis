import sys, time, os
import argparse
from datetime import datetime

def main( ):
    start_time = datetime.now()
    input_file = args_.input_argument
    output_file = args_.output_argument

    with open(output_file[0], 'w') as f_out:
        with open(input_file[0], 'r') as f_in_codon: #open file
            coord_flag = '-1'
            for line in f_in_codon:
                line_attrs = line.strip('\n').split('\t')
                coord = line_attrs[0]
                pos = int(line_attrs[2])
                mut_AA = line_attrs[3]
                freq_count = int(line_attrs[4])
                if coord != coord_flag:




    output_txt = 'Mutation' + '\t' + 'Position'+ '\t' + 'Delta_Delta_G' + '\t' + 'Feature' + '\n'

    end_time = datetime.now()
    print(f"PROGRAM IS COMPLETED||Duration||H:M:S||{end_time - start_time}|")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="A script to make make table for foldx energy")
	parser.add_argument('-i', nargs='+', dest="input_argument",default="something.txt", help="something file")
	parser.add_argument('-o', nargs='+', dest="output_argument",default="something_.txt", help="something file")
	args_ = parser.parse_args()
	main(  )
	print('done')