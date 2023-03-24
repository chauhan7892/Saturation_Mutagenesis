import argparse 
import pandas as pd
# from matplotlib.patches import Circle, Wedge
from matplotlib import rc, rcParams
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt  
import matplotlib.ticker as ticker
from matplotlib.collections import PatchCollection

import numpy as np
from collections import OrderedDict
import os
# # Latex features
rc('text', usetex=False)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})   
rc('axes', labelsize=20)
rc('axes', titlesize=22)   # fontsize of the x and y labels
rc('xtick', labelsize=14)    # fontsize of the tick labels
rc('ytick', labelsize=8)    # fontsize of the tick labels
rc('legend', fontsize=12)    # legend fontsize
# # plt.rc('figure', titlesize=25)  # fontsize of the figure title
# rcParams['font.size'] = '16'
# rcParams['text.latex.preamble'] = r'\usepackage{sfmath}'
# rcParams['ps.fonttype'] = 42
# rcParams['svg.fonttype'] = 'none' 
# # rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

AA_col = {
	'R':'green',
	'C':'red',
	'P':'blue',
	'N':'black'
}

amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']


def main( ):
	input_file = args_.input_argument
	output_file = args_.output_argument

	# read foldx avg delta delta G file for all AA
	# define the columns to select
	columns = ['Coord', 'Seq_Index', 'Mut_AA', 'WT_AA', 'AA_Freq', 'WT_AA_Freq']

	# create an iterator over the input file
	df = pd.read_csv(input_file[0], sep='\t', usecols=columns)
	filt_coord = ["134:209", "233:311", "23:111", "332:411", "434:511", "533:662"]
	df['Coord_Start'] = df['Coord'].apply(lambda x: pd.Series(int(str(x).split(":")[0])))
	df_filt = df[df['Coord'].isin(filt_coord)].sort_values(by = ['Coord_Start', 'Seq_Index'], ascending = [True, True]).reset_index(drop=True)
	df_filt['Norm_AA_Freq'] = (df_filt['AA_Freq'] - df_filt['WT_AA_Freq']) /  df_filt['AA_Freq'].sum()

	# Create a dataframe to store the last 'Seq_Index' value for each 'Coord' group
	last_index_by_coord = df_filt.groupby('Coord',sort=False)['Seq_Index'].last()

	# Create a new column 'New_Seq_Index' for the updated 'Seq_Index' values
	df_filt['New_Seq_Index'] = df_filt['Seq_Index']

	balance_index = 0
	# For each 'Coord' group after the first group, update the 'New_Seq_Index' values
	for i in range(len(last_index_by_coord)):
		df_filt.loc[df_filt['Coord'] == last_index_by_coord.index[i], 'New_Seq_Index'] = df_filt.loc[df_filt['Coord'] == last_index_by_coord.index[i], 'Seq_Index'] + balance_index
		balance_index = balance_index + last_index_by_coord.iloc[i] + 1

	# Drop the original 'Seq_Index' column and rename 'new_seq_index' to 'Seq_Index'
	df_filt = df_filt.drop(['Coord', 'Coord_Start', 'Seq_Index', 'AA_Freq'], axis=1)

	print(df_filt)
	empty_data = np.zeros((len(amino_acids), int(df_filt['New_Seq_Index'].iat[-1])+1), dtype=float)
	df_index = amino_acids
	df_columns = list(range(int(df_filt['New_Seq_Index'].iat[-1])+1))
	df_final = pd.DataFrame(empty_data, index = df_index, columns = df_columns)

	# # iterate over each residues fraction range 
	# # # new dataframe filled with residues of interest data subset
	wt_amino_acid = []
	wt_check = -1
	for ii in range(df_filt.shape[0]):
		subs_AA = df_filt['Mut_AA'][ii]
		wt_AA = df_filt['WT_AA'][ii]
		subs_Pos = df_filt['New_Seq_Index'][ii]
		subs_value = df_filt['Norm_AA_Freq'][ii]
		# subs_value = df_filt['AA_Freq'][ii]
		if subs_AA == 'X':continue
		df_final.at[subs_AA,subs_Pos]= subs_value
		if subs_Pos != wt_check:
			wt_amino_acid.append(wt_AA)
			wt_check = subs_Pos

	fig,ax = plt.subplots()
	ax.set_xticks((np.arange(0, len(df_columns), 20)))
	ax.tick_params('x',  which='major')

	ax.set_yticks((np.arange(len(df_index))))
	ax.tick_params('y', which='major')
	ax.set_yticklabels(df_index)

	coord_dict = {}
	val = 0
	for res in amino_acids:
		coord_dict[res] = val
		val += 1

	im = plt.imshow(df_final,interpolation='nearest', origin = 'lower',alpha=0.7)
	## mark WT
	for i in range(df_final.shape[1]):
		plt.scatter(i, coord_dict[wt_amino_acid[i]], color='red', s=5)

	figure = plt.gcf() # get current figure
	figure.set_size_inches(len(df_columns)/1.5, len(df_index)/1.5)
	# figure.set_size_inches(20, 10)
	plt.subplots_adjust(bottom=0.4)

	fig.colorbar(im, orientation="horizontal", pad=0.2, shrink = 0.3)

	plt.savefig(output_file[0], dpi = 300,bbox_inches='tight')
	plt.show()

	    # write the result to the output file
	df_final.to_csv(output_file[1], sep='\t')

	tmp = output_file[0]
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="A script to make plot")
	parser.add_argument('-i', nargs='+', dest="input_argument",default="something.txt", help="something file")
	parser.add_argument('-o', nargs='+', dest="output_argument",default="something_.txt", help="something file")
	args_ = parser.parse_args()
	main(  )
	print('done')