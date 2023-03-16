
import argparse
import matplotlib.pyplot as plt  
import pandas as pd
# from matplotlib.patches import Circle, Wedge
from matplotlib import rc, rcParams
import numpy as np
from collections import OrderedDict
import os
# Latex features
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})   
rc('axes', labelsize=20)
rc('axes', titlesize=22)   # fontsize of the x and y labels
rc('xtick', labelsize=14)    # fontsize of the tick labels
rc('ytick', labelsize=14)    # fontsize of the tick labels
rc('legend', fontsize=12)    # legend fontsize
# plt.rc('figure', titlesize=25)  # fontsize of the figure title
rcParams['font.size'] = '16'
rcParams['text.latex.preamble'] = r'\usepackage{sfmath}'
rcParams['ps.fonttype'] = 42
rcParams['svg.fonttype'] = 'none' 
# rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

AA_col = {
	'R':'green',
	'C':'red',
	'P':'blue',
	'N':'black'
}

amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def grid_plot(df1_final, feature_col, fig_dir, fig_specific_chunk, file):

	fig_name = ('').join([fig_specific_chunk, file,'.pdf'])
	fig_file = os.path.join(fig_dir, fig_name) 

	row,col = df1_final.shape

	y_label = df1_final.index.values
	y_label[:] = y_label[::-1]
	x_lablel = df1_final.columns.values

	mat =np.ones([row, col, 3], dtype=np.uint8)*255
	for i in range(row):
		for j in range(col):
			try:
				mat[i,j,:]= np.asarray(df1_final.iloc[i,j],dtype=np.uint8)
			except ValueError:
				continue

	leg1_loc = 'lower center'
	fig, ax = plt.subplots(figsize=(15,8))
	ax.set_xticks(0.5 +(np.arange(0,col,5)))
	ax.tick_params('x', labelsize="x-small", which='major')
	ax.set_xticklabels(x_lablel[::5])

	ax.set_yticks(0.5 +(np.arange(row)))
	ax.tick_params('y', labelsize="small", which='major')
	ax.set_yticklabels(y_label)
	plt.imshow(mat,extent=[0,col,0,row])
	
	plt.subplots_adjust(top=0.7)
	plt.tight_layout()
	plt.savefig(fig_file, dpi = 300)
	# show the plot
	plt.show()

def main( ):
	input_file = args_.input_argument
	output_file = args_.output_argument

 	# read foldx avg delta delta G file for all AA
	foldx_data = pd.read_csv(input_file[0], sep = '\t', header = 0)

	fig_dir = output_file[0] # directory for figures output
	fig_specific_chunk = output_file[1]


	# Fake amino acid coordinate dictionary for plotting
	val = 0
	coord_AA_dict = {}
	sorted_keys = sorted(CN_AA, key=CN_AA.get)
	for res in sorted_keys:		
		coord_AA_dict[res] = val
		val += 1

	# iterate over each residues fraction range 

	for file in input_file[1:]:
		new_res_info = file.split('-')
		new_res_start = int(new_res_info[0])-1 # current residues fraction start pos
		new_res_end = int(new_res_info[1]) # residues fraction end pos
		# foldx_data2 = foldx_data.copy()
		df1 = foldx_data[(foldx_data['Position'] > new_res_start) & (foldx_data['Position'] <= new_res_end)]
		df1 = df1.reset_index(drop=True)
		# Empty character matrix to store values
		new_res_len = new_res_end-new_res_start

		# empty_data = np.array((20,new_res_len))
		empty_data = np.empty((20,new_res_len),dtype='str')

		# new dataframe indexed by AA name and Columns by AA Position
		columns1 = range(new_res_start+1,new_res_end+1,1)
		index = sorted_keys
		df1_final = pd.DataFrame(empty_data, index=index, columns=columns1)

		# # new dataframe filled with residues of interest data subset
		for ii in range(df1.shape[0]):
			subs_AA = df1['Mutation'][ii][-1]
			subs_Pos = df1['Position'][ii]
			subs_Feat = df1['Feature'][ii]
			df1_final.at[subs_AA,subs_Pos]= feature_col[subs_Feat]# subs_Feat

		# grid plot for the residues of interest
		grid_plot(df1_final, feature_col, fig_dir, fig_specific_chunk, file)
	tmp = output_file[0]
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="A script to make plot")
	parser.add_argument('-i', nargs='+', dest="input_argument",default="something.txt", help="something file")
	parser.add_argument('-o', nargs='+', dest="output_argument",default="something_.txt", help="something file")
	args_ = parser.parse_args()
	main(  )
	print('done')