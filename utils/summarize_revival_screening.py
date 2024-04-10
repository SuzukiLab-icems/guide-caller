#!/usr/bin/env python3
###########################################################################
#guide_caller/utils/summarize_revival_screening.py
#
#	 Copyright (c) 2024, Noguchi Yuki (Jun Suzuki lab)
#	 This software is released under the MIT License, see LICENSE (https://opensource.org/license/mit/).
#    @citation: Noguchi, Y., Maruoka, M., Suzuki, J. 2024. STAR Protocols.
#    @author:  Noguchi Yuki
#    @contact: nyuhki21@gmail.com,jsuzuki@icems.kyoto-u.ac.jp
#
###########################################################################

import glob
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import warnings
warnings.simplefilter('ignore', FutureWarning)
pd.options.mode.chained_assignment = None

print('ex. python ./utils/summarize_revival_screening.py -m meta_data_for_revival_screening.csv')
print('ex. python ./utils/summarize_revival_screening.py -m meta_data_for_revival_screening.csv -n 25 -f False')

def parser_setting():
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--meta_data', type=str, help='Specify your meta data (csv format).')
	parser.add_argument('-d', '--your_directory', type=str, help='Specify your directory.')
	parser.add_argument('-n', '--number_of_candidates', type=int, help='Specify the number of candidates for plot (Max=180, Default=45).')
	parser.add_argument('-f', '--filter_low_count_sgRNA', type=str, help='Filter low count sgRNAs (Default: True)')
	args = parser.parse_args()
	return vars(args)
	
def _to_log2_cpm_plus_one(data):
	return np.log2(1 + 1000000 * data['Count'] / data['Count'].sum())

def integrate_gene_count_file(count_files, meta_data):
	df_init = pd.DataFrame(columns = [type,'Count','Log2(CPM+1)','Pool','Round'])
	for f in count_files:
		sample_id = f.split('/')[-1].split('.')[0]
		tmp_df = pd.read_table(f,sep='\t').rename(columns = {'sample1':'Count'})
		tmp_df = tmp_df.groupby('Gene').sum().reset_index()
		tmp_df['Pool'] = meta_data[meta_data['sample_id'] == sample_id]['Pool'].iloc[0]
		tmp_df['Round'] = meta_data[meta_data['sample_id'] == sample_id]['Round'].iloc[0]
		tmp_df['Log2(CPM+1)'] = _to_log2_cpm_plus_one(tmp_df) #Normalization
		df_init = pd.concat([df_init,tmp_df])
	return df_init

def _to_enrichment_score(data):
	"""About EnrichmentScore　deriviation
	In the following sgRNA count matrix,

		A(sample)		B(input)
	1	A1				B1
	2	A2				B2

	We calculate ES according to the following processes:
	a. modified FoldChange(mFC) = CPM(A1) / (CPM(B1) + 1)
	b. Enrichment Score(ES) = log2(1 + mFC)
	"""
	data['ES(Before enrichment vs. Input)'] = np.log2(
		1
		+ 1000000 * (data['1'] / data['1'].sum()) #CPM(A1)
		/ (1 + 1000000 * (data["Input"] / data["Input"].sum())) #CPM(B1) + 1
	)
	data['ES(After enrichment vs. Input)'] = np.log2(
		1
		+ 1000000 * (data['2'] / data['2'].sum()) #CPM(A1)
		/ (1 + 1000000 * (data["Input"] / data["Input"].sum())) #CPM(B1) + 1
	)
	data['ES(After enrichment vs. Before enrichment)'] = np.log2(
		1
		+ 1000000 * (data['2'] / data['2'].sum()) #CPM(A1)
		/ (1 + 1000000 * (data['1'] / data['1'].sum())) #CPM(B1) + 1
	)
	return data

def visualize_enriched_candidates(data,pool,out_dir):
	'''Configuration'''
	sns.set(
		style="ticks",
		font_scale=2.5,
		font='arial',
		rc = {'figure.figsize':(12,10),'axes.linewidth':2}
	)
	font_size = 32
	'''Plot'''
	df_enriched = pd.DataFrame(columns=['Rank','ES','Comparison'])
	df_summary = pd.DataFrame(columns=['Rank','ES','Comparison'])
	for c in data.columns:
		out_file = out_dir + '/Enriched Genes in ' + c.split('(')[1].split(')')[0] + ' from Pool' + pool + '.png'
		d = data.loc[:,[c]]
		d['Rank'] = d[c].rank(method='first', ascending=False, na_option='bottom')
		#ES = 2 approximatelly equals to Log2FC = 1.5
		high, low = d[d[c]>=2], d[d[c]<2]
		print(f'Pool{pool}: {len(high.index)} Genes were enriched in {c}')
		sns.scatterplot(data=high, x="Rank", y=c, color="orchid", s=200, edgecolor='orchid', alpha=0.5)
		sns.scatterplot(data=low, x="Rank", y=c, color="whitesmoke", s=75, edgecolor='silver', alpha=0.5)
		plt.title('Enriched Genes: ' + str(len(high.index)) + '\n' + c.split("(")[1].split(")")[0])
		plt.xlabel("Enrichment Rank" ,fontsize=font_size)
		plt.ylabel("Entichment Score",fontsize=font_size)
		plt.xticks(fontsize=font_size,rotation=0)
		plt.yticks(fontsize=font_size,rotation=0)
		plt.tight_layout()
		sns.despine()
		plt.axhline(y=2.0, linewidth=1.5, c='grey', ls='--')
		plt.savefig(out_file, dpi=300)
		plt.close()
		#Saving the results..
		d['ES'], d['Comparison'] = d[c], c
		df_summary = pd.concat([df_summary,d.loc[:,['Rank','ES','Comparison']]])
	df_summary.to_csv(out_dir + '/Summary_of_Enrichment_Score_in_Pool_' + pool + '.csv')
	return df_summary

def visualize_top_candidates(data, pool, out_dir, num_candidates):
	df_plot = data.loc[data.head(n=num_candidates).index.tolist()].reset_index().rename(columns={'index':'Gene'})
	'''Configuration'''
	sns.set(
		style="ticks",
		font_scale=2.5,
		font='arial',
		rc = {'figure.figsize':(int(0.6*num_candidates),int(0.3*num_candidates)),'axes.linewidth':2}
		)
	plt.rcParams['figure.subplot.bottom'] = 0.5
	plt.rcParams['figure.subplot.right'] = 0.9
	plt.rcParams['figure.subplot.left'] = 0.05
	plt.rcParams['figure.subplot.top'] = 0.9
	color_dict = dict({1:'orchid',0:'thistle'})
	font_size = 36
	kwargs = {"edgecolor":"black", "linewidth":1.5}
	'''Plot'''
	#Plot
	fig,ax = plt.subplots()
	ax = sns.barplot(
		data=df_plot,
		x="Gene", y="Log2(CPM+1)",
		hue="Enrichment",
		palette=["royalblue", "midnightblue"],
		linewidth=3,
		edgecolor="black",
		zorder=0)
	ax2 = ax.twinx()
	ax2.axhline(c='mediumvioletred', lw=3, ls='--', alpha = 0.6, zorder=1)
	sns.scatterplot(
		data=df_plot,
		x="Gene", y="Log2FoldChange",
		palette=color_dict,
		hue="Enriched",
		size="Log2FoldChange",
		sizes=(250,1000),
		legend=False,
		**kwargs,
		zorder=2,
		ax=ax2)
	labels = df_plot[df_plot['Enrichment'] == 'Before']
	ax.set_ylabel('Log2(CPM+1)', fontsize=font_size)
	ax2.set_ylabel('Log2FoldChange', fontsize=font_size, rotation=270, labelpad=40)
	ax.set_xticklabels(labels=labels["Gene"], rotation=315, ha ='left')
	ax.set_ylim([1.05*min(df_plot['Log2(CPM+1)']),1.05*max(df_plot['Log2(CPM+1)'])])
	plt.ylim([0.8*min(df_plot['Log2FoldChange']),1.2*max(df_plot['Log2FoldChange'])])
	ax.legend(
		bbox_to_anchor=(1.02, 1.15),
		loc='upper left',
		frameon=False,
		fontsize=font_size,
		fancybox=False,
		edgecolor="black")
	plt.title(f'Top {num_candidates} genes',fontsize=font_size,y=1.02)
	plt.tight_layout()
	plt.savefig(f'{out_dir}/Top_Candidates_in_Pool_{pool}.png',format='png',dpi=300)
	plt.close()

def enrichment_analysis(data_matrix, out_dir, num_candidates, filter_low_count_sgRNA):
	for p in np.unique(data_matrix['Pool']).tolist():
		#1.For enrichment plot
		dtx = data_matrix[data_matrix['Pool'] == p]\
			.pivot_table(
				values='Count',
				columns=['Round'],
				index=['Gene'],
				aggfunc=np.sum #Apply total count for analysis
			).astype('int')
		dtx = _to_enrichment_score(dtx)
		summary_table = visualize_enriched_candidates(
			data = dtx.loc[:,dtx.columns[::-1][:3]],
			pool = p,
			out_dir=out_dir
		)
		#2.For extracting top genes
		#Extract enriched genes in `Before Enrichment` screening. In Noguchi et al. 2024. Cell Genomics, Log2FC=1.5 (about ES=2.0) was applied for extracting candidate genes.
		if filter_low_count_sgRNA == True: candidate_genes = dtx[dtx['ES(Before enrichment vs. Input)'] > 2].index.tolist()
		else: candidate_genes = dtx.index.tolist()
		dtx = data_matrix[data_matrix['Pool'] == p]\
			.pivot_table(
				values='Log2(CPM+1)',
				columns=['Round'],
				index=['Gene'],
				aggfunc=np.sum #Apply total count for analysis
			).loc[candidate_genes].loc[:,['1','2']]
		dtx.columns = ['Before','After']
		dtx['Log2FoldChange'] = dtx['After'] - dtx['Before']
		dtx = dtx.sort_values(by='Log2FoldChange',ascending=False)
		dtx.to_csv(out_dir + '/Summary_of_enriched_genes_in_Pool_' + p + '.csv')
		summary = pd.DataFrame(columns=['Enrichment','Log2(CPM+1)','Log2FoldChange'])
		for c in ['Before','After']:
			tmp = dtx.loc[:,[c,'Log2FoldChange']]
			tmp['Enrichment'] = c
			tmp = tmp.rename(columns = {c:'Log2(CPM+1)'})
			summary = pd.concat([summary,tmp])
		summary['Enriched'] = summary['Log2FoldChange'].apply(lambda x: 1 if x > 0 else 0)
		visualize_top_candidates(
			data = summary,
			pool = p,
			out_dir = out_dir,
			num_candidates = num_candidates
		)
		
def main(meta_data, your_directory, number_of_candidates, filter_low_count_sgRNA):
	if number_of_candidates is None: number_of_candidates = int(45)
	if filter_low_count_sgRNA is None: filter_low_count_sgRNA = True
	meta_data = pd.read_csv(f'./{your_directory}/{meta_data}').dropna(axis=0)
	count_files = glob.glob(f'./{your_directory}/*/*/*.count.txt')
	data = integrate_gene_count_file(count_files=count_files, meta_data=meta_data)
	enrichment_analysis(data_matrix=data, out_dir=your_directory, num_candidates=number_of_candidates, filter_low_count_sgRNA=filter_low_count_sgRNA)

def execute():
    main(**parser_setting())

if __name__ == "__main__":
	execute()
