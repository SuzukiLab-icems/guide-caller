#!/usr/bin/env python3
###########################################################################
#guide_caller/utils/sumamrize_clone_coverage.py
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

def parser_setting():
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--meta_data', help='Specify your meta data (csv format).')
	parser.add_argument('-d', '--your_directory', help='Specify your directory.')
	args = parser.parse_args()
	return vars(args)
	
def _to_log10_cpm_plus_one(data):
	data['log10(CPM+1)'] = np.log1p(1000000 * data['count'] / data['count'].sum())
	return data

def preprocessing(matrix_files,meta_data):
	sgRNA_init = pd.DataFrame(columns = ['sgRNA','count','log10(CPM+1)','day','number'])
	Gene_init = pd.DataFrame(columns = ['Gene','count','log10(CPM+1)','day','number'])
	for f in matrix_files:
		day = meta_data[meta_data['sample_id'] == f.split('/')[2]]['treatment(day)'].iloc[0]
		num = meta_data[meta_data['sample_id'] == f.split('/')[2]]['treatment(number)'].iloc[0]
		data = pd.read_table(f, sep='\t').rename(columns={'sample1':'count'})
		#for aggregating sgRNA
		df_sgRNA = _to_log10_cpm_plus_one(data.loc[:,['sgRNA','count']].groupby('sgRNA').mean()).reset_index()
		df_sgRNA['day'], df_sgRNA['number'] = day, num
		#for aggregating Gene
		df_Gene = _to_log10_cpm_plus_one(data.loc[:,['Gene','count']].groupby('Gene').mean()).reset_index()
		df_Gene['day'], df_Gene['number'] = day, num
		#finalize
		sgRNA_init = pd.concat([sgRNA_init,df_sgRNA[df_sgRNA['log10(CPM+1)']>0]],axis=0)
		Gene_init = pd.concat([Gene_init,df_Gene[df_Gene['log10(CPM+1)']>0]],axis=0)
		sgRNA_init['Type'], Gene_init['Type'] = 'sgRNA', 'Gene'
		analysis_table = pd.concat([
			sgRNA_init.rename(columns={'sgRNA':'id'}),
			Gene_init.rename(columns={'Gene':'id'})]
			)
	return analysis_table

def cumulative_plot(data,out_dir):
	#configuration
	sns.set(style="ticks", font_scale=2.5, font='arial',
		rc = {'figure.figsize':(12.5,8), 'axes.linewidth':1.5})
	sns.color_palette("viridis_r", as_cmap=True)
	#ecdfplot
	sns.ecdfplot(data=data,
		palette='viridis_r',
		x='log10(CPM+1)',
		hue='number',
		lw=3,
		stat='count',
		legend=False)
	plt.tick_params(labelsize = 30, width = 1.5 )
	plt.ylim([0,np.round(1.05*len(np.unique(data['id']).tolist()))])
	plt.xlim([1.2*min(data['log10(CPM+1)']),1.2*max(data['log10(CPM+1)'])])
	plt.xlabel("log10(CPM+1)", fontsize=30)
	plt.ylabel("count", fontsize=30)
	plt.legend(labels=["INPUT", "9 Testis","3 Testis"],
		bbox_to_anchor=(1.02, 1),
		loc='upper left',
		borderaxespad=0,
		frameon=False,
		fontsize=30,
		fancybox=False,
		edgecolor="black")
	plt.tight_layout()
	sns.despine()
	plt.savefig(out_dir,dpi=300)
	plt.close()

def evaluate_library_coverage(your_directory, meta_data,analysis_table):
	for t in ['sgRNA','Gene']:
		data = analysis_table[analysis_table['Type'] == t]
		day = np.unique(data['day']).tolist()
		for d in day:
			if d == 'Input': continue
			out_dir = your_directory + '/' + t + '_cumulative_plot_in_' + d + '.png'
			dtx = data[(data['day'] == d) | (data['day'] == 'Input')]
			cumulative_plot(data=dtx, out_dir=out_dir)
			print(f'[SUMMARY: {t} clone coverage in {d}]')
			print(f'{d}:\n{dtx.loc[:,["number","id"]].groupby("number").count()}')
			print('')

def visulize_sgRNA_multiplicity(data,sample_name,out_dir):
	df_count = data.groupby(by="Gene").count().rename(columns={"sample1":"count"})
	count = df_count.reset_index().groupby(by="count").count()
	print(f'[Number of genes targeted by multiple sgRNAs in {sample_name}]')
	print(count.reset_index().rename(columns={'count':'Number of sgRNAs','Gene':'Number of targeted genes'}).set_index('Number of sgRNAs'))
	print('')
	sns.countplot(data=df_count, x="count", palette="viridis_r")
	plt.xlabel("sgRNA Counts per Gene" ,fontsize=38)
	plt.ylabel("Number of Genes",fontsize=38)
	plt.xticks(fontsize=38,rotation=0)
	plt.yticks(fontsize=38,rotation=0)
	plt.tight_layout()
	plt.savefig(f'{out_dir}/{sample_name}/{sample_name}_sgRNA_multiplicity.png', dpi=300)
	plt.close()
	
def visulize_targeted_gene_distribution(data,sample_name,out_dir):
	df_count = data.groupby(by="Gene").count().rename(columns={"sample1":"count"})
	df_sum = data.groupby(by="Gene").sum().rename(columns={"sample1":"total_count"})
	df_dist = df_sum.merge(df_count.loc[:,["count"]].reset_index(),how="left",on="Gene")
	df_dist["Log10(CPM+1)"] = np.log1p(1000000 * df_dist["total_count"] / df_dist["total_count"].sum())
	sns.histplot(data=df_dist,x="Log10(CPM+1)",hue="count", palette="viridis_r")
	plt.xlabel("Log10(CPM+1)" ,fontsize=38)
	plt.ylabel("Number of Genes",fontsize=38)
	plt.xticks(fontsize=38,rotation=0)
	plt.yticks(fontsize=38,rotation=0)
	plt.tight_layout()
	plt.axvline(x=2, linewidth=1.5, c='mediumorchid')
	plt.savefig(f'{out_dir}/{sample_name}/{sample_name}_count_histgram.png', dpi=300)
	plt.close()

def evaluate_library_quality(your_directory, matrix_files):
	for file in matrix_files:
		sample=file.split("/")[2]
		df = pd.read_table(file, index_col=0)
		df = df[df["sample1"]>0]
		#sgRNA multiplicity against targeted gene
		visulize_sgRNA_multiplicity(data=df,sample_name=sample, out_dir=your_directory)
		#distribution
		visulize_targeted_gene_distribution(data=df,sample_name=sample, out_dir=your_directory)

def main(meta_data,your_directory):
	meta_data = pd.read_csv(f'./{your_directory}/{meta_data}')
	matrix_files = glob.glob(f'./{your_directory}/*/*/*.count.txt')
	#summarize_clone_coverage
	analysis_table = preprocessing(matrix_files = matrix_files, meta_data = meta_data)
	evaluate_library_coverage(your_directory = your_directory, meta_data = meta_data, analysis_table = analysis_table)
	#evaluate_distribution
	evaluate_library_quality(your_directory = your_directory, matrix_files = matrix_files)
	#summarization..
	analysis_table[analysis_table['Type'] == 'sgRNA'].pivot_table(values='log10(CPM+1)',columns=['day','number'],index=['id']).fillna('undetected').to_csv(your_directory + '/sgRNA_count_summary.csv')
	analysis_table[analysis_table['Type'] == 'Gene'].pivot_table(values='log10(CPM+1)',columns=['day','number'],index=['id']).fillna('undetected').to_csv(your_directory + '/Gene_count_summary.csv')

def execute():
    main(**parser_setting())

if __name__ == "__main__":
	execute()
