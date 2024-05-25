#!/usr/bin/env python3
###########################################################################
#	matrix_shaper.py
#	version.3: Update for handling multiple sgRNA library.
#	Copyright (c) 2024, Noguchi Yuki (Jun Suzuki lab)
#	This software is released under the MIT License, see LICENSE (https://opensource.org/license/mit/).
#   @citation: Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki., J. 2024. In vivo CRISPR screening directly targeting testicular cells. Cell Genomics.
#   @author:  Noguchi Yuki
#   @contact: nyuhki21@gmail.com,jsuzuki@icems.kyoto-u.ac.jp
#
###########################################################################

import numpy as np
import pandas as pd
import sys
import os
import argparse
import warnings
warnings.simplefilter('ignore', FutureWarning)
pd.options.mode.chained_assignment = None

def parser_setting():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--count_file',
                        action='store',
                        type=argparse.FileType('r'),
                        default='-',
                        help='File path of count file')
	parser.add_argument('-o', '--output_dir', help='output directry')
	args = parser.parse_args()
	return vars(args)

def stats(data):
	df_counted, df_uncounted = data[data['count'] > 0], data[data['count'] == 0]
	df_counted['mapped'], df_uncounted['mapped'] = 1, 0
	df_stats = pd.concat([df_counted, df_uncounted]).groupby('Gene').sum().reset_index()
	df_stats = df_stats.rename(columns={'count':'total_count'})
	df_stats['Log10(CPM+1)'] = np.log10(1 + 1000000 * df_stats['total_count']/df_stats['total_count'].sum())
	data['count'],data['sgRNA']= data['count'].astype(str) + ',',data['sgRNA'].astype(str) + ','
	df_matrix = data.groupby('Gene').sum().reset_index()
	no_col = data.groupby(by="Gene").sum()["count"].str.split(",",expand=True).columns.max()
	df_sgRNA_count = data.groupby(by="Gene").sum()["count"].str.split(",",expand=True).drop([no_col],axis=1).fillna('')
	df_sgRNA_count['num_sgRNAs'] = df_sgRNA_count.apply(lambda row: sum(row != ''), axis=1)
	df_stats = df_stats.merge(df_sgRNA_count.reset_index().loc[:,['Gene','num_sgRNAs']],how="left", on="Gene")
	df_stats["%mapped"] = 100 * df_stats['mapped'] / df_stats['num_sgRNAs']
	return df_stats, df_matrix

def aggr2matrix(matrix,stats):
	df_result = stats.loc[:,['Gene', 'Log10(CPM+1)', 'total_count', 'num_sgRNAs', 'mapped', '%mapped']].merge(matrix, on="Gene")
	no_col = stats['num_sgRNAs'].max()
	count_cols = [f'UID.{n}(count)' for n in np.arange(no_col)+1]
	id_cols = [f'UID.{n}' for n in np.arange(no_col)+1]
	result_cols = ",".join([f'UID.{n},UID.{n}(count)' for n in np.arange(no_col)+1]).split(',')
	df_result[count_cols] = df_result["count"].str.split(',', expand=True).drop([no_col],axis=1)
	df_result[id_cols] = df_result["sgRNA"].str.split(',', expand=True).drop([no_col],axis=1)
	df_result = df_result.loc[:,
		['Gene', 'Log10(CPM+1)', 'total_count', 'num_sgRNAs', 'mapped', '%mapped']
		+ result_cols
		]
	df_result = df_result.sort_values(by=["total_count"], ascending=False)
	df_result = df_result.sort_values(by=["%mapped"], ascending=False)
	return df_result

if __name__ == "__main__":
	print('start generating sgRNA count matrix..')
	args = parser_setting()
	"""formatting"""
	df = pd.read_table(args['count_file'], sep = '\t').rename(columns = {'sample1':'count'})
	"""stats"""
	df_stats, df_matrix = stats(data=df)
	"""finalization"""
	df_result = aggr2matrix(matrix=df_matrix,stats=df_stats)
	df_result.to_csv(args["output_dir"] + "_result.csv")
	print("matrix has been generated!")
	
