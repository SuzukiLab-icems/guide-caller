#!/usr/bin/env python3
"""
guide-caller matrix_shaper program
Copyright (c) 2023 Yuki NOGUCHI, Jun Suzuki lab
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@version: v2.0.0
@author:  Yuki NOGUCHI
@contact: jsuzuki AT icems.kyoto-u.ac.jp
"""

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
	df_stats['Log2(CPM+1)'] = np.log2(1 + 1000000 * df_stats['total_count']/df_stats['total_count'].sum())
	data['count'],data['sgRNA']= data['count'].astype(str) + ',',data['sgRNA'].astype(str) + ','
	df_matrix = data.groupby('Gene').sum().reset_index()
	df_stats = df_stats.merge(
		data.groupby(by="Gene").sum()["count"].str.split(",",expand=True).drop([6],axis=1).count(axis=1).reset_index().rename(columns={0:"num_sgRNAs"}),
		how="left", on="Gene"
	)
	df_stats["%mapped"] = 100 * df_stats['mapped'] / df_stats['num_sgRNAs']
	return df_stats, df_matrix

def aggr2matrix(matrix,stats):
	df_result = stats.loc[:,['Gene', 'Log2(CPM+1)', 'total_count', 'num_sgRNAs', 'mapped', '%mapped']].merge(matrix, on="Gene")
	df_result[["UID.1(count)","UID.2(count)","UID.3(count)","UID.4(count)","UID.5(count)","UID.6(count)"]] = df_result["count"].str.split(',', expand=True).drop([6],axis=1)
	df_result[["UID.1","UID.2","UID.3","UID.4","UID.5","UID.6"]] = df_result["sgRNA"].str.split(',', expand=True).drop([6],axis=1)
	df_result = df_result.loc[:,[
		'Gene', 'Log2(CPM+1)', 'total_count', 'num_sgRNAs', 'mapped', '%mapped',
		'UID.1','UID.1(count)',
		'UID.2','UID.2(count)',
		'UID.3','UID.3(count)',
		'UID.4','UID.4(count)',
		'UID.5','UID.5(count)',
		'UID.6','UID.6(count)'
		]]
	df_result = df_result.sort_values(by=["total_count"], ascending=False)
	df_result = df_result.sort_values(by=["%mapped"], ascending=False)
	return df_result

if __name__ == "__main__":
	args = parser_setting()
	"""formatting"""
	df = pd.read_table(args['count_file'], sep = '\t').rename(columns = {'sample1':'count'})
	"""stats"""
	df_stats, df_matrix = stats(data=df)
	"""finalization"""
	df_result = aggr2matrix(matrix=df_matrix,stats=df_stats)
	df_result.to_csv(args["output_dir"] + "_result.csv")
	print("matrix has been generated!")
	
