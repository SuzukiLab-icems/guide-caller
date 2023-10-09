##################
#matrix.py [1] [2]
#1: File Directry (HN00171248/F2-1_1/mageck_result)
#2: File Name (F2-1_1)
#Created by NOGUCHI,Yuki
##################

import pandas as pd
import numpy as np
import sys

args = sys.argv
print("Directry=" + args[1])
print("File Name=" + args[2])

df = pd.read_table(args[1] + "/" + args[2] + ".count.txt", sep = '\t')
df = df.rename(columns = { df.columns[2]: 'count' })

df_counted = df[df['count'] > 0]
df_counted['mapped'] = 1
df_uncounted = df[df['count'] == 0]
df_uncounted['mapped'] = 0
df_summary = pd.concat([df_counted, df_uncounted]).groupby('Gene').sum().reset_index()
df_summary = df_summary.rename(columns = { df.columns[2]: 'total_count' })

df['count'] = df['count'].astype(str) + ','
df['sgRNA'] = df['sgRNA'].astype(str) + ','
df_matrix = df.groupby('Gene').sum().reset_index()

df_result = df_summary.merge(df_matrix, on="Gene")
df_result[["UID.1(count)","UID.2(count)","UID.3(count)","UID.4(count)","UID.5(count)","UID.6(count)"]] = df_result["count"].str.split(',', expand=True).drop([6],axis=1)
df_result[["UID.1","UID.2","UID.3","UID.4","UID.5","UID.6"]] = df_result["sgRNA"].str.split(',', expand=True).drop([6],axis=1)
df_result = df_result.loc[:,[
	'Gene', 'total_count', 'mapped',
	'UID.1','UID.1(count)',
	'UID.2','UID.2(count)',
	'UID.3','UID.3(count)',
	'UID.4','UID.4(count)',
	'UID.5','UID.5(count)',
	'UID.6','UID.6(count)'
	]]
df_result = df_result.sort_values(by=["total_count"], ascending=False)
df_result = df_result.sort_values(by=["mapped"], ascending=False)
df_result.to_csv(args[1] + "/" + args[2] + "_result.csv")
