#!/usr/bin/env python

import pandas as pd
import sys

file_path = sys.argv[1]
df = pd.read_csv(file_path, sep='\t', header=None, names=['Chr', 'Start', 'End', 'Value'])

df['Rounded_Value'] = df['Value'].round()

mode_value = df['Rounded_Value'].value_counts().idxmax()

filtered_df = df[df['Value'] > 1.5 * mode_value]

total_difference_all = (df['End'] - df['Start']).sum()

total_difference_filtered = (filtered_df['End'] - filtered_df['Start']).sum()

print("CollapseRatio:",total_difference_filtered/total_difference_all)

#Output collapse regions#

base_name = os.path.splitext(file_path)[0]

output_filename = f"{base_name}.collapse.bed" 

filtered_df.to_csv(output_filename, sep='\t', index=False, header=False, columns=['Chr', 'Start', 'End', 'Value'])


