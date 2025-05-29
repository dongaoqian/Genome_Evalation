#!/usr/bin/env python
from collections import Counter
import argparse

parser = argparse.ArgumentParser(description='统计多个agp文件中每个染色体的等位基因分布')
parser.add_argument('-b', '--blast_out', required=True, help='blsatn生成的比对文件')
args = parser.parse_args()

blast_out=args.blast_out

#根据score对windows进行过滤
filter_score={}
filter_score_list={}
with open(blast_out) as data:
    for line in data :
        line=line.strip().split("\t")
        if line[0] in filter_score:
            filter_score[line[0]].append(float(line[-1]))
        else:
            filter_score[line[0]]=[]
            filter_score[line[0]].append(float(line[-1]))#储存window和对应的score

filtered_1_num=0

for i in filter_score:
    filter_score[i].sort()
    if  len(filter_score[i]) > 1 and filter_score[i][-1]/filter_score[i][-2] >= 1.5 and filter_score[i][-1] >=1000:
        filter_score_list[i]=str(filter_score[i][-1])#筛选最佳比对得分大于1000且最佳比对是次优比对的1.5倍的window和score
    else:
        filtered_1_num+=1

#筛选可以用于switch error 检测的reads        
seq_dic={}
seq=[]
with open(blast_out) as data:
    for line in data :
        line=line.strip().split("\t")
        if line[0] in filter_score_list and str(float(line[-1])) == filter_score_list[line[0]]:
            line.insert(0,line[0][:14])
            seq.append(line)
            if line[0] in seq_dic:
                seq_dic[line[0]].append(line[2])
            else:
                seq_dic[line[0]]=[]
                seq_dic[line[0]].append(line[2])#满足条件窗口对应的序列，和窗口比对到的染色体

for seq in seq_dic:
    seq_dic[seq]=Counter(seq_dic[seq])

sum_window=0
sum_filtered_window=0
error_window=0
filterde_2_num=0
for seq in seq_dic:
    items=list(seq_dic[seq].items())
    items.insert(0, [])
    for s in items[1:]:
        items[0].append(int(s[-1]))
    sum_window+=sum(items[0])
    if sum(items[0]) >5 and max(items[0])/sum(items[0]) >= 0.5: #至少有5个有效窗口,且有超过50%窗口映射到同一单倍体
        sum_filtered_window+=sum(items[0])
        if len(items) >2:
           items[0].sort()
           error_window+=sum(items[0][:-1])

print(f"\n\ttotal number of windows :{len(filter_score)}".center(40))
print(f"\tthe number of windows filtered by condition 1 :{filtered_1_num}".center(40))
print(f"\tthe number of windows filtered by condition 2 :{sum_window- sum_filtered_window}".center(40))
print(f"\tthe number of windows available :{sum_filtered_window}".center(40))
print(f"\tthe number of windows with switch_error: {error_window}".center(40))
print("\tswitch_error ratio: {:.3%}\n".format(error_window/sum_filtered_window).center(40))
