import pandas as pd
import numpy as np
import csv


#输入上一步合并好的peaks
fbigmer = pd.read_csv('merged.csv',header = 0)
#输入基因注释文件（已经去掉gff文件中前缀的）
geneanno = pd.read_table('../gencode.vM17.annotation_table.txt',header = None)
geneanno_clean = geneanno[[0, 2, 3,4,8]].copy()
geneanno_clean.columns=['chr','type','start','end','anno']
ganno_clean = geneanno_clean.loc[geneanno_clean['type'] == 'gene']
ganno_sorted = ganno_clean.sort_values(by=['chr', 'start'], ascending=True)

#坐标缩小10000倍，用于检索时候快速进行初步定位
fbigmer_co = fbigmer.copy()
fbigmer_co['hashs'] = fbigmer_co['start']//10000
fbigmer_co['hashe'] = fbigmer_co['end']//10000
#遍历基因注释文件，将基因坐标上下1000bp的enhancer标注
chrlist = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
for chri in chrlist:
    ganno_temp = ganno_sorted.loc[ganno_sorted['chr'] == chri ]
    bigmer_temp = fbigmer_co.loc[fbigmer_co['chr'] == chri ]
    print(chri)
    for i,row in ganno_temp.iterrows():
        brostart = ganno_temp.loc[i,'start'] - 1000
        broend = ganno_temp.loc[i,'end'] + 1000
        hashs = brostart//10000
        hashe = broend//10000
        fbig_temp2 = bigmer_temp.loc[(bigmer_temp['hashe']+1>=hashs)&(bigmer_temp['hashs']<=hashe)]
        if fbig_temp2.empty == True:continue
        for j,row2 in fbig_temp2.iterrows():
            start = fbig_temp2.loc[j,'start']
            end = fbig_temp2.loc[j,'end']
            if (end <= brostart)|(start >= broend):continue
            else: fbigmer_co.loc[j,'MergeTo'] = 1
#删去标注的基因并输出
filted = fbigmer_co.loc[fbigmer_co['MergeTo'] == -1]
filted.to_csv('filtered_mm10.csv',index=False)