import pandas as pd
import numpy as np
import csv

#peaks的初步处理，peaks格式中要有'chr','start','end','length','abs_summit','source'，且为csv或tab格式
def clean_table(pwd):
    listchr=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
    if "csv" in pwd:
        df1_all = pd.read_csv(pwd,header= 0)
    else:
        df1_all = pd.read_table(pwd,header= 0)
    df1_all['source']=pwd
    df1=df1_all[['chr','start','end','length','abs_summit','source']]
    df1_clean = df1.loc[df1['chr'].isin(listchr)]
    df1_without146 = df1_clean.loc[df1_clean['length'] > 146 ]
    return df1_without146
def index(i,df):
    k = df.loc[i,'MergeTo']
    while k != -1:
        i=k
        k= df.loc[i,'MergeTo']
    return i

##pwdlist中输入所有要合并的peaks文件的路径来源
tablelist=[]
pwdlist=["whyte_mm10.csv","EnhancerAtlas2.0_mm10.csv","ATAC_L2C.txt","ATAC_4C.txt","ATAC_8C.txt","ATAC_ICM.txt","DNase_2C.txt","DNase_4C.txt","DNase_8C.txt","DNase_MII-Oocyte.txt","DNase_morula.txt","ATAC_E2C.csv","ATAC_E2C2.csv","H3K27ac_2C.csv","H3K27ac_8C.csv","H3K27ac_Oocytes.csv"]
for i in pwdlist:
    tablelist.append(clean_table(i))
df_without146=pd.concat(tablelist, axis=0, ignore_index=True)
df_sorted = df_without146.sort_values(by=['chr', 'abs_summit'], ascending=True)
df_sorted.insert(6,'MergeTo',np.zeros(len(df_sorted))-1)
df_sorted.set_axis(range(0,len(df_sorted)), inplace=True)

flag = 'chr1'
listchr=[0]
for i in range(0,len(df_sorted)):
    if df_sorted.loc[i,'chr'] != flag:
        flag = df_sorted.loc[i,'chr'] 
        listchr.append(i)
listchr.append(len(df_sorted))

#合并峰之间距离小于100的,并改start和end
for i in range(0,21):
    Start = listchr[i]
    End = listchr[i+1]
    for j in range(Start,End):
        for k in range(j+1,End):
            if abs( df_sorted.loc[k,'abs_summit'] - df_sorted.loc[j,'abs_summit'])<= 100:
                df_sorted.loc[k,'MergeTo'] = index(j,df_sorted)
            else: break
listarray=np.array(df_sorted['MergeTo'])
uniquearray = np.unique(listarray)
uniquelist=uniquearray.tolist()
df_merged=df_sorted[df_sorted['MergeTo'] == -1]
merged_co = df_merged.copy()
for i in range(1,len(uniquelist)):
    df_temp = df_sorted.loc[df_sorted ['MergeTo'] == uniquelist[i]]
    add = df_sorted.loc[uniquelist[i],]
    dicto = add.to_dict()
    df_temp = df_temp.append(dicto,ignore_index=True)
    start = min(df_temp['start'])
    end = max(df_temp['end'])
    alls = list(np.unique(df_temp['source']))
    merged_co.loc[uniquelist[i],'start'] = start
    merged_co.loc[uniquelist[i],'end'] = end
    merged_co.loc[uniquelist[i],'source'] = ";".join(alls)
#输出结果
merged_co.to_csv('merged.csv',index=False)

