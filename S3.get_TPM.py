import numpy as np
import pandas as pd
import csv

#输入上一步过滤后的enhancer信息文件
df1 = pd.read_csv('../merge_peaks/filtered_mm10.csv',header = 0)
#输入基因注释文件（这里根据原本的基因注释文件处理过了，提取了chr,start,end,geneid）
df2 = pd.read_csv('../vM17_geneannotation.csv',header = 0)
gene = df2[['chr','start','end','gene_id']]
gene.set_index('gene_id',inplace=True)
gene["length"] = gene['end']-gene['start']
enh = df1[['chr','start','end','length']]
enh.index = list(range(1,len(enh)+1))

#打开enh,gene count文件并过滤不需要的数据形成表格，这里以XW的数据为例
enh_count0 = pd.read_table("enhancer_XW_counts.txt",index_col=0)
enh_count = enh_count0[enh_count0['gene_id'].str.contains("Enh")] 
gene_count0 = pd.read_table("gene_XW_counts.txt")
gene_count = gene_count0[gene_count0['gene_id'].str.contains("ENSMUSG")]  
gene_count.set_index('gene_id',inplace=True)
#输入count文件中的时期
#timelist_GSR=['2cell', '4cell', '8cell', 'EpiE65', 'ExeE65', 'ICM','MIIOocyte', 'morula', 'TE']
timelist_XW=['E2C','ICM','L2C','M4C','M8C','MIIOocyte','Zygote']
#把位置信息和不同时期count信息合并
gene2 = pd.concat([gene,gene_count],axis=1)
enh2 = pd.concat([enh,enh_count],axis=1)

#用于计算TPM的函数
def cal_TPM(df,timelist):
    df_co = df.copy()
    for i in timelist:
        df_co[i]=df[i]*1000/df['length']
    alstreads=[]
    for i in timelist:
        alstreads.append(sum(df_co[i]))
    k = 0
    df_co2=df_co.copy()
    for i in timelist:
        df_co2[i] = (df_co[i]/alstreads[k])*1000000
        k=k+1
    return df_co2

#计算并输出结果
gene_TPM = cal_TPM(gene2,timelist_XW)
gene_TPM.to_csv("gene_XW_TPM.csv")
enh_TPM = cal_TPM(enh2,timelist_XW)
enh_TPM.to_csv("enh_XW_TPM.csv")
