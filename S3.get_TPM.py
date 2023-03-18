import numpy as np
import pandas as pd
import csv

#Input: the filtered enhancer files
df1 = pd.read_csv('../merge_peaks/filtered_mm10.csv',header = 0)
#Input: gene annotation files
df2 = pd.read_csv('../vM17_geneannotation.csv',header = 0)
gene = df2[['chr','start','end','gene_id']]
gene.set_index('gene_id',inplace=True)
gene["length"] = gene['end']-gene['start']
enh = df1[['chr','start','end','length']]
enh.index = list(range(1,len(enh)+1))

#Open and filter count files
enh_count0 = pd.read_table("enhancer_XW_counts.txt",index_col=0)
enh_count = enh_count0[enh_count0['gene_id'].str.contains("Enh")] 
gene_count0 = pd.read_table("gene_XW_counts.txt")
gene_count = gene_count0[gene_count0['gene_id'].str.contains("ENSMUSG")]  
gene_count.set_index('gene_id',inplace=True)
#Get period columns according to the early embryonic development
#timelist_GSR=['2cell', '4cell', '8cell', 'EpiE65', 'ExeE65', 'ICM','MIIOocyte', 'morula', 'TE']
timelist_XW=['E2C','ICM','L2C','M4C','M8C','MIIOocyte','Zygote']
#merge gene count files and enhancer count files accroding to the development time period
gene2 = pd.concat([gene,gene_count],axis=1)
enh2 = pd.concat([enh,enh_count],axis=1)

#Calculte TPM
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

#Output results
gene_TPM = cal_TPM(gene2,timelist_XW)
gene_TPM.to_csv("gene_XW_TPM.csv")
enh_TPM = cal_TPM(enh2,timelist_XW)
enh_TPM.to_csv("enh_XW_TPM.csv")
