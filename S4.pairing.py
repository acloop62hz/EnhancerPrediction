import pandas as pd
import numpy as np
import csv


#Input: the TPM files
gene = pd.read_csv('gene_GSR_TPM.csv',header = 0)
listchr=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
timelist=list(gene.columns[5:])
gene['de']=0
#delete the genes with no expression
for k,row in gene.iterrows():
    if gene.loc[k,timelist].any()==0:
        gene.loc[k,'de'] = -1
gene_clean=gene[(gene['de']==0)&(gene['chr'].isin(listchr))]
gene_clean.drop('de',axis=1,inplace=True)
gene_clean.to_csv("gene_GSR_clean.csv",index=False)

enh = pd.read_csv('enh_GSR_TPM.csv',header = 0,index_col=0)
enh['de']=0
for k,row in enh.iterrows():
    if enh.loc[k,timelist].any()==0:
        enh.loc[k,'de'] = -1
enh_clean=enh[(enh['de']==0)&(enh['chr'].isin(listchr))]
enh_clean.drop('de',axis=1,inplace=True)
enh_clean.set_index('chr',inplace=True)

#go over gene_GSR_TPM files，find enhancers that located +-1Mbp of the gene location
with open('gene_GSR_clean.csv', 'r') as f:
    lines1 = f.read().splitlines()
#output enh-gene pairs files
f4 = open('gene_1Mbp_enh_GSR.csv','w',newline='')
writer = csv.writer(f4)
writer.writerow(['chr','gene_id','enh_id','corr'])

chr0='fdaff'
for genel in lines1:
    if "gene_id" in genel:
        continue
    gene=genel.split(',')
    gid=gene[0]
    chrg=gene[1]
    gs=int(gene[2])
    ge=int(gene[3])
    gene_tpm=list(map(float, gene[5:]))
    # Set chr as indices for faster indexing
    if chrg!=chr0:
        dft=enh_clean.loc[chrg,]
        dft.to_csv(chrg)
        chr0=chrg
    # Go over the enhancers in the chr and fetch enhancers that located +-1Mbp of the gene position
    with open(chrg, 'r') as f2:
        linese = f2.read().splitlines()
    for enhl in linese:
        if "start" in enhl:
            continue
        enh=enhl.split(',')
        if (int(enh[1])>=gs-1000000) and (int(enh[2])<=ge+1000000):
            enh_tpm=list(map(float, enh[5:]))
            data=pd.DataFrame({'e':enh_tpm,'g':gene_tpm})
            datar=data.corr()
            if(np.isnan(datar.iloc[0,1])) :
                r=-1
            else:
                r=datar.iloc[0,1]
            eid=enh[4]
            writer.writerow([chrg,gid,eid,r])
f4.close()