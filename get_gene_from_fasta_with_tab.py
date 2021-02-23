from Bio import SeqIO
import pandas as pd
df=pd.read_csv("/home/hzhu/ciliates/Paramecium_tetraurelia/ptet_HGT_gene.tab",sep="\t")
pdict={}
for i in df.index:
	pdict.update({df.iloc[i]['id']:1})
ll=[]
for rec in SeqIO.parse("/home/hzhu/ciliates/Paramecium_tetraurelia/ptetraurelia_mac_annotation_v1.protein.fa","fasta"):
	if str(rec.id) in pdict:
		ll.append(rec)
SeqIO.write(ll,"HGT_gene.fasta","fasta")