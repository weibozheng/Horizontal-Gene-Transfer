from Bio import SeqIO
import re
import sys
import os

fafile=sys.argv[1]
outfolder=sys.argv[2]
# mobj=re.match(".*/(.*)",fafile)
# basename=mobj.group(1)
outfile=open(outfolder+"/shreddered.fq","w")

phase=10
for rec in SeqIO.parse(fafile,"fasta"):
	this_id=str(rec.id)
	this_seq=str(rec.seq)
	this_length=len(this_seq)
	tear_num=this_length//150
	for j in range(tear_num):
		for i in range(150//phase):
			start=j*150+i*phase
			end=(j+1)*150+i*phase
			this_seg=this_seq[j*150+i*phase:(j+1)*150+i*phase]
			this_qual=["F"]*len(this_seg)
			print("@"+this_id+"_"+str(start+1)+'_'+str(end)+"\n"+this_seg+"\n"+"+"+"\n"+''.join(this_qual),file=outfile)
outfile.close()
