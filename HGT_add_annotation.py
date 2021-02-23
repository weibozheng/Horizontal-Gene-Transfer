import re,sys,os
import pandas as pd
target_gff=open(sys.argv[1])
modfile=open(sys.argv[2])
outfile=open(sys.argv[3],"w")
dta=pd.read_csv('/home/hzhu/ciliates/Paramecium_tetraurelia/ptetraurelia_mac_annotation_v1.tab',sep="\t")
def region_belong(tttdict,tag,ostart,oend):
	for tstart in tttdict[tag]:
		tstr=tttdict[tag][tstart]
		mobjx=re.match("(.*?)_(.*)",tstr)
		tend=int(mobjx.group(1))
		gid=mobjx.group(2)
		#print(gid)
		ol=oend-ostart+1
		overlap=0
		if ostart<tstart and oend>tstart and oend<tend:
			overlap=oend-tstart+1
		if ostart<tstart and tend<oend:
			overlap=tend-tstart+1
		if tstart<ostart and tend>ostart and tend<oend:
			overlap=tend-ostart+1
		if tstart<ostart and oend<tend:
			overlap=oend-ostart+1
		if overlap!=0 and overlap/ol>0.5:
			return True,gid
	return False,False
def build_dict(tdict,ttag,tstart,tend,tgid):
	if ttag not in tdict:
		tdict.update({ttag:{int(tstart):tend+'_'+tgid}})
	else:
		tdict[ttag].update({int(tstart):tend+'_'+tgid})
gene_dict=dict()
CDS_dict=dict()
des_dict=dict()
for i in dta.index:
	desc_str=dta.loc[i,'CROSS_REFERENCES']+"\t"+str(dta.loc[i,'DESCRIPTION'])
	tid=dta.loc[i,'ID']
	des_dict.update({tid:desc_str})
for line in target_gff:
	if line.count('#')==0:
		ll=line.strip().split("\t")
		chr_tag=ll[0]
		ttype=ll[2]
		start=ll[3]
		end=ll[4]
		att=ll[8]
		if ttype=='gene':
			mobj=re.search("Name=(.*)",att)
			gid=mobj.group(1)
			build_dict(gene_dict,chr_tag,start,end,gid)
		if ttype=='CDS':
			mobj=re.search("Parent=(.*)",att)
			gid=mobj.group(1)
			build_dict(CDS_dict,chr_tag,start,end,gid)
		
for line in modfile:
	ll=line.split("\t")
	target_tag=ll[0]
	start=ll[1]
	end=ll[2]
	this_type='noncoding'
	if target_tag in gene_dict:
		(x,y)=region_belong(gene_dict,target_tag,int(start),int(end))
		if x:
			att=y
			(a,b)=region_belong(CDS_dict,target_tag,int(start),int(end))
			if a:
				this_type='CDS'
			else:
				this_type='intron'
		else:
			this_type='noncoding'
			att='None'
	if att in des_dict:
		print(line.strip()+"\t"+this_type+"\t"+att+"\t"+des_dict[att],file=outfile)
	else:
		print(line.strip()+"\t"+this_type+"\t"+att+"\t"+'None'+'\t'+'None',file=outfile)