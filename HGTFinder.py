import sys,getopt,re,os
import pandas as pd
def build_table_single(this_inf,thisfh,this_pd):
	content_list=this_inf.readlines()
	for j in range(1,len(content_list)-1):
		line=content_list[j]
		line_list=line.split("\t")
		read_tag=line_list[0]
		flag=int(line_list[1])
		tag=line_list[2]
		pos=line_list[3]
		cigar=line_list[5]
		mrnm=line_list[6]
		read_seq=line_list[9]
		read_qual=line_list[10]
		total_l=0
		flist=re.findall("\d+\S",cigar)
		for s in flist:
			mobjx=re.match("(\d+)(\S)",s)
			m1=mobjx.group(1)
			S1=mobjx.group(2)
			if S1!="M":
				total_l+=int(m1)
			else:
				ml=int(m1)
				break
		if ml>=50:
			mobj1=re.match("(.*)_(\d+)_(\d+)",read_tag)
			chr_tag=mobj1.group(1)
			chr_start=int(mobj1.group(2))+total_l
			chr_end=chr_start+ml-1
			spe_name=this_pd.loc[tag,'species_name']
			strain_name=this_pd.loc[tag,'strain_name']
			target_start=int(pos)+total_l
			target_end=target_start+ml-1
			print(chr_tag+"\t"+str(chr_start)+"\t"+str(chr_end)+"\t"+tag+"\t"+str(target_start)+"\t"+str(target_end)+"\t"+spe_name+"\t"+str(strain_name),file=thisfh)
			# try:
				# print(chr_tag+"\t"+str(chr_start)+"\t"+str(chr_end)+"\t"+tag+"\t"+str(target_start)+"\t"+str(target_end)+"\t"+spe_name+"\t"+strain_name,file=thisfh)
			# except:
				# print(type(chr_start),type(chr_end),type(target_start),type(target_end))
				# sys.exit()
help_str='HGTFinder.py -g <genome> -o <output folder>'
genomefile=''
outfolder=''
virus_bwa_ref='/home/hzhu/db/virus_refseq/virus_refseq_idx'
bac_bwa_ref_list=[]
for i in range(10):
	bac_bwa_ref_list.append('/home/hzhu/db/bac_refseq_db/'+'bac'+str(i+1)+'_idx')
try:
	opts,args=getopt.getopt(sys.argv[1:],"hg:o:",["help"])
except getopt.GetoptError:
	print (help_str)
	sys.exit(2)
for opt,value in opts:
	if opt in ("-h","--help"):
		print (help_str)
		sys.exit()
	if opt in ("-o"):
		outfolder=value
	if opt in ("-g"):
		genomefile=value

if not(outfolder and genomefile):
	print (help_str)
	sys.exit()
else:
	virus_outfolder=outfolder+'/virus'
	bac_outfolder=outfolder+'/bac'
	bac_combine_outfolder=bac_outfolder+'/bac_combine'
	readsfile=outfolder+"/shreddered.fq"
	
	os.system('mkdir '+ outfolder)
	os.system('mkdir '+ virus_outfolder)
	os.system('mkdir '+ bac_outfolder)
	os.system('mkdir '+ bac_combine_outfolder)
	os.system('mkdir '+ bac_combine_outfolder+'/virus_reads')
##############reads produce################
if "shreddered.fq" not in os.listdir(outfolder):
	os.system("python ~/tools/HGTFinder/HGT_reads_produce.py "+genomefile+" "+outfolder)
#########################################
virus_desc=virus_bwa_ref+'.des'
vpd=pd.read_csv(virus_desc,sep="\t",index_col=0)
tlistdir=os.listdir(virus_outfolder)
samfile=virus_outfolder+'/virus.sam'
outfile=virus_outfolder+'/virus_HGT.tab'
if "virus.sam" not in tlistdir:
	cmd=r"bwa mem -t 20 "+virus_bwa_ref+" "+readsfile+r" | samtools view -F 4 -o "+samfile
	print(cmd)
	os.system(cmd)
else:
	print("virus.sam already in "+ virus_outfolder+", skip")
inf=open(samfile)
outfh=open(outfile,"w")
build_table_single(inf,outfh,vpd)
inf.close()
outfh.close()
##########################################
combine_count=bac_combine_outfolder+'/bac_HGT_combine.tab'
for i in range(10):
	this_bac_outfolder=bac_outfolder+'/bac'+str(i+1)
	this_sam=this_bac_outfolder+'/bac.sam'
	this_outfile=this_bac_outfolder+'/bac_HGT.tab'
	cpd=pd.read_csv(bac_bwa_ref_list[i]+'.des',sep="\t",index_col=0)
	if not os.path.exists(this_bac_outfolder):
		os.system("mkdir "+this_bac_outfolder)
	if 'bac.sam' not in os.listdir(this_bac_outfolder):
		bac_cmd=r"bwa mem -t 20 "+bac_bwa_ref_list[i]+" "+readsfile+r" | samtools view -F 4 -o "+this_sam
		print(bac_cmd)
		os.system(bac_cmd)
	else:
		print("bac.sam already in "+ this_bac_outfolder+", skip")
	if "bac_HGT.tab" not in this_bac_outfolder:
		inf=open(this_sam)
		outfh=open(this_outfile,"w")
		build_table_single(inf,outfh,cpd)
		inf.close()
		outfh.close()
	else:
		print("bac_HGT.tab already in "+this_bac_outfolder+", skip")
	combine_cmd="cat "+this_outfile+" >> "+combine_count
	print(combine_cmd)
	os.system(combine_cmd)
modified_combine_count=combine_count+'.mod'
merge_cmd="python ~/tools/HGTFinder/HGT_tab_merge.py "+combine_count+" "+modified_combine_count
print(merge_cmd)
os.system(merge_cmd)