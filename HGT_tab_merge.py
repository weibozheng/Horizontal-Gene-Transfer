import re,sys
def try_merge_line(il,al):
	ilist=il.split("\t")
	alist=al.split("\t")
	itag=ilist[0]
	atag=alist[0]
	istart=int(ilist[1])
	iend=int(ilist[2])
	astart=int(alist[1])
	aend=int(alist[2])
	b_istart=int(ilist[4])
	b_iend=int(ilist[5])
	b_astart=int(alist[4])
	b_aend=int(alist[5])
	
	if itag==atag:
		if (not ( istart>aend or astart>iend )) and (not ( b_istart>b_aend or b_astart>b_iend )):
			if astart<istart:
				istart=astart
			if aend>iend:
				iend=aend
			if b_astart<b_istart:
				b_istart=b_astart
			if b_aend>b_iend:
				b_iend=b_aend
			return(itag+"\t"+str(istart)+"\t"+str(iend)+"\t"+ilist[3]+"\t"+str(b_istart)+"\t"+str(b_iend)+"\t"+ilist[6]+"\t"+ilist[7])
		else:
			return('')
	else:
		return('')
inf=open(sys.argv[1])
outfile=open(sys.argv[2],"w")
dict_all=dict()
for line in inf:
	line_list=line.split("\t")
	line_tag=line_list[3]
	if line_tag not in dict_all:
		dict_all.update({line_tag:[line[:-1]]})
	else:
		dict_all[line_tag].append(line[:-1])
for key in dict_all:
	this_ll=dict_all[key]
	this_count=0
	while this_count<len(this_ll)-1:
		this_line=this_ll[this_count]
		skip_this_seq=1
		for i in range(this_count+1,len(this_ll)):
			that_line=this_ll[i]
			merge_result=try_merge_line(this_line,that_line)
			if merge_result!='':
				this_ll[this_count]=merge_result
				this_ll.pop(i)
				skip_this_seq=0
				break
		if skip_this_seq==1:
			this_count+=1
	dict_all[key]=this_ll
count_l=[]
count_spe=0
for key in dict_all:
	for line in dict_all[key]:
		this_spe=line.split("\t")[6]
		print(line,file=outfile)
		if this_spe not in count_l:
			count_l.append(this_spe)
			count_spe+=1
print(count_spe)