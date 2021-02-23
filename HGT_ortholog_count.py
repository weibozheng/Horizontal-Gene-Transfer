infc=open('/home/hzhu/ciliates/Paramecium_tetraurelia/ptet_protein_blastp_total_ciliate_paramecium_free.tab')
infb=open('/home/hzhu/ciliates/Paramecium_tetraurelia/ptet_protein_blastp_total_bac.tab')
outf=open('/home/hzhu/ciliates/Paramecium_tetraurelia/ptet_HGT_gene.tab',"w")
print("id\tbac\tciliate",file=outf)
qd={}
for line in infc:
	ll=line.split("\t")
	qid=ll[0]
	qlike=float(ll[2])
	if qid not in qd:
		qd.update({qid:qlike})
	else:
		if qlike > qd[qid]:
			qd[qid]=qlike
qdb={}
for line in infb:
	ll=line.split("\t")
	qid=ll[0]
	qlike=float(ll[2])
	if qid not in qdb:
		qdb.update({qid:qlike})
	else:
		if qlike > qdb[qid]:
			qdb[qid]=qlike
hgt_count=0
for key in qd:
	if key in qdb:
		if qdb[key]>qd[key] and qdb[key]>50:
			hgt_count+=1
			print("%s\t%f\t%f"%(key,qdb[key],qd[key]),file=outf)
outf.close()
print(hgt_count)