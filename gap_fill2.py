#!/usr/bin/python

import os
import pandas as pd
import re
import collections as cl
import sys
import getopt
import multiprocessing as mul


outputpath=''
refpath=''
querypath='./ngaps/'
opts,args=getopt.getopt(sys.argv[1:],"i:o:r:")
for op, value in opts:
	if op=='-i':
		querypath=value[:-1]+value[-1].replace('/','')+'/'
	if op=='-o':
		outputpath=value[:-1]+value[-1].replace('/','')+'/'
	if op=='-r':
		refpath=value

size0=100
alignment0=2000
insert0=50
match=98
match0=49	

tempfolder=outputpath+'temp00'


try:
	os.mkdir(outputpath+'/')
except:
	pass	


try:
	os.mkdir(tempfolder)
except:
	pass

manager = mul.Manager()


def distractgap(x,nongapcount=9):

	t=x.replace('-','N')
	a1=re.finditer(r'N+',t)
	a1=[a.span() for a in a1]
	if len(a1) == 1:
		return [a1[0][0]],[a1[0][1]]
	if len(a1) == 0:
		return [],[]
	a2=[[a1[i][1]-a1[i][0],a1[i+1][0]-a1[i][1],a1[i+1][1]-a1[i+1][0]] for i in range(0,len(a1)-1)]
	nonignore=[]
	for i in range(len(a2)/2):
		m1=i
		m2=len(a2)-1-i
		s0=(a2[m1][1])*nongapcount
		s1=(a2[m2][1])*nongapcount
		if s0 > a2[m1][0] or s0 > a2[m1][2]:
			nonignore.append(m1)
			if s1 > a2[m2][0] or s1 > a2[m2][2]:
				nonignore.append(m2)
				continue
			a2[m2-1][2]+=a2[m2][2]-s0
			continue
		if s1 > a2[m2][0] or s1 > a2[m2][2]:
			nonignore.append(m2)
			a2[m2-1][2]+=a2[m2][2]-s0
			continue
		a2[m2-1][2]+=a2[m2][2]-s0
		a2[m1+1][0]+=a2[m1][0]-s1
	if len(a2)/2*2 != len(a2):
		m1=len(a2)-len(a2)/2-1
		s0=(a2[m1][1])*nongapcount
		if s0 > a2[m1][0] or s0 > a2[m1][2]:
			nonignore.append(m1)
	nonignore.sort()
	p1=[a1[a][1] for a in nonignore]
	p1.append(a1[-1][1])
	p0=[a1[a+1][0] for a in nonignore]
	p0.insert(0,a1[0][0])
	
	return p0, p1


def blastncheck(query,ref,thetitle,lr):
	
	def align_combine(qstart,qend,rchr,rstd,rstart,rend,identities):
	

		def eachgroup_combine(index,rstart,rend,d):

			if len(index)<2:

				return [index]

			sortindex=sorted(range(len(rstart)),key=lambda x: rstart[x])

			index0=[index[i] for i in sortindex]

			rstart0=[rstart[i] for i in sortindex]

			rend0=[rend[i] for i in sortindex]

			e1=rend0[0]

			combine_index=[]
			combine_index0=[index0[0]]
			for i in xrange(1,len(sortindex)):

				dis1=rstart0[i]-e1

				if dis1<d:

					combine_index0.append(index0[i])


				else:

					combine_index.append(combine_index0)

					combine_index0=[index0[i]]

				e1=max(e1,rend0[i])

			combine_index.append(combine_index0)


			return combine_index

	
		groups=cl.defaultdict(list)
		
		for i,chr0 in enumerate(rchr):
			
			groups[chr0+':'+rstd[i]].append(i)
					
		alignments=[]
		for key in groups.keys():
			
			g_index=groups[key]
			
			g_rstart=[rstart[i] for i in g_index]
			
			g_rend=[rend[i] for i in g_index]
			
			alignments_index0=eachgroup_combine(g_index,g_rstart,g_rend,10000)
			
			
			for each_index in alignments_index0:

				
				each_qstart=[qstart[i] for i in each_index]
				
				each_qend=[qend[i] for i in each_index]
				
				
				each_index_q=eachgroup_combine(range(len(each_qstart)),each_qstart,each_qend,20)
				
				each_qstart=[min([each_qstart[x] for x in  index0]) for index0 in each_index_q]
				
				each_qend=[max([each_qend[x] for x in  index0]) for index0 in each_index_q]

				chr0=key.split(':')[0]
				
				std0=key.split(':')[1]
				
				if std0=='+':
							
					each_rstart=min([rstart[i] for i in each_index])
					
					each_rend=max([rend[i] for i in each_index])
					
				
				else:
					
					each_rstart=max([rstart[i] for i in each_index])
					
					each_rend=min([rend[i] for i in each_index])


				identity=sum([identities[i]*abs(rstart[i]-rend[i]) for i in each_index])/(1+sum([abs(rstart[i]-rend[i]) for i in each_index]))

				identity=float("%.2f"%identity)
					
				size0=sum([abs(each_qend[i]-each_qstart[i]) for i in xrange(len(each_qend))])
				
				alignments.append([chr0,std0,size0,each_qstart,each_qend,each_rstart,each_rend,identity*size0])
		
		return alignments

	outfile=query+'_out'
	#cmd='blastn -query {:s} -db {:s} -outfmt 10 -out {:s} -num_threads 1  -qcov_hsp_perc 80 -perc_identity 85  -max_target_seqs 15'.format(query, ref, outfile)	
	cmd='minimap2 -x ava-pb -splice -m 200 -N 10 -g 1000  --dual=yes -t 3 {:s} {:s} > {:s} 2>&-'.format(query, ref, outfile)
	#print cmd
	
	#os.system(cmd)
	
	header=['query', 'qlength', 'qstart', 'qend', 'strand', 'ref', 'rlength', 'rstart', 'rend', 'nmatch', 'mismatch', 'quality']	

	try:
		t=pd.read_csv(outfile,sep='\s+',names=header)
	except:
		return  []

	if len(t)<1:

		return []


	t=t.loc[(t['match']>200) & (t['query']==thetitle)]

	#t=t.loc[(t['match']>200)]
	
	qstart,qend,rchr,rstart,rend,std=list(t['qstart']),list(t['qend']),list(t['ref']),list(t['rstart']),list(t['rend']),list(t['strand'])

	identities=[1.0*nmatch/mismatch for nmatch,mismatch in zip(list(t['nmatch']),list(t['mismatch']))]
	
	#rstd=['+' if rstart[i]<rend[i] else '-' for i in xrange(len(t))]

	alignments=align_combine(qstart,qend,rchr,rstd,rstart,rend,identities)

	t=pd.DataFrame.from_records(alignments)

	#t=t.loc[t[7]>70]


	if len(t)<1:

		return []	


	results=[]
	for i in xrange(len(t)):


		contig=list(t[0])[i]

		size0=lr

		std=list(t[1])[i]
		
		qstart=min(list(t[3])[i])
		
		qend=max(list(t[4])[i])
		
		score=list(t[7])[i]/size0

		strand=list(t[1])[i]
		
		rstart1=int(list(t[5])[i])
		
		rend1=int(list(t[6])[i])

				
		results.append([contig,strand,rstart1,rend1,qstart,size0-qend,score])
	
	return sorted(results, key=lambda x: x[-1], reverse=False)




def findgap(refpath,tempfolder,filepath,twofiles,sizes,cuts,infor):
	

	frontfind,endfind=[],[]
	
	if filepath[0]!='':
		
		frontfind=blastncheck(filepath[0],refpath,twofiles[0],sizes[0])
		
			
	if filepath[1]!='':
		
		endfind=blastncheck(filepath[1],refpath,twofiles[1],sizes[1])
		
	allsamechre=[]
	for i,frontfind0 in enumerate(frontfind):

		chr0,std,s0,e0=frontfind0[0],frontfind0[1],frontfind0[2],frontfind0[3]
	
		samechre=[a for a in endfind if a[0]==chr0]

		dis=[]
		for samechre0 in samechre:

			cutstart=min(s0,e0,samechre0[2],samechre0[3])

			cutend=max(s0,e0, samechre0[2],samechre0[3])
				
			dis0=cutend-cutstart-abs(s0-e0)-abs(samechre0[2]-samechre0[3])

			dis.append(dis0)
				
		samechreindex=sorted([x for x in xrange(len(samechre)) if dis[x]<100000],key=lambda x : dis[x])

		samechre=sorted([samechre[x] for x in samechreindex],key=lambda x: x[6],reverse=True)

		if len(samechre)==0:

			continue

		bothscore=frontfind0[6]+samechre[0][6]

		allsamechre.append([frontfind0,samechre[0],bothscore])
		
	if len(allsamechre)>0:
		
		allsamechre=sorted(allsamechre,key=lambda x: x[-1],reverse=True)
	
		print allsamechre

		for allsamechre0 in allsamechre:
	
			pd.DataFrame.from_records([[twofiles[0]]+allsamechre0[0], [twofiles[1]]+allsamechre0[1]]).to_csv('allblastnpair.csv', header=None,sep=',',index=False, mode='a')


	for result in allsamechre:

		ref_chr=result[0][0]

		front_break_ref=result[0][3]
		front_break_query=cuts[0][1]-result[0][5]
	
		back_break_ref=result[1][2]
		back_break_query=cuts[1][0]+result[1][4]

		score=result[-1]

		strands=result[0][1]+result[1][1]
	
		out=infor+sizes+[front_break_query,back_break_query,ref_chr,strands,front_break_ref,back_break_ref,score]
	
		with open('allngaps_find.csv', mode='a') as f:
		
			f.write(','.join([str(x) for x in out])+'\n')
		
		f.close()
	
	
	return 0

	
def prepare(querypath,outputpath):
	
	def cutfasta(file0,tempfolder):
		
		zipped=0
		if '.gz' in file0:
			zipped=1
		
		if zipped==1:
			os.system('gunzip %s'%file0)
		
		filename=file0.split('/')[-1]
		
		with open(file0,mode='r') as f:
			
			reads=f.read().split('>')[1:]
		
		f.close()
		
		for i,read in enumerate(reads):
			
			outputfile0='{:s}/{:s}~{:d}'.format(tempfolder,filename,i)
			
			with open(outputfile0, mode='w') as f:
				f.write('>'+read)
			f.close()
		
		if zipped==1:
			os.system('gzip %s'%file0)
		
	
	tempfolder=outputpath+'/temp00'
	
	try:
		os.mkdir(outputpath+'/')
	except:
		pass	
	
	
	try:
		os.mkdir(tempfolder+'/')
	except:
		print "warnning, temp00 exits"
		
	if os.path.isfile(querypath):
	
		allfiles=[querypath]
		
	else:
		
		allfiles=[querypath+'/'+x for x in os.listdir(querypath)]
	
	for file0 in allfiles:
		
		cutfasta(file0,tempfolder)
	
	return tempfolder

def organizefiles(tempfolder,outputpath):
	
	allfiles=os.listdir(tempfolder)
	
	original_infor=cl.defaultdict(list)
	
	for file0 in allfiles:
		
		prefix=file0.split('+')[0]
		
		original_infor[prefix].append(file0)
		
	
	for eachfile in original_infor.keys():
		
		allcontigs=sorted(original_infor[eachfile], key=lambda x: x.split('+')[-1])
		
		for contigfile in allcontigs:
			
			os.system('cat {:s} >> {:s}'.format(tempfolder+'/'+contigfile, outputpath+'/'+eachfile))


class runminimap:
	
	def __init__(self,tempfolder,querypath,refpath):
		

		self.temp=tempfolder
		
		self.refpath=refpath

		self.query=querypath.split('/')[-1]
		
	
		
	def infor(self,i, j,title, cuts):	
		
		self.title=title

		self.i=i

		self.j=j
		
		self.cuts=cuts
		
		self.twofiles=['','']
		
		self.sizes=[abs(cuts[0][0]-cuts[0][1]), abs(cuts[1][0]-cuts[1][1])]
		
		#frontfile=file0+'~'+str(j)+'f'
		#backfile=file0+'~'+str(j)+'b'
	
		frontfile=self.temp+'/'+self.query+'_ngap.fa'
		backfile=self.temp+'/'+self.query+'_ngap.fa'

		self.path=['','']

		if self.sizes[0]>500:

			self.twofiles[0]='>%s_%d-%d'%(self.title,self.cuts[0][0],self.cuts[0][1])

			self.path[0]=self.temp+'/'+'%d_%d_f'%(self.i ,self.j)+'.fa'
			self.path[0]=frontfile	

		if self.sizes[1]>500:

			self.twofiles[1]='>%s_%d-%d'%(self.title,self.cuts[1][0],self.cuts[1][1])

			self.path[1]=self.temp+'/'+'%d_%d_b'%(self.i ,self.j)+'.fa'
			self.path[1]=backfile


		return self

	def writeanchor(self,query):


		frontfile=self.temp+'/'+self.query+'_ngap.fa'
		backfile=self.temp+'/'+self.query+'_ngap.fa'


		if self.sizes[0]>500:
			
			with open(self.path[0],mode='a') as f:
				f.write(self.twofiles[0]+'\n'+query[self.cuts[0][0]:self.cuts[0][1]]+'\n')
			f.close()
			
		if self.sizes[1]>500:
			
			with open(self.path[1],mode='a') as f:
				f.write(self.twofiles[1]+'\n'+query[self.cuts[1][0]:self.cuts[1][1]]+'\n')
			f.close()

		return self

	def findgap0(self,path=''):

		if path !='':
			self.path=path

		findgap(self.refpath,self.temp,self.path,self.twofiles,self.sizes,self.cuts,[self.query,self.title,self.j,self.cuts[0][1],self.cuts[1][0]])



def runngap(ngap,path):

	ngap.findgap0(path)	




def correctloop(tempfolder,querypath,refpath,allgaps=''):
	

	with open(querypath,mode='r') as f:
	
		allcontigs=f.read().split('>')[1:]

	f.close()



	allngaps=[]	
	minimap=runminimap(tempfolder,querypath,refpath)
	for i,x in enumerate(allcontigs):
		
		print('start NO.'+str(i)+'th turn\n\nstart calling minimap\n')
		
		print('start correction\n')
		
		query=x.splitlines()
	
		title=query[0].split(' ')[0].strip()
		if len(allgaps[title])==0:
			continue

		query=''.join(query[1:])
		
		if allgaps=='':
			allgaps0=zip(*distractgap(query))
		else:		
			allgaps0=allgaps[title]
		#print allgaps0

		anchorsize=5000	
		allcuts=[((max(0, gap[0]-anchorsize), gap[0]),(gap[1], min(len(query), gap[1]+anchorsize)  )) for gap in allgaps0 if abs(gap[1]-gap[0])>=100]
	
		for j,cuts in enumerate(allcuts):

			allngaps.append(minimap.infor(i,j,title, cuts).writeanchor(query))

	
	allngap_files=[tempfolder+'/'+x for x in os.listdir(tempfolder) if '.fa' in x and x+'_out' not in os.listdir(tempfolder)]


	for ngap_file in allngap_files:

		cmd='minimap2 -x ava-pb -splice -m 200 -N 10 -g 1000 -G 20 --dual=yes -t 32 {:s} {:s} > {:s} 2>&-'.format(ngap_file, refpath, ngap_file+'_out')
		print cmd
	
		os.system(cmd)
	
	#p=mul.Pool(processes=16)
	for ngap in allngaps:
		runngap(ngap,'')
		#p.apply_async(runminimap,(ngap,''))
	
	
	#p.close()
	#p.join()
		
				
		

		
	return 0



def run(querypath,refpath,outputpath,gapfile,tempfolder):

	print querypath, refpath, outputpath, gapfile

	#tempfolder=prepare(querypath,outputpath)
	
	t=pd.read_csv(gapfile, header=None, sep='\t')
	
	t=t[t[3]>=500]

	allgaps=cl.defaultdict(list)

	for title, start, end in zip(list(t[0]),list(t[1]),list(t[2])):
		allgaps[str(title)].append([start,end])
	

	correctloop(tempfolder,querypath,refpath,allgaps)

	#organizefiles(tempfolder,outputpath)

	try:
		os.system('rm -rf {:s}/{:s}'.format(outputpath,tempfolder))
	except:
		print('warning: unable to delect temporary folder, please do manually\n')

	print 'finished program\n'
	

def main(querypath,refpath, outputpath):
	
	allinputfiles=os.listdir(querypath)
	
	for i0,thefile in enumerate(allinputfiles[:1]):

		tempfolder=outputpath+'temp'+str(i0)

		try:
			os.mkdir(tempfolder)
		except:
			pass
		
		prefix=thefile.split('_')[0]
		
		path=allpaths[prefix]+'/'+thefile.split('ngap')[0]+'fasta'
	
		run(path,refpath,outputpath,querypath+thefile,tempfolder)

	
	


#main fuction to run step by step 

if __name__=='__main__':


	print 'starting'
	
	main(querypath,refpath,outputpath)
