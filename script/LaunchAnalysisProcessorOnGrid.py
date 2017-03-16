import os
import sys
import string
from Ganga import *

j=Job()
j.application = Executable( exe=File('/gridgroup/ilc/garillot/SDHCALMarlinProcessor/script/AnalysisProcessorOnGrid.py'), args=['arg1','arg2','arg3','arg4'] )
j.backend = 'CREAM'
j.backend.CE = 'lyogrid07.in2p3.fr:8443/cream-pbs-calice'
j.comment = 'Analysis'

def frange(x, y, jump) :
    while y >= x - 1e-10 :
        if abs( round(y,0)-y ) < 1e-10 :
            yield int( round(y,0) )
        else :
            yield y
        y -= jump


#particle = ['pi-']
#energy = [90,80,70,60,50,40,30,20,10]
#model = ['FTFP_BERT_HP' , 'QGSP_BERT_HP' , 'FTF_BIC']
#version = ['9.6']

#par=[]
#args = [ [ p , str(e) , m , str(v) ] for p in particle for e in energy for m in model for v in version ]


par=[]
args=[]

version = ['9.6']

particle = ['pi-']
energy = [90,80,70,60,50,40,30,20,10]
model = ['FTFP_BERT_HP' , 'QGSP_BERT_HP' , 'FTF_BIC' , 'FTFP_BERT_EMY' , 'QGSP_BERT_EMY']

args.extend( [ [ p , str(e) , m , str(v) ] for p in particle for e in energy for m in model for v in version ] )

particle2 = ['e-']
energy2 = [85 , 50 , 40 , 30 , 25 , 20 , 15 , 10]
model2 = ['FTFP_BERT' , 'QGSP_BERT' , 'FTFP_BERT_EMY' , 'QGSP_BERT_EMY']

args.extend( [ [ p , str(e) , m , str(v) ] for p in particle2 for e in energy2 for m in model2 for v in version ] )

for i in args:
	par.append(i)


par
s = ArgSplitter()
s.args=par
j.splitter=s
j.submit()

for subj in j.subjobs :
	subjComment = ""
	for arg in subj.application.args :
		subjComment = subjComment + " " + str(arg)
	subj.comment = str(subjComment)
