import os
import sys
import string
from Ganga import *

j=Job()
j.application = Executable(exe=File('/gridgroup/ilc/garillot/SDHCALMarlinProcessor/script/EfficiencyProcessorOnGrid.py'), args=['arg1','arg2','arg3'])
j.backend='CREAM'
j.backend.CE='lyogrid07.in2p3.fr:8443/cream-pbs-calice'
j.comment = 'eff processor'

def frange(x, y, jump) :
    while y >= x - 1e-10 :
        if abs( round(y,0)-y ) < 1e-10 :
            yield int( round(y,0) )
        else :
            yield y
        y -= jump


#qbar = [ frange(0.5 , 10 , 0.5) ]
#delta = [ frange(0.5 , 6 , 0.5) ]
#d = [ 0.05 , 0.1 , 0.2 , 0.4 , 0.7 , 1 ]

qbar = [ 5 ]
delta = [ 3 ]
d = [ 0.2 ]

par=[]
args = [ [ str(q) , str(de) , str(dd) ] for q in qbar for de in delta for dd in d ]

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
