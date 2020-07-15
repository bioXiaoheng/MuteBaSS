#Feb 22, 2018: add in NCD1b, NCD1opt, and NCDmid; optimized NCDopt
#Oct 13, 2017: major revisions: add in NCDs
#Oct 4, 2017: Do folded only, and modify spec to [0.3, 0.5]
#script of new NCDopt: evaluate every tf and pick the lowest, and NCDb: take the central SNP frequency as tf
import numpy as np
import sys
spec = [.01*i for i in range(30,51)]
s_spec = [.01*i for i in range(1,51)]

#return true for substitutions
def subst(l):
    SN = int(len(l[1:])/2)
    if SN > 1: #trans input
        f = [ l[2*i-1]/l[2*i] for i in range(1,SN+1)]
        #if abs(l[0]-l[1]) == 1:
        if set(f) == set([0,1]):
            return True
    elif float(l[1])/l[2] in [0,1]: 
        return True
    return False


#Basic unit to calc NCD for an array of frequencies.
def NCD(f, tf):
	N = len(f)
	f = np.array(f)
	vmin = np.vectorize(min)
	f = vmin(1.-f, f) #minor allele frequency
	#print f
	tf = min(1.-tf,tf)
	var = (f-tf)**2
	return((np.sum(var)/N)**.5)

#original NCD2
#when windowB is given instead of window, it automatically becomes NCD1
def NCD2(window,tf):
	if len(window) == 0:
		return('NA')
	SN = int(len(window[0])/2) #number of species
	#for more than one species
	F = [[l[i*2-1]/l[i*2] for i in range(1,SN+1)] for l in window]
	#convert to mAF
	MAF = [[min(x, 1.-x) for x in l] for l in F]
	f=[max(l) for l in MAF]
	return(NCD(f,tf))

#all numbers in window should have been converted to float
#NCDopt (aka NCDopt) evaluate every tf and pick the lowest
def NCDopt(window):
	SN = int(len(window[0][1:])/2)#number of species
	#for more than one species
	F = [[l[i*2-1]/l[i*2] for i in range(1,SN+1)] for l in window]
	#convert to mAF
	MAF = [[min(x, 1.-x) for x in l] for l in F]
	f=[max(l) for l in MAF]
	#for l in MAF:
	#	l = [x for x in l if x != 0]
	#	f.append(max(l))
	NCDopt = [1,0.5]
	for tf in spec:
		ncd = NCD(f,tf)
		if ncd <= NCDopt[0]:
			NCDopt = [ncd,tf]
	return(NCDopt)

#take subs and polys separately
def NCDs(window):
	SN = int(len(window[0][1:])/2) #number of species
	#convert to frequency
	F = [[l[i*2-1]/l[i*2] for i in range(1,SN+1)] for l in window]
	#convert to mAF
	MAF = [[min(x, 1.-x) for x in l] for l in F]
	f=[max(l) for l in MAF]
	#separate poly and subs
	fpoly = [x for x in f if x != 0]
	fsubs = [x for x in f if x == 0]
	N = len(f) #total # of sites
	fpoly = np.array(fpoly)
	#print f
	NCDS = 1; f0=0
	for tf in s_spec:
		var = (fpoly-tf)**2
		ncd = ((np.sum(var) + len(fsubs))/(4*N))**.5
		if ncd <= NCDS:
			NCDS=ncd; f0=tf
	return(NCDS,f0)

#NCDb take central SNP freqency as tf
def NCDb(window, f0):
	#print window
	N = len(window)
	SN = int(len(window[0][1:])/2) #number of species
	if N <= 1:
		return 1
	#if there are more than 1 site
	F = [[l[i*2-1]/l[i*2] for i in range(1,SN+1)] for l in window]
	ncd = NCD(F,f0)
	return(ncd)

#NCDbs calculate NCD1b+NCDsub(0.5)
def NCDmid(window, f0):
	windowB = [l for l in window if not subst(l)]
	SN = int(len(window[0][1:])/2) #number of species
	#convert to frequency
	F = [[l[i*2-1]/l[i*2] for i in range(1,SN+1)] for l in window]
	#convert to mAF
	MAF = [[min(x, 1.-x) for x in l] for l in F]
	f=[max(l) for l in MAF]
	#separate poly and subs
	fpoly = [x for x in f if x != 0]
	fsubs = [x for x in f if x == 0]
	N = len(f) #total # of sites
	fpoly = np.array(fpoly)
	var = (fpoly-f0)**2
	#calculate weighted NCD1b+NCDsub(0.5)
	SB = float(len(fsubs))/(4*N) + len(fpoly)/N * np.sum(var)
	return(SB**.5)





