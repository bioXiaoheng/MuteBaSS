#2019.3.18 Update: fixed the chisq module to be compatible with empty windows
#2018.3.17 Update: modified HKApara() to accomodate Ballette output
import numpy as np

#cfg = {(sample size): [allpoly, allsub]}
#n is the number of species
def HKApara(n, configfile):
    cfg = {}
    with open(configfile,'r') as config:
        for l in config:
            l = l.strip().split('\t')
            size = tuple([int(l[i]) for i in range(n)]) #first n for sample size
            poly = float(l[n]); sub = float(l[n+1])
            #poly = sum([float(l[i]) for i in range(n,2*n)]) #next n for poly
            #if n==1:
            #    sub = float(l[2])
            #else:
            #    sub = sum([float(l[i]) for i in range(2*n, 4*n-3)])
            cfg[size] = [poly,sub]
    return cfg

#return true for substitutions
def subst(l):
    SN = int(len(l[1:])/2)
    if SN > 1: #trans input
        f = [ l[2*i-1]/l[2*i] for i in range(1,SN+1)]
        #if abs(l[0]-l[1]) == 1:
        if set(f) == set([0,1]):
            return True
    #single species output
    elif float(l[1])/l[2] in [0,1]: 
        return True
    return False

#HKA calculation for n species, n=1,2,...,n
#window format: pos x1 n1 x2 n2 ...; x and n should already converted
def HKA(n, window, cfg):
    Ocounts = {}
    #default to ignore sample sizes no included in cfg
    for l in window:
        #obtain sample size
        size = tuple([l[2*i+2] for i in range(n)])
        if size not in cfg:
            print('line',' '.join([str(x) for x in l]),'have unrecognized sample size.')
            continue
        if size not in Ocounts:
            Ocounts[size] = [0,0,0] #poly, sub, # of sites
        elif subst(l):
            Ocounts[size][1] += 1
            Ocounts[size][2] += 1
        else:
            Ocounts[size][0] += 1
            Ocounts[size][2] += 1
    #obtain expected counts
    N=len(window)
    Xsq = 0
    for size in Ocounts:
        O = Ocounts[size]

        weight = float(O[2])/N
        #calc chisq
        E = [float(O[2])*frac for frac in cfg[size]]
        try:
            chisq = [(E[i]-O[i])**2 / E[i] for i in range(2)]
        except:
            #print E, O
            continue
        #add sign:
        if E[0] <= O[0]: #no less poly than expected
            Xsq += weight*sum(chisq)
        else:
            Xsq -= weight*sum(chisq)

    return(Xsq)

