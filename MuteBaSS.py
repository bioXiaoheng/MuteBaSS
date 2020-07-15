#python script to perform sliding windows for summary stats calculation
#2018.3.17 Update: 1. Modified HKAconfig so it's compatible with new Ballette's output
                ##2. add in -c commands for the column indice where pos & x1 start
#2018.3.16 Update: 1. incorporated more NCDs (as of 2018.2)
            #2. fixed a bug in slideWin() scanning pipeline (as of 3.16)
#2018.1.29 Update: 1. re-write the fixed-size SNP-centered window.
                #2. substitute "while l!=['']" to "for l in file"
import sys, os, optparse, re
import HKAtrans, transNCDs

#return true for substitutions
def subst(l):
    SN = int(len(l[1:])/2)
    if SN > 1: #trans input
        f = [ float(l[2*i-1])/l[2*i] for i in range(1,SN+1)]
        #if abs(l[0]-l[1]) == 1:
        if set(f) == set([0,1]):
            return True
    #single species output
    elif float(l[1])/l[2] in [0,1]: 
        return True
    return False

#step is the number of snps between test sites
#r for radius, the number of snps flanking the test locus on each side
#one window contains 2*r+1 lines
def slideSNP(infile, pi, xi, outfile, step, r,configfile,**keyword):
    methods = []
    if 'HKA' in keyword:
        if keyword['HKA'] == True:
            methods.append('HKA')
    if 'NCD' in keyword:
        if keyword['NCD']==True:
            methods.append('NCD')
            tf =keyword['tf']
    if 'NCDopt' in keyword:
        if keyword['NCDopt']==True:
            methods.append('NCDopt')
    if 'NCDsub' in keyword:
        if keyword['NCDsub']==True:
            methods.append('NCDsub')

    if 'NCDopt' in methods or 'NCDsub' in methods:
        header = 'sitePos\t'+'\t'.join(methods)+'\toptF\n'
    else:
        header = 'sitePos\t'+'\t'.join(methods)+'\n'

    score = open(outfile,'w')
    score.write(header)
    snp = open(infile,'r')
    l=snp.readline().strip().split('\t') #header
    SN=int(len(l[xi:])/2)
    #print('Finding balancing selection affecting', SN , 'species.')
    #obtaining HKA parameters
    if 'HKA' in methods:
        config = HKAtrans.HKApara(SN, configfile)
    print('Start scanning %s-species input with %s.\nWindow centered on every %s sites, covering %s sites on either side.' % (SN,', '.join(methods), step,r))

    window=[]
    for l in snp:
        l = l.strip().split('\t')
        l = [l[pi]]+[float(x) for x in l[xi:]]
        #first window:
        if len(window) < 2*r+1:
            window.append(l)
            continue
        else: #calculate & write down the scores
            snpPos = window[r][0]
            Score={}; f0=[]
            if 'HKA' in methods:
                try:
                    Score['HKA'] = HKAtrans.HKA(SN, window,config)
                except:
                    print(config)
            if 'NCD' in methods:
                Score['NCD'] = transNCDs.NCD2(window, tf)
            if 'NCDopt' in methods:
                ncda= transNCDs.NCDopt(window)
                Score['NCDopt'] = ncda[0]; f0.append(ncda[1])
            if 'NCDsub' in methods:
                ncds= transNCDs.NCDs(window)
                Score['NCDsub'] = ncds[0]; f0.append(ncds[1])
            if 'NCDopt' in methods or 'NCDsub' in methods:
                ol='%s\t'%(snpPos)+'\t'.join([str(Score[m]) for m in methods])+'\t%s\n'%( ','.join([str(f) for f in f0]))
            else:
                ol='%s\t'%(snpPos)+'\t'.join([str(Score[m]) for m in methods])+'\n'
            score.write(ol)
            #take the next step
            window = window[int(step):]
    if len(window) > 2*r+1-step:
        snpPos = window[r][0]
        Score={}; f0=[]
        if 'HKA' in methods:
            try:
                Score['HKA'] = HKAtrans.HKA(SN,window,config)
            except:
                print(config)
        if 'NCD' in methods:
            Score['NCD'] = transNCDs.NCD2(window, tf)
        if 'NCDopt' in methods:
            ncda= transNCDs.NCDopt(window)
            Score['NCDopt'] = ncda[0]; f0.append(ncda[1])

        if 'NCDsub' in methods:
            ncds= transNCDs.NCDs(window)
            Score['NCDsub'] = ncds[0]; f0.append(ncds[1])

        if 'NCDopt' in methods or 'NCDsub' in methods:
            ol='%s\t'%(snpPos)+'\t'.join([str(Score[m]) for m in methods])+'\t%s\n'%( ','.join([str(f) for f in f0]))
        else:
            ol='%s\t'%(snpPos)+'\t'.join([str(Score[m]) for m in methods])+'\n'
        score.write(ol)

    score.close()
    #print 'Pipeline finished'

#fixed size
def slideWin(infile,pi, xi, outfile, step, w, configfile, **keyword):
    methods = []
    if 'HKA' in keyword:
        if keyword['HKA'] == True:
            methods.append('HKA')
    if 'NCD' in keyword:
        if keyword['NCD']==True:
            methods.append('NCD')
            tf =keyword['tf']
    if 'NCDopt' in keyword:
        if keyword['NCDopt']==True:
            methods.append('NCDopt')
    if 'NCDsub' in keyword:
        if keyword['NCDsub']==True:
            methods.append('NCDsub')

    if 'NCDopt' in methods or 'NCDsub' in methods:
        header = 'midpos\t'+'\t'.join(methods)+'\tnumSites\toptF\n'
    else:
        header = 'midpos\t'+'\t'.join(methods)+'\tnumSites\n'
    score = open(outfile,'w')
    score.write(header)
    snp = open(infile, 'r')
    l=snp.readline().strip().split('\t') #header
    SN = int(len(l[xi:])/2)
    #print('Finding balancing selection affecting', SN , 'species.')
    if 'HKA' in methods:
        config = HKAtrans.HKApara(SN,configfile)

    print('Start scanning %s-species input with %s.\nWindow size %sbp, step size %sbp.' % (SN,', '.join(methods),w, step))
    window = []
    start = 0; end = start + w; midpos = start + w/2
    l=snp.readline().strip().split('\t')
    while l != ['']:
        if float(l[pi]) < end:
            l = [l[pi]]+[float(x) for x in l[xi:]]
            window.append(l)
            l=snp.readline().strip().split('\t')
        elif float(l[pi]) >= end:
            if len(window)>0:
                Score={}; f0=[]
                if 'HKA' in methods:
                    Score['HKA'] = HKAtrans.HKA(SN,window,config)
                if 'NCD' in methods:
                    Score['NCD'] = transNCDs.NCD2(window, tf)
                if 'NCDopt' in methods:
                    ncda= transNCDs.NCDopt(window)
                    Score['NCDopt'] = ncda[0]; f0.append(ncda[1])
                if 'NCDsub' in methods:
                    ncds= transNCDs.NCDs(window)
                    Score['NCDsub'] = ncds[0]; f0.append(ncds[1])

                if 'NCDopt' in methods or 'NCDsub' in methods:
                    ol='%s\t'%(int(midpos))+'\t'.join([str(Score[m]) for m in methods])+'\t%s\t%s\n'%(len(window), ','.join([str(f) for f in f0]))
                else:
                    ol='%s\t'%(int(midpos))+'\t'.join([str(Score[m]) for m in methods])+'\t%s\n'%(len(window))
                score.write(ol)
            #take a new step
            start += step; end += step; midpos += step
            window = [l  for l in window if float(l[0])>=start]
            l=snp.readline().strip().split('\t')

    if len(window)>1:
        Score={}; f0=[]
        if 'HKA' in methods:
            Score['HKA'] = HKAtrans.HKA(SN,window,config)
        if 'NCD' in methods:
            Score['NCD'] = transNCDs.NCD2(window, tf)
        if 'NCDopt' in methods:
            ncda= transNCDs.NCDopt(window)
            Score['NCDopt'] = ncda[0]; f0.append(ncda[1])
        if 'NCDsub' in methods:
            ncds= transNCDs.NCDs(window)
            Score['NCDsub'] = ncds[0]; f0.append(ncds[1])

        if 'NCDopt' in methods or 'NCDsub' in methods or 'NCD1a' in methods:
            ol='%s\t'%(int(midpos))+'\t'.join([str(Score[m]) for m in methods])+'\t%s\t%s\n'%(len(window), ','.join([str(f) for f in f0]))
        else:
            ol='%s\t'%(int(midpos))+'\t'.join([str(Score[m]) for m in methods])+'\t%s\n'%(len(window))
        score.write(ol)
    score.close()

#assuming input is poly+sub mixed. centered on every polymorphic site
#r is the number of sites flanking the centered SNP
def slideSNPct(infile, pi, xi, outfile,r,**keyword):
    methods = []
    if 'NCD' in keyword:
        if keyword['NCD']==True:
            methods.append('NCD')
            tf =keyword['tf']
    if 'NCDmid' in keyword:
        if keyword['NCDmid']==True:
            methods.append('NCDmid')
    header = 'snpPos\t'+'\t'.join(methods)+'\tnumSNPs\tfc\n'
    score = open(outfile[:-4]+'_snpCT.txt', 'w')
    score.write(header)
    snp = open(infile,'r')
    l = snp.readline() #header
    windowAll = []
    l = snp.readline().strip().split('\t') #first line to get ready for while loop
    SN = int(len(l[xi:])/2)
    print('Starting scanning %s-species input with %s.\nWindow centered on every SNP, containing %s sites on either side.' % (SN, ' '.join(methods), r))
    #print header, methods
    #read the first window
    for i in range(r*2):
        l = [l[pi]]+[float(x) for x in l[xi:]]
        windowAll.append(l)
        l = snp.readline().strip().split('\t')
        if l == ['']: 
            #print 'e.o.f.'
            break
    while l != ['']:
        l = [l[pi]]+[float(x) for x in l[xi:]]
        windowAll.append(l)
        if len(windowAll) < 2*r+1:
            continue
        windowAll = windowAll[-(2*r+1):]
        #if mid-snp is substitution, keep reading till it reaches a poly
        while subst(windowAll[r]): 
            l = snp.readline().strip().split('\t')
            if l == ['']: 
                break
            l = [l[pi]]+[float(x) for x in l[xi:]]
            windowAll.append(l)
            windowAll = windowAll[-(2*r+1):]
        f0 = windowAll[r]; del windowAll[r]
        snpPos = f0[0]
        f0 = [float(f0[i*2-1])/float(f0[i*2]) for i in range(1,SN+1)]
        #remove subs
        f0 = [x for x in f0 if x not in [0,1]]
        try:
            f0=max(f0)
        except: #when l=['']
            f0='NA'
            continue
        #record scores
        Scores = {}
        if 'NCD' in methods:
            Scores['NCD'] = transNCDs.NCD2(windowAll, tf)
        if 'NCDmid' in methods:
            Scores['NCDmid'] = transNCDs.NCDmid(windowAll, f0)
        windowB = [l for l in windowAll if not subst(l)]
        score.write('%s\t%s\t%s\t%s\n'%(snpPos, '\t'.join([str(Scores[m]) for m in methods]), len(windowB), f0))
        #take a new step
        l = snp.readline().strip().split('\t')
    score.close()
    snp.close()


#assuming input is poly+sub mixed. centered on every polymorphic site
#fixed window size w, with w/2 to either side of center
def readMore(filevar, pi, xi, All, SNPs, Pos, numL):
    eof = False
    for i in range(numL):
        l = filevar.readline().strip().split('\t')
        if l == ['']:
            eof = True
            #print 'e.o.f'
            break
        l = [float(l[pi])]+[float(x) for x in l[xi:]]
        All.append(l) ; Pos.append(l[0])
        if not subst(l):
            SNPs.append(l)
    return(filevar, All, SNPs, Pos, eof)

#all folded?
def slideSizeCT(infile,pi, xi, outfile, w, **keyword):
    methods = []
    if 'NCD' in keyword:
        if keyword['NCD']==True:
            methods.append('NCD')
            tf =keyword['tf']
    if 'NCDmid' in keyword:
        if keyword['NCDmid']==True:
            methods.append('NCDmid')

    header = 'snpPos\t'+'\t'.join(methods)+'\tnumSites\tnumSNPs\tfc\n'#tstartSNP\tendSNP\tstartSite\tendSite\
    score = open(outfile[:-4]+'_fixSizeCT.txt', 'w')
    score.write(header)
    #read input file
    sites = open(infile,'r'); eof = False
    l = sites.readline().strip().split(); SN = int(len(l[xi:])/2) #header
    print('Start scanning %s-species input with %s.\nWindow centered on each SNP, spanning %sbp.' % (SN,', '.join(methods),w))
    #Initiating
    Stored = []; SNPs = []; Pos =[]
    #read first 1000 lines, check if larger than window size
    while eof == False:
        (sites, Stored, SNPs, Pos, eof) = readMore(sites,pi,xi, Stored, SNPs, Pos, 1000)
        #read more if stored region is smaller than window size
        if Stored[-1][0] - Stored[0][0] < w:
           continue
        #Grab sites given each snp center and output results
       #print SNPs[:5]
        i=0; start_i = 0
        for snp in SNPs:
            i+=1
            start = max(snp[0]-w/2, 0); end = snp[0]+w/2
            #read more if last window not fully covered
            if end > Stored[-1][0]:
                break
            #initiating
            try:
                f0 = [snp[2*j-1]/snp[2*j] for j in range(1,SN+1)]
            except:
                print(j, range(1,SN+1), snp)
                sys.exit()
            try:
                f0 = [x for x in f0 if x not in [0.,1.]][0]
            except:
                print(f0, snp) #something with the parsing
                continue
            Window = []
            #grab sites within the window
            for j in range(start_i, len(Pos)):
                if Pos[j] < start:
                    start_i += 1; continue
                elif Pos[j] >= start and Pos[j] <= end:
                    Window.append(Stored[j])
                else: #pos > end
                    break
            #calculate stuff
            Scores = {}
            if 'NCD' in methods:
                Scores['NCD'] = transNCDs.NCD2(Window, tf)
            if 'NCDmid' in methods:
                Scores['NCDmid'] = transNCDs.NCDmid(Window, f0)

            windowB = [l for l in Window if not subst(l)]
            score.write('%s\t%s\t%s\t%s\t%s\n'%(int(snp[0]), '\t'.join([str(Scores[m]) for m in methods]), len(Window), len(windowB), f0))#,windowB[0][0], windowB[-1][0],windowAll[0][0], windowAll[-1][0]\t%s\t%s\t%s\t%s
        #take a new step, start from (i-1)th snp in SNP
        SNPs = SNPs[(i-1):]
        Stored = [x for x in Stored if x[0] >= SNPs[0][0]-w/2]
        Pos = [x[0] for x in Stored]

    #print len(windowAll)
    #print windowAll
    score.close()
    sites.close()

def tsp(l):#l: pos x1 n1 x2 n2 etc
    SN = int((len(l)-1)/2)
    f = [l[2*i+1]/l[2*i+2] for i in range(SN)]
    poly = [x>0 and x<1 for x in f]
    if sum(poly) > 1: 
    #multi poly
        return 1
    elif sum(poly) == 1 and len(set([f[i] for i in range(SN) if not poly[i]])) >= 2:
    #poly+0+1
        return 2
    else:
    #within-sp poly
        return 0

#convert e.g. (((1,2),3),4) to the list of compatible patterns
def getTrees(Tree):
    #get the # of species (i.e. # of commas+1)
    SN = max([int(s) for s in re.findall(r'\d+', Tree)])
    trees = []
    for i in range(SN):
        t = SN*[0]
        t[i] = 1
        trees.append(tuple(t))
        t = [1-x for x in t]
        trees.append(tuple(t))
    #parse the string
    #left = [i for i in range(len(Tree)) if Tree[i]=='(']
    left = [ i.start() for i in re.finditer('\(',Tree)]
    #right = [i for i in range(len(Tree)) if Tree[i]==')']
    right = [ i.start() for i in re.finditer('\)',Tree) ]
    assert len(left)==len(right)
    pairs=[]
    for n1 in left[::-1]:
        i=0
        while right[i] <= n1:
            i+=1
        pairs.append([n1,right[i]])
        right.remove(right[i])
    #generate patterns
    for p1,p2 in pairs:
        #extract the numbers
        seg = Tree[p1+1:p2]
        sp = [int(s) for s in re.findall(r'\d+', seg)]
        t = SN*[0]
        for p in sp:
            t[p-1] = 1
        trees.append(tuple(t))
        t = [1-x for x in t]
        trees.append(tuple(t))
    trees = list(set(trees)) #remove re-occurrence
    if SN*[0] in trees:
        trees.remove(SN*[0])
    if SN*[1] in trees:
        trees.remove(SN*[1])
    return trees

#use when given more than 3 species
def checkTree(l,trees):
    SN = int((len(l)-1)/2)
    if SN <= 3:
        return True
    else:
        sub = [l[2*i+1]/l[2*i+2] for i in range(SN)]
        if tuple(sub) in trees:
            return True
        else:
            return False


def check(infile, indice, Tree):
    if indice.count(',')!= 1:
        print('Error. Please indicate column indices in the correct format: p,x')
        sys.exit()
    else:
        indice=indice.split(',')
        pi=int(indice[0])-1; xi=int(indice[1])-1

    #Checking input file
    Nline = 1; pos=0
    with open(infile) as f:
        l = next(f).strip().split('\t') #header (or not)
        colnum = len(l)
        if colnum < 3:
            print('Error: input format. Please make sure your file is tab-delimited, and includes at least 3 columns.')
            sys.exit()
        if len(l[xi:])%2 != 0:
            print('Error: please include observed and total counts (two columns) for every species.')
            sys.exit()
        SpNum=int(len(l[xi:])/2)
        print('Input file includes %s species.'%(SpNum))
        if SpNum > 3 and Tree == 1:
            print('Error: for more than 3 species, please include the species tree with \"--tree\" command, e.g. (((1,2),3),4)')
            sys.exit()
        elif SpNum > 3:
            print('Tree topology: %s'%(Tree))
            trees = getTrees(Tree)
            if trees == False:
                print('Error: please make sure your tree matches the number of species included in your input file.')
                sys.exit()
        #read file
        for l in f:
            Nline += 1
            lp = l
            l = l.strip().split('\t')
            if len(l) != colnum:
                print('Error: please make sure every row has the same number of columns. The file should be tab-delimited.')
                sys.exit()
            #check the position
            try: 
                l[pi]=float(l[pi])
            except(ValueError):
                print('line %s: %s' % (Nline, lp))
                print('Error: make sure positions are represented in numbers.')
                sys.exit()
            if float(l[pi]) <= pos and pos > 0:
                print('Error: make sure the positions are listed in ascending order.')
                sys.exit()
            pos = float(l[pi])
            #Extract the line
            l = [l[pi]]+l[xi:]
            try:
                l = [float(x) for x in l]
            except(ValueError):
                print('line %s: %s' % (Nline, lp))
                print('Error: allele counts can only be numbers.')
                sys.exit()
            #for multiple species, check sub & poly
            if SpNum>1:
                freq = [l[i*2+1]/l[i*2+2] for i in range(SpNum)]
                if subst(l):
                    if len(set(freq))==1:
                        print('line %s: %s' % (Nline, lp))
                        print('Error: please do not include monomorphical sites.')
                        sys.exit()
                    elif not checkTree(l,trees):
                        print('line %s: %s' % (Nline, lp))
                        print('Error: please only include substitutions that are compatible with the tree topology.')
                        sys.exit()
                else:
                    if tsp(l)==1:
                        print('line %s: %s' % (Nline, lp))
                        print('Error: please do not include trans-species polymorphic sites.')
                        sys.exit()
                    elif tsp(l)==2:
                        print('line %s: %s' % (Nline, lp))
                        print('Error: do not include within-species polymorphic sites that are monomorphically different in other species.')
                        sys.exit()
            elif SpNum==1 and l[1]==l[2]: #subs in single species
                print('line %s: %s' % (Nline, lp))
                print('Error: for single-species input, please provide ancestral allele frequency, and do not include monomorphically ancestral sites.')
                sys.exit()
    print('Error check finished. File format correct')

def getConfig(inputfile, pi, xi, configfile):
    #assuming the format is correct
    Counts = {}; #total=0
    with open(inputfile,'r') as sites:
        for l in sites:
            #lp=l
            l=l.strip().split('\t')
            try:
                l = [l[pi]]+[int(s) for s in l[xi:]]
                #total+=1
            except(ValueError):
                continue #skip header
            SN = int((len(l)-1)/2)
            size = tuple([l[2*i+2] for i in range(SN)])
            if size not in Counts:
                Counts[size]=[0,0] #poly, sub

            if subst(l):
                Counts[size][1] += 1
            else:
                Counts[size][0] += 1
    #print Counts
    with open(configfile,'w') as config:
        for size in Counts:
            total = sum(Counts[size])
            frac = [float(Counts[size][i])/total for i in [0,1]]
            ol = '%s\t%s\t%s\n' % ( '\t'.join([str(n) for n in list(size)]) , frac[0], frac[1])
            config.write(ol)
    print('Finished generating configuration file.')


def main():
    #parsing arguments
    parser = optparse.OptionParser(usage='python %prog -i <input file> -c <p,x> [--check] [--tree <tree>] [--fixSize] -w <window size> -s <step size> [--HKA][--config <config file>] [--getConfig] [--NCD][--tf <tf>] [--NCDopt] [--NCDsub] [--NCDmid] -o <output file>')
    parser.add_option('-i','--input', dest='infile', help = 'Path and name of your input file.\n')
    parser.add_option('-o','--output', dest='outfile', help = 'Path and name of your output file.\n')
    parser.add_option('-c','--indices', dest='index', help = 'Index numbers for the columns of locus positions p, and allele counts of the first species x, respectively. Format of this argument should be \"p,x\", without space.\n')

    parser.add_option('--fixSize', action='store_true', dest = 'size', default = False, help = 'Option to fix the size of scanning windows. When true, provide the length of window in bp with \"-w\" command.\n')

    parser.add_option('-w','--window', dest='r', type = 'int', help='Number of sites flanking the test locus on either side. When choose to fix window size (\"--fixSize\"), input the length of window in bp.\n')
    parser.add_option('-s','--step', dest='step', type = 'float', help='For windows with fixed number of sites, each step is the number of sites between each neighboring test site. For windows with fixed length, each step is the length, in bp, between neighboring test sites.\n')
    
    
    parser.add_option('--HKA', action = 'store_true', dest='HKA', default=False, help = 'Option to perform HKA scan. User must provide (with \"--config\") configuration file, which records the fractions of substitutions and polymorphisms among all informative sites.\n')
    parser.add_option('--config', dest = 'configfile', help = 'Path and name of your configuration file for HKA scan. File should be tab-delimited, with the first k columns showing sample sizes for each species, then two columns for the corresponding fractions of within-species polymorphisms and substitutions, respectively.\n')
    parser.add_option('--getConfig', action = 'store_true', dest = 'getConfig', help = 'Option to generate configuration file from concatenated whole-genome (neutral) input file.\n')
    
    parser.add_option('--NCD', action = 'store_true', default=False, help = 'Option to perform NCD scan. Default target frequency (tf) is 0.5. User can customize the value with \"--tf\".\n')
    parser.add_option('--tf', dest = 'tf', default = 0.5, help = 'Target frequency for NCD scan. Default for tf is 0.5. tf must be no greater than 0.5.\n')
    
    parser.add_option('--NCDopt', dest = 'NCDopt', default = False, action = 'store_true', help = 'Option to perform NCDopt scan.')
    parser.add_option('--NCDsub', dest = 'NCDsub', default = False, action = 'store_true', help = 'Option to perform NCDsub scan.')
    parser.add_option('--NCDmid', dest = 'NCDmid', default = False, action = 'store_true', help = 'Option to perform NCDmid scan.')
        
    parser.add_option('--check', action='store_true', default=False, help='Option to check the format of input file.\n')
    parser.add_option('--tree', dest = 'tree', default = 1, help = 'Tree topology if given more than 3 species. Only needed when checking input format.\n')

    opt,v = parser.parse_args(sys.argv[1:])

    #print 'opt:',opt
    #print 'v',v
    #print sys.argv[1:]
    #print opt.p, opt.tf
    if sys.argv[1:] == []:
        parser.print_help()
        sys.exit()

    print('Scanning for balancing selection on %s. Output written to %s.' %(opt.infile, opt.outfile))

    if opt.check:
        print('Checking format of input file', opt.infile)
        check(opt.infile,opt.index,opt.tree)
        sys.exit()

    #Assume format is currect from now on, obtain indices
    indice = opt.index.split(',')
    pi = int(indice[0])-1 ; xi = int(indice[1])-1
    print('Allele position/name on column #%s, allele counts start at column #%s.' % (pi+1, xi+1))

    if opt.getConfig:
        print('Generating configuration file from %s...'%(opt.infile))
        getConfig(opt.infile,pi,xi,opt.configfile)

    #print opt
    if opt.NCDmid:
        print('Performing NCDmid scan...')
        if opt.size:
            print('with fixed-sized window centered on SNPs')
            slideSizeCT(opt.infile, pi, xi, opt.outfile, opt.r, NCD = opt.NCD, tf = float(opt.tf), NCDmid = opt.NCDmid)
        else:
            print('with fixed number of sites flanking center SNPs')
            slideSNPct(opt.infile, pi, xi, opt.outfile, opt.r, NCD = opt.NCD, NCDmid = opt.NCDmid, tf = float(opt.tf))

    if opt.NCD or opt.NCDopt or opt.HKA or opt.NCDsub:
        methods = ['HKA','NCD','NCDopt','NCDsub']
        yes = [opt.HKA, opt.NCD, opt.NCDopt, opt.NCDsub]
        methods = [methods[i] for i in range(4) if yes[i]]
        print('Performing %s scan/s...'%(', '.join(methods)))
        if not opt.size:
            slideSNP(opt.infile, pi, xi, opt.outfile, opt.step, opt.r, configfile = opt.configfile, NCD = opt.NCD, NCDopt = opt.NCDopt, tf = float(opt.tf), HKA = opt.HKA, NCDsub = opt.NCDsub)
        else:
            slideWin(opt.infile, pi, xi, opt.outfile, opt.step, opt.r, configfile = opt.configfile, NCD = opt.NCD, NCDopt = opt.NCDopt, tf = float(opt.tf), HKA = opt.HKA, NCDsub = opt.NCDsub)



if __name__ == '__main__':
    main()