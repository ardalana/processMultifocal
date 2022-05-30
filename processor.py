## Takes intersected mutations from all tumors in this group (patient), removes positions mutating nonspecifically in different tumors of the group, finds denovo mutations (non-dbsnp) in the group, finds non-repeat-region denovo mutations in the group and prints them for each tumor in inclusive and exclusive lists, also separately for indels and snps.
import sys,glob,pandas

rawCalls=glob.glob('*.'+sys.argv[1])
dbsnp=pandas.read_csv(sys.argv[2],sep='\t',header=None,names=['chrom','pos'],dtype={0:'object',1:'object'})
repeats=sys.argv[3]
## Make together all present mutations from all tumors in this patient, call as "total mutations"
allMutsIndices=[]
# allMutsLines=[]
for file in rawCalls:
    with open(file,'r') as fin:
        for line in fin:
            if not line.startswith('#') and (line.split('\t')[0],line.split('\t')[1],line.split('\t')[3],line.split('\t')[4]) not in allMutsIndices:
                allMutsIndices.append((line.split('\t')[0],line.split('\t')[1],line.split('\t')[3],line.split('\t')[4]))
                # allMutsLines.append((line.split('\t')[0],line.split('\t')[1],line.split('\t')[2],line.split('\t')[3],line.split('\t')[4],line.split('\t')[5],line.split('\t')[6],line.split('\t')[7],line.split('\t')[8],line.split('\t')[9],line.split('\t')[10].strip()))
print('groups total mutations: '+str(len(allMutsIndices)))
## Make a dataframe out of this list of tuples and then completely remove repeated positions, i.e positions that have mutated in more than one way in different tumors (and thus will have the same chr and pos values in allMuts tuples). Call the remaining allMutsClean or "monomorphic mutations". Though there were no non-uniquely evolved positions though in this case, the code was shown to work on a test set.
allMuts=pandas.DataFrame(allMutsIndices)
allMutsClean=allMuts.drop_duplicates(subset=(0,1),keep=False)
print('groups monomorphic mutations: '+str(len(allMutsClean)))

allMutsClean.columns=['chrom','pos','ref','alt']
# allMutsClean.columns=['chrom','pos','id','ref','alt','qual','filter','info','format','tumor','normal']
# print(allMutsClean.dtypes) ## shows all these columns to be 'object', so should modify dtype in dbsnp input (read_csv above) as 'object' too, so these are consistant and it becomes possible to do 'merge' for subtracting allMutsClean from dbsnp.

allMutsCleanList=[]
## Method 1 for converting a pandas dataframe into a list of tuples. Method 2 (below) did not work here for unkown reason. Can do both in one Function(def) in future improvements. Also need to add getting raw calls directly here (currently from R), and allMuts in pandas (currently from base python)
for i in range((allMutsClean.shape[0])):
    allMutsCleanList.append(tuple(allMutsClean.iloc[i,:]))
#print(len(allMutsCleanList))

## get those variants that are present in dbsnp, call them dbsnpVars
dbsnpVars=pandas.merge(allMutsClean,dbsnp,on=['chrom','pos'],how='inner')
print('input dbsnp sites: '+str(len(dbsnp)))
print('groups variants on dbsnps: '+str(len(dbsnpVars)))

dbsnpVarsList=[]
## Method 2 for converting a pandas dataframe into a list of tuples: it did not work with allMutsClean for unknown reason. Did not try for it itertuples(), exactly the same way as below, which may work.
for index, row in dbsnpVars.iterrows():
    rowTuple=(row.chrom,row.pos)
    dbsnpVarsList.append(rowTuple)

denovoMutsList=[]
for i in allMutsCleanList:
    if (i[0],i[1]) not in dbsnpVarsList:
        denovoMutsList.append(i)
print('groups denovo mutations: '+str(len(denovoMutsList)))
## Strangely did not work directly, so had to write denovoMutsList to file and then open the file again and make a list of it.
with open('denovoMutations','w') as fout:
    for i in sorted(denovoMutsList):
        print(i[0],i[1],i[2],i[3],sep='\t',file=fout)

denovoMuts=[]
with open('denovoMutations','r') as fin:
    for line in fin:
        denovoMuts.append((line.split('\t')[0],line.split('\t')[1],line.split('\t')[2],line.split('\t')[3].strip()))

denovoMutsOnRepeats=[]
repsNumLines=0
with open(repeats,'r') as fin:
    for line in fin:
        repsNumLines+=1
        for i in denovoMuts:
            if line.split('\t')[0]==i[0] and int(line.split('\t')[1])<=int(i[1])<=int(line.split('\t')[2].strip()):
                denovoMutsOnRepeats.append(i)
print('input repeat sites: '+str(repsNumLines))
print('denovo mutations on repeats: '+str(len(denovoMutsOnRepeats)))

denovoMutsOffRepeats=list(set(denovoMuts)-set(denovoMutsOnRepeats))
print('denovo mutations off repeats: '+str(len(denovoMutsOffRepeats)))

with open('denovoMutationsOffRepeats','w') as fout:
    for i in sorted(denovoMutsOffRepeats):
        print(i[0],i[1],i[2],i[3],sep='\t',file=fout)

mutsDict={}
for i in sorted(denovoMutsOffRepeats):
        mutsDict[i[0]+':'+i[1]]=i[0],i[1],i[2],i[3]

for file in rawCalls:
    linesDict={}
    with open(file,'r') as fin, open(file+'.profile.group','w') as fout1, open(file+'.profile.group.snp','w') as fout2, open(file+'.profile.group.indel','w') as fout3, open(file+'.profile.private','w') as fout4, open(file+'.profile.private.snp','w') as fout5, open(file+'.profile.private.indel','w') as fout6:

        print('#CHROM','POS','REF','ALT',sep='\t',file=fout1);print('#CHROM','POS','REF','ALT',sep='\t',file=fout2);print('#CHROM','POS','REF','ALT',sep='\t',file=fout3)

        for line in fin:
            if line.startswith('#'):
                print(line.strip(),file=fout4); print(line.strip(),file=fout5); print(line.strip(),file=fout6)

            else:
                linesDict[line.split('\t')[0]+':'+line.split('\t')[1]]=line.split('\t')[0],line.split('\t')[1],line.split('\t')[2],line.split('\t')[3],line.split('\t')[4],line.split('\t')[5],line.split('\t')[6],line.split('\t')[7],line.split('\t')[8],line.split('\t')[9],line.split('\t')[10].strip()

        for i in mutsDict:
# for each key in the final mutations dictionary, if it is found in the lines dictionary of this file, print the four values of the 'muts' key to a group list of variants for this file, and if not found, print the first three values and repeat the third one (Ref, so means no change or the position not called as mutated for this tumor), so produce outputs of equal length for all tumors with both thier private mutations and mutations in other members of the group (patient) listed properly. In case the 'muts' key IS found in the 'lines' dictionary, which means that this is a private (!=unique) mutation for this tumor, print all values (full line) also to a separate 'private variants' file, also separately for snpos and indels.

            if i in linesDict:
                print(mutsDict[i][0],mutsDict[i][1],mutsDict[i][2],mutsDict[i][3],sep='\t',file=fout1)
                print(linesDict[i][0],linesDict[i][1],linesDict[i][2],linesDict[i][3],linesDict[i][4],linesDict[i][5],linesDict[i][6],linesDict[i][7],linesDict[i][8],linesDict[i][9],linesDict[i][10],sep='\t',file=fout4)
                if len(mutsDict[i][2])==len(mutsDict[i][3])==1:
                    print(mutsDict[i][0],mutsDict[i][1],mutsDict[i][2],mutsDict[i][3],sep='\t',file=fout2)
                    print(linesDict[i][0],linesDict[i][1],linesDict[i][2],linesDict[i][3],linesDict[i][4],linesDict[i][5],linesDict[i][6],linesDict[i][7],linesDict[i][8],linesDict[i][9],linesDict[i][10],sep='\t',file=fout5)
                else:
                    print(mutsDict[i][0],mutsDict[i][1],mutsDict[i][2],mutsDict[i][3],sep='\t',file=fout3)
                    print(linesDict[i][0],linesDict[i][1],linesDict[i][2],linesDict[i][3],linesDict[i][4],linesDict[i][5],linesDict[i][6],linesDict[i][7],linesDict[i][8],linesDict[i][9],linesDict[i][10],sep='\t',file=fout6)

            else:
                print(mutsDict[i][0],mutsDict[i][1],mutsDict[i][2],mutsDict[i][2],sep='\t',file=fout1)
                if len(mutsDict[i][2])==len(mutsDict[i][3])==1:
                    print(mutsDict[i][0],mutsDict[i][1],mutsDict[i][2],mutsDict[i][2],sep='\t',file=fout2)
                else:
                    print(mutsDict[i][0],mutsDict[i][1],mutsDict[i][2],mutsDict[i][2],sep='\t',file=fout3)
