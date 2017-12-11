import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('classic')
pdSummary = pd.read_csv('1_Summary.csv')#, index_col = 0)
trimmedgene = []
for i in range(len(pdSummary['Sequence number'] - 1)):
    conv = int(i)
    typecheck = isinstance(pdSummary['V-GENE and allele'][conv], float)
    if typecheck == False:
        splstr = pdSummary['V-GENE and allele'][conv].split(' ')
        trimmedgene.append(splstr[1])
    else:
        trimmedgene.append('')
pdSummary['Trimmed V-GENE'] = trimmedgene


ProdFilter = pdSummary[pdSummary['V-DOMAIN Functionality'] == 'productive']
#Filter Unproductive/Bad Seq
varbool = ProdFilter['AA JUNCTION'].duplicated()
invbool = ~varbool #We need to invert so unique Abs are "True"
UniqAbs = ProdFilter[invbool]
 #If a duplicate, remove from dataset
SeqDuplic = 0
s = UniqAbs['AA JUNCTION'].str.len().sort_values().index
print(UniqAbs[['Trimmed V-GENE','Sequence number','AA JUNCTION']].reindex(s))
while SeqDuplic != '':
    SeqDuplic = input('Enter Seq Number of Suspected Duplicate to Remove:')
    if SeqDuplic != '':
        xin = int(SeqDuplic)
        UniqAbs = UniqAbs[UniqAbs['Sequence number'] != xin]
        s = UniqAbs['AA JUNCTION'].str.len().sort_values().index
        print(UniqAbs[['Trimmed V-GENE','Sequence number', 'AA JUNCTION']].reindex(s))
    else:
        break
IgRep = UniqAbs.groupby('Trimmed V-GENE')['Sequence ID'].nunique()
print(IgRep)
def absolute_value(val):
    a = np.round(val/100.* IgRep.sum(), 0)
    a = int(a)
    return a
IgRep.plot.pie(autopct = absolute_value, fontsize = 12, figsize=(6,6))
plt.show()
