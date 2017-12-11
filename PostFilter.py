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



unprod = pdSummary[pdSummary['V-DOMAIN Functionality'] != 'productive']
print(unprod['Sequence ID'])
IgRep = pdSummary.groupby('Trimmed V-GENE')['Sequence ID'].nunique()
def absolute_value(val):
    a = np.round(val/100.* IgRep.sum(), 0)
    a = int(a)
    return a
IgRep.plot.pie(autopct = absolute_value, fontsize = 12, figsize=(6,6))
plt.show()
