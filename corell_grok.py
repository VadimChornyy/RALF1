import pandas as pd
import numpy as np
import csv
import scipy.interpolate as interp
import yfinance as yf

Lo = 1
# Original WhO list with duplicate ICP (handled in processing)
WhO = [
    'XRP', "AAVE","ANKR","AVAX","AXS","BAL","BTC","BCH","ADA","CGLD","LINK","COMP",
    "CRV","ICP","MANA","DOGE","EOS","ETH","ETC","GRT","LTC","LRC","MKR","NU",
    "OMG","ORN","MATIC","SHIB","SKL","SOL","XLM","SNX","XTZ","UNI"
]


WhO= ['MANA', 'ICP', 'DASH', 'ADA','SNX', 'SHIB', 'LTC', 'UNI', 'MKR', 'EOS']
wrkdir = r""
GetSCV = 1
aDecm = 7
andecim=1
Lengt0=200



def decimat(adat_):
    if Lo:
        adat_ = np.log(adat_)
    k = 0
    adat__ = np.zeros(int(len(adat_) / aDecm), float)
    for i in range(int(len(adat_) / aDecm)):
        adat__[k] = np.median(adat_[i * aDecm:i * aDecm + aDecm])
        k = k + 1

    arrray=adat__[1:len(adat__)].copy()
    
    indices = np.arange(len(arrray))
    non_nan_mask = ~np.isnan(arrray)
    nan_mask = np.isnan(arrray)
    
    interpolated_values = np.interp(indices[nan_mask], indices[non_nan_mask], arrray[non_nan_mask])
    
    arrray[nan_mask] = interpolated_values
    
    if Lo:
        return np.exp(arrray)
    else:
        return (arrray)
    
def getcsv(WhO, xYears, wrkdir):
    for i in range(len(WhO)):
        nm=WhO[i]
        try:         
            df = yf.download("%s-USD"% (nm), period="5y", interval="1d")
            arrr = df['Close']
            arrr=arrr[::-1]
            arrr.to_csv(wrkdir + nm + '.csv', index=False)
        except:
            nm = nm
if GetSCV:
    getcsv(WhO, 5, wrkdir)

def loaddata(aLengt, ticker1, key):
    adat_ = []
    with open(wrkdir + ticker1 + '.csv', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        dat = []
        i = 0
        for row in spamreader:
            if i > 0:
                dat.append(row[0])
            i = i + 1
    dat = dat[len(dat)::-1]
    try:
        dat = np.asarray(dat, float)
        arrr = np.asarray(dat, float)
    except:
        dat = dat
        arrr = np.zeros(2, int)
        
    if andecim>=1:
        sz=len(arrr)
        xold=np.asarray(range(sz),int)
        f=interp.interp1d(xold, arrr,'cubic')
        xnew=np.asarray(range((sz-1)*andecim),int)/andecim
        arrr=f(xnew)
        sz=len(arrr)
        arrr=arrr[sz-aLengt*aDecm:sz]
        # %varexp --plot arrr 
        
    return arrr, adat_

arrrxxR = []
nams = []
lenar = []

for i in range(len(WhO)):
    i0 = 0
    while i0 < 6 and not i0 < 0:
        try:
            arrrxx1, adat1_ = loaddata(Lengt0, WhO[i], 1)
            arrrxxR.append(arrrxx1)
            lenar.append(len(arrrxx1))
            nams.append(WhO[i])
            i0 = -1
        except:
            i0 = i0 + 1
            #tm.sleep(10)

llar = int(np.median(np.asarray(lenar, int)))
nnams_ = []
aaer = []
for i in range(len(nams)):
    if len(arrrxxR[i]) >= llar:
        aer = decimat(arrrxxR[i])
        aaer.append(aer[len(aer) - int(llar / aDecm) + 1:].copy())
aaer = np.asarray(aaer, float)

arrrxxR_ = []
ii = 0
for i in range(len(nams)):
    if len(arrrxxR[i]) >= llar:
        # arrrxxR_.append(np.diff(np.log(aaer[ii])))
        zzzz=np.log(aaer[ii])
        zzzz=(zzzz-np.mean(zzzz))/np.std(zzzz)
        zzzz_=zzzz
        zzzz_[0]=0
        zzzz_[1:]=np.diff(zzzz)
        arrrxxR_.append(zzzz+zzzz_)        
        nnams_.append(nams[ii])
        ii = ii + 1
arrrxxR_ = np.asarray(arrrxxR_, float)

srarr = np.median((arrrxxR_), axis=0)
arrrxxRR_ = np.asarray(arrrxxR_, float) * 0
for i in range(len(nnams_)):
    arrrxxRR_[i] = (arrrxxR_[i] - srarr)

# Compute the correlation matrix
correlation_mat = np.corrcoef(arrrxxRR_)

# Sort based on average absolute correlation
avg_correlations = np.mean(np.abs(correlation_mat), axis=0)
sorted_indices = np.argsort(avg_correlations)  # Sort ascending (lowest to highest correlation)
correlation_mat_sorted = correlation_mat[sorted_indices][:, sorted_indices]
nnams_sorted = [nnams_[i] for i in sorted_indices]

# Extract the most independent cryptocurrencies (top 25% lowest average correlation)
num_independent = max(1, len(nnams_sorted) // 4)  # At least 1, up to 25% of the list
most_independent_indices = sorted_indices[:num_independent]
most_independent_coins = [nnams_sorted[i] for i in range(num_independent)]
print("Most independent cryptocurrencies:", most_independent_coins)

# Visualize the sorted correlation matrix
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(font_scale=2)
plt.figure(figsize=(12, 10))
sns.heatmap(correlation_mat_sorted, 
            xticklabels=nnams_sorted,
            yticklabels=nnams_sorted,
            cmap='RdPu',
            vmin=-0.2, vmax=1.0,
            cbar_kws={'label': 'Correlation'})
plt.title('Correlation Matrix of Cryptocurrencies (Sorted by Average Correlation)\n')
plt.show()