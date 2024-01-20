import urllib.request
import time as tm

def getcsv(WhO,xYears,wrkdir):
    tm0=int(tm.time())
    tm0=tm0-24*60*60
    tm1=tm0-xYears*365*24*60*60
    
    for i in range(len(WhO)):
        nm=WhO[i]
        anm="https://query1.finance.yahoo.com/v7/finance/download/%s?period1=%s&period2=%s&interval=1d&events=history&includeAdjustedClose=true"%(nm,tm1,tm0)   
        
        webUrl  = urllib.request.urlopen(anm)
        data = webUrl.read().decode()
        with open( wrkdir+nm+'.csv', 'w' ) as output:
            output.write( data )