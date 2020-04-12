import multiprocessing as mp
import concurrent.futures
import pylab as plt
import numpy as np
import urllib.request, json
import pandas as pd
from PIL import Image
import msvcrt
import cv2 as cv
import time as tm

from scipy import stats as scp
import dateutil.parser
from operator import itemgetter
from scipy.signal import savgol_filter

from RALf1FiltrVID import RALf1FiltrV,filterFourierV

api_key = 'ONKTYPV6TAMZK464' 

ticker ="BTC-USD"#"GLD"#"DJI","LOIL.L"#""BZ=F" "LNGA.MI" #"BTC-USD"#"USDUAH"#"LTC-USD"#"USDUAH"#
#interv="15min"
interv="Daily"
#url_string =  "https://www.alphavantage.co/query?function=TIME_SERIES_INTRADAY&symbol=%s&interval=%s&outputsize=full&apikey=%s"%(ticker,interv,api_key)        
url_string =  "https://www.alphavantage.co/query?function=TIME_SERIES_DAILY&symbol=%s&outputsize=full&apikey=%s"%(ticker,api_key)

#INTRADAY
#d_intervals = {"1min","5min","15min","30min","60min"}
#from scipy.signal import savgol_filter

Lengt=1000
Ngroup=1
Nproc=Ngroup*mp.cpu_count()
Lo=1
aTmStop=5
NIt=2
NIter=20
DT=0.25
Nf_K=3
    
anamef="fralf.tmp"
fo = open(anamef, "w")
fo.write(str(0)+'\n')
fo.close()  

def loaddata(aLengt,key):
    adat_=[]
    if key>0:  
        data = json.loads(urllib.request.urlopen(url_string).read().decode())['Time Series (%s)'%(interv)]
        df = pd.DataFrame(columns=['Date','Low','High','Close','Open'])
        arrr=[]
        adate=[]
        adt=[]
        for k,v in data.items():
            date = dateutil.parser.parse(k)
            data_row = [date.date(),float(v['3. low']),float(v['2. high']),
                                float(v['4. close']),float(v['1. open'])]
            df.loc[-1,:] = data_row
            rr=(np.asarray(data_row)[1]+np.asarray(data_row)[2])/2
            if rr!=0:
                adate.append(date.timestamp())
                adt.append(k)
                arrr.append(rr)
            df.index = df.index + 1
    #        if np.asarray(arrr,int).size>=aLengt:#495:1023:
    #            break
        aa=[[] for i in range(3)]
        aa[0]=adate
        aa[1]=arrr  
        aa[2]=adt
        m=aa
        aa=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]   
        aa.sort(key=itemgetter(0))
        m=aa
        aa=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]     
        ada=list(aa)[2]
        arrr=list(aa)[1]
        sz=np.asarray(arrr).size
        ln=min(sz,aLengt)
        arr=np.asarray(arrr).copy()
        arrr=[]        
        for i in range(ln-1):
            arrr.append(arr[sz-ln+i])
            adat_.append(ada[sz-ln+i])
    else:
        file=open(ticker+ '.txt','r')
        arrr=np.asarray(file.readlines(),float).copy()
        if len(arrr)>aLengt:
            arrr=arrr[len(arrr)-aLengt:]
        file.close()
    return arrr,adat_ 

if __name__ == '__main__':   
    ImApp=[]
    wrkdir = r"e:\Work\\"

    arrrxx,adat_=loaddata(Lengt,1)
    arrrxx=np.asarray(arrrxx,float)
    arrrxx_=[]
    arrrxx_.append(arrrxx)
    
    esz=len(arrrxx_)
    arr_rezDzRez=[[] for j in range(esz)]
    kkk=0
    out=0
    while kkk <esz:        
        arrr=arrrxx_[kkk].copy()    
        arrr=np.asarray(arrr,float)    
        Lengt=arrr.size  
        Nf=Lengt
        
        nn=int(Nf*DT)             
        NNew=int(Nf/2)  
        Nf=Nf+nn        
        ar0=np.asarray(arrr[0:])           
        
        arr_z=np.zeros(Nf,float)
        arr_z[0:Nf-NNew]=arrr[0:Nf-NNew].copy()
        arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]
          
        aname=ticker
        adat0=adat_[Nf-NNew]

        fig = plt.figure()
        axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
        axes_ = fig.add_axes([0, 0, 0.3, 0.3])     
        axes.plot(ar0, 'r.')
        axes.plot(arr_z, 'go-')  #cut data used for model
        axes.text(4, 4, 'Course = %s, start = %s, step = %s'%(aname,adat0,interv),
                verticalalignment='bottom', horizontalalignment='right',
                transform=axes_.transAxes,color='blue', fontsize=14)        
        fig.savefig(wrkdir +'dynamic.png',dpi=300,transparent=False,bbox_inches = 'tight')
        frame=Image.open(wrkdir +'dynamic.png')
        ImApp.append(frame)
        cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)        
        gray_sz1=len(cimg[0])
        gray_sz2=len(cimg)
        aDur=4
        fourcc = cv.VideoWriter_fourcc(*'MP4V')
        out = cv.VideoWriter(wrkdir + aname+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
        for icl in range(len(ImApp)):
            cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
            out.write(cimgx[0:gray_sz2,0:gray_sz1,:]) 
        out.release()
        plt.show()
       
        arr_rezBz=np.zeros(Nf,float)
        mn=np.mean(arr_z[0:Nf-NNew])
        arr_rezMx=  np.zeros((Ngroup,Nf),float)
        arr_rezMn=  np.zeros((Ngroup,Nf),float)
        Arr_AM=   np.asarray([[[0 for k in range(Nf)] for j in range(NIter)] for i in range(Nproc)],float)
        Arr_BM=   np.asarray([[[0 for k in range(Nf)] for j in range(NIter)] for i in range(Nproc)],float)
        all_rezAz=np.zeros((NIter,Nf),float)
        arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]  
        argss=[[0] for j in range(Nproc)]    
        hh=0
        hhh=0
        hhh_=0
        Koef=-1
        ZZ=0
        key=0
        TKoef=-100
        
        nnn=int(0.5*nn)
        arr_rezCz=np.zeros(int(nnn*(aTmStop-1)),float)
        while hhh_<aTmStop and not key == 13: 
            Aprocess=[]
            if hhh==NIter:
                if hhh_<aTmStop-1:
                    TKoef=-100
                    hhh=0
                    arr_rezCz[hhh_*nnn:(hhh_+1)*nnn]=arr_rezBz[Nf-NNew:Nf-NNew+nnn].copy()
                    arr_z=np.zeros(Nf+nnn,float)
                    arr_z[0:nnn]=arr_rezBz[0:nnn].copy()
                    arr_z[nnn:Nf]=arr_rezBz[nnn:Nf].copy()
                    arr_z[Nf:]=arr_z[Nf-1]
                    Lengt=Lengt+nnn
                    ar0=arr_z[0:Lengt-nnn].copy()
                    Nf=Nf+nnn
                    arr_rezBz=np.zeros(Nf,float)
                    mn=np.mean(arr_z[0:Nf-NNew])
                    arr_rezMx=  np.zeros((Ngroup,Nf),float)
                    arr_rezMn=  np.zeros((Ngroup,Nf),float)
                    Arr_AM=   np.asarray([[[0 for k in range(Nf)] for j in range(NIter)] for i in range(Nproc)],float)
                    Arr_BM=   np.asarray([[[0 for k in range(Nf)] for j in range(NIter)] for i in range(Nproc)],float)
                    all_rezAz=np.zeros((NIter,Nf),float)
                    arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]  
                    hhh_=hhh_+1
                else:
                    hhh_=hhh_+1
                    ZZ=1
            if ZZ==0:                                                          
                for iProc in range(Nproc):
                    arr_A=arr_z.copy()
                    if Lo:
                        arr_A=np.log(arr_z)
        
                    if hh==0 and iProc==0:
                        file=open('datapy%d.tmp'%(iProc+1),'w')
                        for i in range(Nf):
                            file.write(str(arr_A[i])+'\n')
                        file.close()
        
                        file=open("anum.tmp",'wt')
                        file.write("%d\n"%(NNew))
                        file.write("%d"%(iProc+1))
                        file.close()
        
                    program =wrkdir + "RALF1FiltrX_lg.py"
                    NChan=1
                    argss[iProc]=["python", "%s"%NChan, "%s"%NNew, "%s"%NIt]#"%s"%(iProc+1)]
                    for i in range(Nf):
                        argss[iProc].append(arr_A[i])
        
                arezAMx=[]
                # for iProc in range(Nproc):
                #     arezAMx.append(RALf1FiltrV(argss[iProc]))

                with concurrent.futures.ThreadPoolExecutor(max_workers=Nproc) as executor:
                    future_to = {executor.submit(RALf1FiltrV, argss[iProc]) for iProc in range(Nproc)}
                    for future in concurrent.futures.as_completed(future_to):                
                        arezAMx.append(future.result())
        
                Aprocess=[]
                for iProc in range(Nproc):
                    arr_B=arr_z.copy()
                    if Lo:
                        arr_B=np.log(arr_z)
        
                    arr_B=-arr_B
                    
                    program =wrkdir + "RALF1FiltrX_lg.py"
                    NChan=1
                    argss[iProc]=["python", "%s"%NChan, "%s"%NNew, "%s"%NIt]#"%s"%(iProc+1)]
                    for i in range(Nf):
                        argss[iProc].append(arr_B[i])
        
                arezBMx=[]
                with concurrent.futures.ThreadPoolExecutor(max_workers=Nproc) as executor:
                    future_to = {executor.submit(RALf1FiltrV, argss[iProc]) for iProc in range(Nproc)}
                    for future in concurrent.futures.as_completed(future_to):                
                        arezBMx.append(future.result())      
        
                arezAMx= np.asarray(arezAMx,float)
                arezBMx=-np.asarray(arezBMx,float)
                        
                Arr_AAA=np.zeros((Ngroup,int(Nproc*2*(hhh+1)/Ngroup),Nf),float)  
                for iGr in range(Ngroup):
                    for iProc in range(int(Nproc/Ngroup)):                
                        Arr_AM[iProc+int(iGr*(Nproc/Ngroup))][hhh]= np.asarray(arezAMx[iProc+int(iGr*(Nproc/Ngroup))],float)
                        Arr_BM[iProc+int(iGr*(Nproc/Ngroup))][hhh]= np.asarray(arezBMx[iProc+int(iGr*(Nproc/Ngroup))],float)
                        Arr_AAA[iGr][0*(hhh+1)+iProc*2*(hhh+1):1*(hhh+1)+iProc*2*(hhh+1)]=Arr_AM[iProc+int(iGr*(Nproc/Ngroup))][0:hhh+1].copy()
                        Arr_AAA[iGr][1*(hhh+1)+iProc*2*(hhh+1):2*(hhh+1)+iProc*2*(hhh+1)]=Arr_BM[iProc+int(iGr*(Nproc/Ngroup))][0:hhh+1].copy()
                    
                    Arr_BBB=Arr_AAA[iGr].transpose()
                    for i in range(Nf):
                        arr_rezMx[iGr][i]=max(Arr_BBB[i])
                        arr_rezMn[iGr][i]=min(Arr_BBB[i])
                        
                aMx=arr_rezMx.transpose()
                aMn=arr_rezMn.transpose()                 
                for i in range(Nf-NNew,Nf):
                    arr_rezBz[i]=max(aMx[i])+min(aMn[i])                      
                
                if Lo:           
                    arr_rezBz=filterFourierV(arr_rezBz,np.log(arr_z),NNew,1)
                else:
                    arr_rezBz=filterFourierV(arr_rezBz,arr_z,NNew,1)            
                  
                if Lo:
                    arr_rezBz[0:Nf-NNew]=np.log(ar0[0:Nf-NNew].copy())
                    # arr_rezBz[Nf-NNew:Nf]=arr_rezBz[Nf-NNew:Nf]*np.std(np.log(ar0[Nf-NNew:]))/np.std(arr_rezBz[Nf-NNew:Nf])
                    # arr_rezBz[Nf-NNew:Nf]=arr_rezBz[Nf-NNew:Nf]-np.mean(arr_rezBz[Nf-NNew:len(ar0)])+np.mean(np.log(ar0[Nf-NNew:]))
                    arr_rezBz[Nf-NNew:Nf]=arr_rezBz[Nf-NNew:Nf]-arr_rezBz[Nf-NNew]+np.log(ar0[Nf-NNew-1])
                    #arr_rezBz= savgol_filter(arr_rezBz, 11, 5)        
                    for i in range(Nf):
                        arr_rezBz[i]=np.exp(arr_rezBz[i]) 
                else:
                    arr_rezBz[0:Nf-NNew]=ar0[0:Nf-NNew].copy()  
                    # arr_rezBz[Nf-NNew:Nf]=arr_rezBz[Nf-NNew:Nf]*np.std(ar0[Nf-NNew:])/np.std(arr_rezBz[Nf-NNew:Nf])                           
                    # arr_rezBz[Nf-NNew:Nf]=arr_rezBz[Nf-NNew:Nf]-np.mean(arr_rezBz[Nf-NNew:len(ar0)])+np.mean(ar0[Nf-NNew:])
                    arr_rezBz[Nf-NNew:Nf]=arr_rezBz[Nf-NNew:Nf]-arr_rezBz[Nf-NNew]+ar0[Nf-NNew-1]
                    #arr_rezBz= savgol_filter(arr_rezBz, 11, 5)  
                arr_rezBz[0:Nf-NNew]=ar0[0:Nf-NNew].copy()  
                       
                all_rezAz[hhh]=arr_rezBz.copy()        
                all_rezAz_=all_rezAz[0:hhh+1].transpose()
                
                for i in range(Nf-NNew,Nf):
                    if Lo:
                        arr_rezBz[i]=np.exp(np.mean(np.log(all_rezAz_[i])))
                    else:
                        arr_rezBz[i]=np.mean(all_rezAz_[i]) 

                mm1=ar0[Nf-NNew:].copy()                            
                mm2=arr_rezBz[Nf-NNew:len(ar0)].copy()        
                Koef=100*scp.pearsonr(mm1,mm2)[0]
                if (Koef+100)> 0.62*(TKoef+100):                                        
                    fig = plt.figure()
                    axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
                    axes_ = fig.add_axes([0, 0, 0.3, 0.3])   
                    axes.plot(all_rezAz_,'oy',alpha=0.1)
                    axes.plot(ar0, 'ro-', alpha=0.1)
                    axes.plot(arr_rezBz,'cx-', alpha=0.5)
                    axes.text(4, 4, 'Course = %s, start = %s, step = %s'%(aname,adat0,interv),
                            verticalalignment='bottom', horizontalalignment='right',
                            transform=axes_.transAxes,color='blue', fontsize=14)        
                    axes_.plot(mm1,mm2, 'ok', markersize=3, alpha=0.1)               
                    axes_.text(0.2, 0.6, '%d'%int(np.floor(Koef)),
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes_.transAxes,color='green', fontsize=30)
                    fig.savefig(wrkdir +'dynamic.png',dpi=300,transparent=False,bbox_inches = 'tight')
                    frame=Image.open(wrkdir +'dynamic.png')
                    cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)        
                    gray_sz1=min(gray_sz1,len(cimg[0]))
                    gray_sz2=min(gray_sz2,len(cimg))
                    ImApp.append(frame)
                    out = cv.VideoWriter(wrkdir + aname+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
                    for icl in range(len(ImApp)):
                        cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
                        out.write(cimgx[0:gray_sz2,0:gray_sz1,:]) 
                    out.release()
                    plt.show()
                    hh=hh+1
                    hhh=hhh+1
                    print (hh)
                    #arr_z=arr_rezBz.copy()
                if msvcrt.kbhit():              
                    key = ord(msvcrt.getch())  
    
        arr_rezDz=[]
        file=open(wrkdir + aname+'_rezpy.txt','w')
        for i in range(len(ar0)):
            arr_rezDz.append(ar0[i])
        for i in range(len(arr_rezCz)):
            arr_rezDz.append(arr_rezCz[i])
        k=i
        for i in range(Nf-NNew,Nf):    
            arr_rezDz.append(arr_rezBz[i])            
        arr_rezDz=np.asarray(arr_rezDz,float)
        aln=len(arr_rezDz)-NNew+nnn
        for i in range(aln):     
            file.write(str(arr_rezDz[i])+'    ')
            file.write("\n")
        file.close()
        arr_rezDzRez[kkk]=arr_rezDz.copy()
        df = pd.DataFrame(arr_rezDzRez)
        df.to_excel (wrkdir +r'export_traces.xlsx', index = None, header=False) 
        
        fig = plt.figure()
        axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
        axes_ = fig.add_axes([0, 0, 0.3, 0.3])
        axes.plot(arrr, 'ro-')
        axes.plot(arr_rezDz, 'cx-', alpha=0.5) #predicted data
        axes.text(4, 4, 'Course = %s, start = %s, step = %s'%(aname,adat0,interv),
                verticalalignment='bottom', horizontalalignment='right',
                transform=axes_.transAxes,color='blue', fontsize=14)                
        axes_.plot(mm1,mm2, 'ok', markersize=3, alpha=0.1)               
        fig.savefig(wrkdir +'dynamic.png',dpi=300,transparent=False,bbox_inches = 'tight')
        frame=Image.open(wrkdir +'dynamic.png')
        cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)        
        gray_sz1=min(gray_sz1,len(cimg[0]))
        gray_sz2=min(gray_sz2,len(cimg))
        for icl in range(10):
            ImApp.append(frame)
        out = cv.VideoWriter(wrkdir + aname+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
        for icl in range(len(ImApp)):
            cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
            out.write(cimgx[0:gray_sz2,0:gray_sz1,:])       
        out.release()
        plt.show()
        kkk=kkk+1
        
        