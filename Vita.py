import concurrent.futures
import pylab as plt
import numpy as np
import pandas as pd
from PIL import Image 
import cv2 as cv
#import time as tm

from scipy import stats as scp
import dill 
#from scipy.signal import savgol_filter

from RALf1FiltrVID import filterFourierQ
from RALf1FiltrVID import RALf1FiltrQ
import scipy.interpolate as interp

wrkdir = r"c:\Work\\"
Lengt=1000
Ngroup=3
Nproc=2*Ngroup#*(mp.cpu_count())
Lo=0
aTmStop=3
NIt=3
NIter=20
DT=0.25
Nf_K=3
    
def loaddata(aLengt,key):
    adat_=[]
    arrr=[]
    excel_data_df = pd.read_excel('bok.xls', sheet_name='1')
    
    dat=np.asarray(excel_data_df, float)
    f=interp.interp1d(dat[:,0], dat[:,1],'cubic')
    xnew=np.asarray(range(75,915),float)[::5]
    arrr.append(f(xnew))
    
    excel_data_df = pd.read_excel('bok.xls', sheet_name='2')
    
    dat=np.asarray(excel_data_df, float)
    f=interp.interp1d(dat[:,0], dat[:,1],'cubic')

    arrr.append(f(xnew))
        
    return arrr,adat_

if __name__ == '__main__':   
    ImApp=[]
    aname='PSS-Geo'
        
    arrrxx,adat_=loaddata(Lengt,1)
    arrrxx=np.asarray(arrrxx,float)
    
    esz=len(arrrxx)
    arr_rezDzRez=[[] for j in range(esz)]
    kkk=0
    out=0       
            
    while kkk <esz:  
        try:
            dill.load_session(wrkdir + aname+".ralf")
        except:
            aname=aname
            arrr=np.asarray(arrrxx[kkk]).copy()  
    
            arrr=np.asarray(arrr,float)    
            Lengt=len(arrr)
            Nf=Lengt
            
            nn=int(Nf*DT)             
            NNew=int(Nf/2)  
            Nf=Nf+nn        
            ar0=np.asarray(arrr[0:])           
            
            arr_z=np.zeros(Nf,float)
            arr_z[0:Nf-NNew]=arrr[0:Nf-NNew].copy()
            arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]
              
            adat0=''
                                   
            Arr_AAA=np.zeros((Ngroup,int(Nproc*NIter/Ngroup),Nf),float) 
            arr_rezBz=np.zeros(Nf,float)
            mn=np.mean(arr_z[0:Nf-NNew])
            
            all_rezAz=np.zeros((NIter,Nf),float)
            arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]  
            all_RezM=np.zeros((Ngroup,NIter,Nf),float)
            argss=[[0] for j in range(Nproc)]    
    
            hh0=0
            hhh=0
            hhh_=0
                    
            Koef=-1
            ZZ=0
            key=0
            TKoef=-100
            
            nnn=int(nn/2)
    
            fig = plt.figure()
            axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
            axes_ = fig.add_axes([0, 0, 0.3, 0.3])     
            axes.plot(ar0, 'r.')
            axes.plot(arr_z, 'go-')  #cut data used for model
            axes.text(4, 4, 'Course = %s, start = %s, step = %s'%(aname,adat0,''),
                    verticalalignment='bottom', horizontalalignment='right',
                    transform=axes_.transAxes,color='blue', fontsize=14)        
            fig.savefig(wrkdir +'dynamic.png',dpi=300,transparent=False,bbox_inches = 'tight')
            frame=Image.open(wrkdir +'dynamic.png')
            ImApp.append(frame)
            cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)        
            gray_sz1=len(cimg[0])
            gray_sz2=len(cimg)
            aDur=2
            fourcc = cv.VideoWriter_fourcc(*'MP4V')
            out = cv.VideoWriter(wrkdir + aname+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
            for icl in range(len(ImApp)):
                cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
                out.write(cimgx[0:gray_sz2,0:gray_sz1,:]) 
            out.release()
            plt.show()
            del(out)
    
        while hhh_<aTmStop and not key == 13: 
            Aprocess=[]
            if hhh==NIter:
                if hhh_<aTmStop-1:
                    TKoef=-100
                    hh0=0
                    hhh=0
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
                    Arr_AAA=np.zeros((Ngroup,int(Nproc*NIter/Ngroup),Nf),float) 
                    all_rezAz=np.zeros((NIter,Nf),float)
                    arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]  
                    all_RezM=np.zeros((Ngroup,NIter,Nf),float)
                    hhh_=hhh_+1
                else:
                    hhh_=hhh_+1
                    ZZ=1
            if ZZ==0:                
                if Lo:
                    arr_A=np.log(arr_z)
                else:
                    arr_A=arr_z.copy()
                    
                Asr=np.mean(arr_A[0:Nf-NNew])
                arr_A=arr_A-Asr
                Klg=np.power(10,np.floor(np.log10(np.max(abs(arr_A)))))
                arr_A=arr_A/Klg
       
                program =wrkdir + "RALF1FiltrX_lg.py"
                NChan=1
                for iProc in range(Nproc):
                    argss[iProc]=["%s"%iProc, "%s"%NChan, "%s"%NNew, "%s"%NIt]#"%s"%(iProc+1)]
                    for i in range(Nf):
                        argss[iProc].append(str("%1.6f"%(arr_A[i])))
                
                ssss=int(Nf-NNew+(len(ar0)-(Nf-NNew)))
                
                arezAMx=[]
                # for iProc in range(Nproc):
                #     arezAMx.append(RALf1FiltrQ(argss[iProc]))

                with concurrent.futures.ThreadPoolExecutor(max_workers=Nproc) as executor:
                    future_to = {executor.submit(RALf1FiltrQ, argss[iProc]) for iProc in range(Nproc)}
                    for future in concurrent.futures.as_completed(future_to):                
                        arezAMx.append(future.result())
                del(future)                        
                del(executor)
                
                arezAMx= np.asarray(arezAMx,float)*Klg+Asr
                arr_RezM=  np.zeros((Ngroup,Nf),float)
                for iGr in range(Ngroup):
                    Arr_AAA[iGr][hhh*int(Nproc/Ngroup):(hhh+1)*int(Nproc/Ngroup)]=(
                        arezAMx[int(iGr*(Nproc/Ngroup)):int((iGr+1)*(Nproc/Ngroup))]).copy()
                        
                    arr_RezM[iGr]=(np.mean(Arr_AAA[iGr][0:(hhh+1)*int(Nproc/Ngroup),:],axis = 0)+
                        np.mean(Arr_AAA[iGr][0:(hhh+1)*int(Nproc/Ngroup),:],axis = 0))/2

                    if Lo:
                        arr_RezM[iGr]=filterFourierQ(arr_RezM[iGr],np.log(arr_z),NNew,1)
                        arr_RezM[iGr][0:Nf-NNew]=np.log(ar0[0:Nf-NNew])                         
                    else:
                        arr_RezM[iGr]=filterFourierQ(arr_RezM[iGr],arr_z,NNew,1)
                        arr_RezM[iGr][0:Nf-NNew]=ar0[0:Nf-NNew].copy()
                        
                    all_RezM[iGr][hhh]=arr_RezM[iGr].copy() 
                    arr_RezM[iGr]=np.mean(all_RezM[iGr][max(0,hhh-int(NIter/2)):hhh+1,:],axis = 0)                  
                    
                arr_rezBz=(np.amax(arr_RezM, axis=0)+np.amin(arr_RezM, axis=0))/2
                
                P=np.zeros(3,float)                    
                if Lo:
                    P[2]=np.mean(np.log(ar0[Nf-NNew:ssss]))
                    P[1]=np.mean(arr_rezBz[Nf-NNew:ssss])                                         
                    P[0]=np.std(arr_rezBz[Nf-NNew:ssss])/np.std(np.log(ar0[Nf-NNew:ssss]))
    
                else: 
                    P[2]=np.mean(ar0[Nf-NNew:ssss])
                    P[1]=np.mean(arr_rezBz[Nf-NNew:ssss])                                            
                    P[0]=np.std(arr_rezBz[Nf-NNew:ssss])/np.std((ar0[Nf-NNew:ssss]))
                    
                arr_rezBz[Nf-NNew:]=(arr_rezBz[Nf-NNew:]-P[1])/P[0] +P[2]                     
                for iGr in range(Ngroup):
                    arr_RezM[iGr][Nf-NNew:]=(arr_RezM[iGr][Nf-NNew:]-P[1])/P[0] +P[2]
                    
                if Lo:   
                    for iGr in range(Ngroup):    
                        arr_RezM[iGr]=np.exp(arr_RezM[iGr]) 
                    arr_rezBz=np.exp(arr_rezBz) 
                    
                mm1=ar0[Nf-NNew:].copy()                            
                mm2=arr_rezBz[Nf-NNew:len(ar0)].copy()   
                if np.std(mm1)>0 and np.std(mm2)>0:
                    Koef=100*scp.spearmanr(mm1,mm2)[0]
                else:
                    Koef=-2  
                if (Koef+100)>=(TKoef+100):  
                    TKoef=Koef                                  
                    fig = plt.figure()
                    axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
                    axes_ = fig.add_axes([0, 0, 0.3, 0.3])   
                    axes.plot(ar0, 'ro-', alpha=0.1)
                    axes.plot(arrr, 'rx-')
                    for iGr in range(Ngroup):
                        axes.plot(arr_RezM[iGr],linewidth=3.,alpha=0.2)
                    axes.plot(arr_rezBz,'cx-', alpha=0.5)
                    axes.text(4, 4, 'Course = %s, start = %s, step = %s'%(aname,adat0,''),
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
                    del(out)
                    hhh=hhh+1
                    #arr_z=arr_rezBz.copy()
                else:
                    try:
                        dill.load_session(wrkdir + aname+".ralf")
                    except:
                        hh0=hh0                  
                hh0=hh0+1
                dill.dump_session(wrkdir + aname+".ralf")   
                if hh0==2*NIter:
                    hhh=NIter 
                print (hhh+10000*hh0)
                
                    
#        df = pd.DataFrame(arr_rezBz)
#        df.to_excel (wrkdir +r'export_traces.xlsx', index = None, header=False) 
        
        mm1=arrr[len(arrr)-int(len(arrr)/2):].copy()                            
        mm2=arr_rezBz[len(arrr)-int(len(arrr)/2):len(arrr)].copy()  
        if np.std(mm1)>0 and np.std(mm2)>0:
            Koef=100*scp.spearmanr(mm1,mm2)[0] 
            
        fig = plt.figure()
        axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
        axes_ = fig.add_axes([0, 0, 0.3, 0.3])
        axes.plot(arrr, 'ro-')
        axes.plot(arr_rezBz, 'cx-', alpha=0.5) #predicted data
        axes.text(4, 4, 'Course = %s, start = %s, step = %s'%(aname,adat0,''),
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
        for icl in range(10):
            ImApp.append(frame)
        out = cv.VideoWriter(wrkdir + aname+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
        for icl in range(len(ImApp)):
            cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
            out.write(cimgx[0:gray_sz2,0:gray_sz1,:])       
        out.release()
        plt.show()
        del(out)
        kkk=kkk+1    
