import concurrent.futures
import pylab as plt
import numpy as np
import pandas as pd
from PIL import Image 
import cv2 as cv
#import time as tm
import hickle as hkl

from scipy import stats as scp
import dill 
#from scipy.signal import savgol_filter
import scipy.interpolate as interp

from RALf1FiltrVID import filterFourierQ
from RALf1FiltrVID import RALf1FiltrQ
from RALf1FiltrVID import RandomQ
import RALF1FilterX as XFilter 

wrkdir = r"c:\Work\\"
aname='lena-Geo'
Lengt=1000
Ngroup=3
Nproc=2*Ngroup#*(mp.cpu_count())
Lo=1
aTmStop=3
NIt=3
NIter=100
DT=0.45
Nf_K=3

def fig2img ( fig ):
    fig.savefig(wrkdir +'dynamic.png',dpi=150,transparent=False,bbox_inches = 'tight')
    frame=Image.open(wrkdir +'dynamic.png')
    # fig.canvas.draw()
    # frame=Image.frombytes('RGB', fig.canvas.get_width_height(),
    #                        fig.canvas.tostring_rgb())
    return frame
    
def loaddata(aLengt,key):
    adat_=[]
    arrr=[]
    excel_data_df = pd.read_excel('lena.xls', sheet_name='1')
    
    dat=np.asarray(excel_data_df, float)
    f=interp.interp1d(dat[:,0], dat[:,1],'cubic')
    xnew=np.asarray(range(3,463),float)[::4]
    arrr.append(f(xnew))
            
    return arrr,adat_

if __name__ == '__main__':   
    global NQRandm

    ImApp=[]
        
    arrrxx,adat_=loaddata(Lengt,1)
    arrrxx=np.asarray(arrrxx,float)
    
    esz=len(arrrxx)
    arr_rezDzRez=[[] for j in range(esz)]
    kkk=0
    out=0       
            
    w=0
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
            NNew=int(Nf*0.6)  
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
            all_RezN=np.zeros((Ngroup,NIter,Nf),float)
            all_RezM=np.zeros((Ngroup,NIter,Nf),float)
            all_RezNM=np.zeros((Ngroup,NIter,Nf),float)
            all_RezMM=np.zeros((Ngroup,NIter,Nf),float)
            argss=[[0] for j in range(Nproc)]    
    
            hh0=0
            hhh=0
            hhh_=0
                    
            Koef_=[]
            ZZ=0
            key=0
            try:
                dill.load_session(wrkdir + aname+".ralf")
            except:    
                fig = plt.figure()
                axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
                axes_ = fig.add_axes([0, 0, 0.3, 0.3])     
                axes__ = fig.add_axes([0.4, 0, 0.3, 0.3])
                axes.plot(ar0, 'r.')
                axes.plot(arr_z, 'go-')  #cut data used for model
                axes.text(4, 4, 'Reflection on Action of Lorentz Forces-1, #2011612714  \n\n Course = %s'%(aname),
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes_.transAxes,color='blue', fontsize=14)        

                frame=fig2img(fig)  
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
            if hhh==int(NIter/1):
                if hhh_<aTmStop-1:
                    nnn=int(nn*0.5)
                    aTmStop=6
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
                    all_RezN=np.zeros((Ngroup,NIter,Nf),float)
                    all_RezM=np.zeros((Ngroup,NIter,Nf),float)
                    all_RezNM=np.zeros((Ngroup,NIter,Nf),float)
                    all_RezMM=np.zeros((Ngroup,NIter,Nf),float)
                    hhh_=hhh_+1
                else:
                    hhh_=hhh_+1
                    ZZ=1
            
            if ZZ==0:                  
                try:
                    [hhha,Arr_BBB]=(hkl.load(wrkdir + aname+".rlf1"))       
                    Arr_AAA[:,0:int(Nproc*hhha/Ngroup),:]=Arr_BBB[:,0:int(Nproc*hhha/Ngroup),:].copy()
                except:            
                    hhha=hhh-1
                
                if hhh>=hhha: 
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
                    
                    for iGr in range(Ngroup):
                        Arr_AAA[iGr][hhh*int(Nproc/Ngroup):(hhh+1)*int(Nproc/Ngroup)]=(
                            arezAMx[int(iGr*(Nproc/Ngroup)):int((iGr+1)*(Nproc/Ngroup))]).copy()
                    
                WrtTodr=0
                if hhh>=hhha-1:    
                    WrtTodr=1
                    aDur=4
                        
                dNIt=10
                arr_RezM=  np.zeros((Ngroup,Nf),float)                  
                for hhhb in range(max(0,hhh-int(NIter/dNIt)),hhh+1):
                    for iGr in range(Ngroup):            
                        arr_RezM[iGr]=(np.amax(Arr_AAA[iGr][max(0,hhhb-int(NIter/dNIt)):(hhhb+1)*int(Nproc/Ngroup),:],axis = 0)+
                                           np.amin(Arr_AAA[iGr][max(0,hhhb-int(NIter/dNIt)):(hhhb+1)*int(Nproc/Ngroup),:],axis = 0))/2
                        
                        all_RezN[iGr][hhhb]=arr_RezM[iGr].copy()
                        nI=(hhhb+1)-max(0,hhhb-int(NIter/dNIt))                        
                        if nI>1:
                            NQRandm=512
                            D=np.std(arr_RezM[iGr])                     
                            aa=RandomQ(Nf,512)                        
                            ss4=np.concatenate((aa, aa, aa))
                            DD=[]
                            for hhhc in range(nI):
                                DD.append(ss4[hhhc:hhhc+Nf])
                            DD=np.asarray(DD,float)                              
                            DD=(DD/np.std(DD)+1e-6)*D/2   
                            
                            mn=np.mean(all_RezN[iGr][max(0,hhhb-int(NIter/dNIt)):hhhb+1])
                            dd=(all_RezN[iGr][max(0,hhhb-int(NIter/dNIt)):hhhb+1]-mn)
                            
                            aa=RandomQ(Nf,512)                        
                            ss4=np.concatenate((aa, aa, aa))
                            liix=np.zeros((nI,Nf),int)
                            for i in range(nI):  
                                liix[i]=ss4[i:i+Nf].copy()
                                dd[i]=dd[i][liix[i]].copy()                                
                                
                            for ii in range(2):
                                dd1=dd[:,int(ii*Nf/2):int((ii+1)*Nf/2)].copy()
                                ddA=dd1*(1-(dd1<0))
                                ddA=ddA+DD[:,int(ii*Nf/2):int((ii+1)*Nf/2)]*(ddA==0)
                                ddB=-dd1*(1-(dd1>0))
                                ddB=ddB+DD[:,int(ii*Nf/2):int((ii+1)*Nf/2)]*(ddB==0)
                                dd[:,int(ii*Nf/2):int((ii+1)*Nf/2)]=mn+(XFilter.RALF1FilterX(  ddA,len(ddA),len(ddA[0]),1,0)-
                                         XFilter.RALF1FilterX(  ddB,len(ddB),len(ddB[0]),1,0))/2 
                        
                            for i in range(nI):
                                dd[i][liix[i]]=dd[i].copy()
                                
                            dd=(np.amax(dd,axis=0)+np.amin(dd,axis=0))/2
                            arr_RezM[iGr]=(np.maximum(arr_RezM[iGr],dd)+np.minimum(arr_RezM[iGr],dd))/2
                            
                        #np.mean(all_RezN[iGr][max(0,hhhb-int(NIter/6)):hhhb+1,:],axis = 0) 
                        all_RezM[iGr][hhhb]=arr_RezM[iGr].copy() 
                        arr_RezM[iGr]=(np.amax(all_RezM[iGr][max(0,hhhb-int(NIter/dNIt)):hhhb+1,:],axis = 0)+
                            np.amin(all_RezM[iGr][max(0,hhhb-int(NIter/dNIt)):hhhb+1,:],axis = 0))/2 
                        
                        if Lo:
                            arr_RezM[iGr]=filterFourierQ(arr_RezM[iGr],np.log(arr_z),NNew,1,1)
                            arr_RezM[iGr][0:Nf-NNew]=np.log(ar0[0:Nf-NNew])  
                            arr_RezM[iGr][Nf-NNew:]=(arr_RezM[iGr][Nf-NNew:]-arr_RezM[iGr][Nf-NNew]) +np.log(ar0[Nf-NNew-1])
                        else:
                            arr_RezM[iGr]=filterFourierQ(arr_RezM[iGr],arr_z,NNew,1,1)
                            arr_RezM[iGr][0:Nf-NNew]=ar0[0:Nf-NNew].copy()
                            arr_RezM[iGr][Nf-NNew:]=(arr_RezM[iGr][Nf-NNew:]-arr_RezM[iGr][Nf-NNew]) +(ar0[Nf-NNew-1])

                        all_RezNM[iGr][hhhb]=arr_RezM[iGr].copy() 
                        arr_RezM[iGr]=(np.amax(all_RezNM[iGr][max(0,hhhb-int(NIter/dNIt)):hhhb+1,:],axis = 0)+
                            np.amin(all_RezNM[iGr][max(0,hhhb-int(NIter/dNIt)):hhhb+1,:],axis = 0))/2                         
    
                        all_RezMM[iGr][hhhb]=arr_RezM[iGr].copy() 
                        arr_RezM[iGr]=np.mean(all_RezMM[iGr][0:hhhb+1],axis = 0)
                        #np.mean(all_RezMM[iGr][max(0,hhhb-int(NIter/2)):hhhb+1,:],axis = 0) 
                
                arr_rezBz=(np.mean(arr_RezM, axis=0)+np.mean(arr_RezM, axis=0))/2            
                    
                if Lo:
                    arr_rezBz[Nf-NNew:]=(arr_rezBz[Nf-NNew:]-arr_rezBz[Nf-NNew]) +np.log(ar0[Nf-NNew-1])                     
                    for iGr in range(Ngroup):                    
                        arr_RezM[iGr][Nf-NNew:]=(arr_RezM[iGr][Nf-NNew:]-arr_RezM[iGr][Nf-NNew]) +np.log(ar0[Nf-NNew-1])
                else:
                    arr_rezBz[Nf-NNew:]=(arr_rezBz[Nf-NNew:]-arr_rezBz[Nf-NNew]) +ar0[Nf-NNew-1]
                    for iGr in range(Ngroup):
                        arr_RezM[iGr][Nf-NNew:]=(arr_RezM[iGr][Nf-NNew:]-arr_RezM[iGr][Nf-NNew]) +ar0[Nf-NNew-1]
                                           
                if Lo:   
                    for iGr in range(Ngroup):    
                        arr_RezM[iGr]=np.exp(arr_RezM[iGr]) 
                    arr_rezBz=np.exp(arr_rezBz) 
                    
                mm1=ar0[Nf-NNew:].copy()                            
                mm2=arr_rezBz[Nf-NNew:len(ar0)].copy()   
                if np.std(mm1)>0 and np.std(mm2)>0:
                    if hhh>=hhha:                               
                        hkl.dump([hhh+1,Arr_AAA], wrkdir + aname+".rlf1")                        
                    Koef_.append(100*scp.pearsonr(mm1,mm2)[0])                               
                    fig = plt.figure()
                    axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
                    axes_ = fig.add_axes([0, 0, 0.3, 0.3]) 
                    axes__ = fig.add_axes([0.4, 0, 0.3, 0.3])                     
                    axes.plot(ar0, 'ro-', alpha=0.1)
                    axes.plot(arrr, 'rx-')
                    for iGr in range(Ngroup):
                        axes.plot(arr_RezM[iGr],linewidth=3.,alpha=0.2)
                    axes.plot(arr_rezBz,'yx-',linewidth=4.,alpha=0.5)
                    axes.text(-0.1, 4, '%s  '%(hhh+1),
                            verticalalignment='bottom', horizontalalignment='left',
                            transform=axes_.transAxes,color='black', fontsize=24) 
                    axes.text(4, 4, 'Reflection on Action of Lorentz Forces-1, #2011612714  \n\n Course = %s'%(aname),
                            verticalalignment='bottom', horizontalalignment='right',
                            transform=axes_.transAxes,color='blue', fontsize=14)     
                    gint=np.polyfit(mm1,mm2,1)
                    axes_.plot(mm1,gint[1]+gint[0]*mm1,'y.',linewidth=2.)                    
                    axes_.plot(mm1,mm2, 'ok', markersize=3, alpha=0.1) 
                    axes_.text(0.2, 0.6, '%d'%int(np.floor(np.asarray(Koef_,float)[::-1][0])),
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes_.transAxes,color='green', fontsize=30)
                    
                    axes__.text(1.8, 0.6, 'Dunning –Kruger\n effect',
                            verticalalignment='bottom', horizontalalignment='center',
                        transform=axes_.transAxes,color='green', fontsize=14)  
                    axes__.plot(Koef_,'y',linewidth=2.)
                    
                    frame=fig2img(fig)  
                    cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)        
                    gray_sz1=min(gray_sz1,len(cimg[0]))
                    gray_sz2=min(gray_sz2,len(cimg))
                    ImApp.append(frame)
                    if WrtTodr>0:
                        out = cv.VideoWriter(wrkdir + aname+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
                        for icl in range(len(ImApp)):
                            cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
                            out.write(cimgx[0:gray_sz2,0:gray_sz1,:]) 
                        out.release()
                        del(out)
                    plt.show()
                    hhh=hhh+1
                    #arr_z=arr_rezBz.copy()
                else:
                    try:
                        dill.load_session(wrkdir + aname+".ralf")
                    except:
                        hh0=hh0                  
                hh0=hh0+1
                if WrtTodr>0:
                    dill.dump_session(wrkdir + aname+".ralf")   
                if hh0==2*NIter:
                    hhh=NIter 
                print (hhh+10000*hh0)
                
                    
#        df = pd.DataFrame(arr_rezBz)
#        df.to_excel (wrkdir +r'export_traces.xlsx', index = None, header=False) 
        
        mm1=arrr[len(arrr)-int(len(arrr)/2):].copy()                            
        mm2=arr_rezBz[len(arrr)-int(len(arrr)/2):len(arrr)].copy()  
        if np.std(mm1)>0 and np.std(mm2)>0:
            Koef=100*scp.pearsonr(mm1,mm2)[0] 
            
        fig = plt.figure()
        axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
        axes_ = fig.add_axes([0, 0, 0.3, 0.3])
        axes__ = fig.add_axes([0.4, 0, 0.3, 0.3]) 
        axes.plot(arrr, 'ro-')
        axes.plot(arr_rezBz, 'cx-', alpha=0.5) #predicted data
        axes.text(4, 4, 'Reflection on Action of Lorentz Forces-1, #2011612714  \n\n Course = %s'%(aname),
                verticalalignment='bottom', horizontalalignment='right',
                transform=axes_.transAxes,color='blue', fontsize=14)     
        gint=np.polyfit(mm1,mm2,1)
        axes_.plot(mm1,gint[0]+gint[1]*mm1,'y',linewidth=2.) 
        axes_.plot(mm1,mm2, 'ok', markersize=3, alpha=0.1)    
        axes_.text(0.2, 0.6, '%d'%int(np.floor(np.asarray(Koef_,float)[::-1][0])),
            verticalalignment='bottom', horizontalalignment='right',
            transform=axes_.transAxes,color='green', fontsize=30)    
        
        axes__.text(1.8, 0.6, 'Dunning –Kruger\n effect',
            verticalalignment='bottom', horizontalalignment='center',
                    transform=axes_.transAxes,color='green', fontsize=14)  
        axes__.plot(Koef_,'y',linewidth=2.)
        
        frame=fig2img(fig) 
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