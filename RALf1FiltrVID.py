import numpy as np
from operator import itemgetter
import time as tm
import RALF1FilterX as XFilter
import sys
import lfib1340 
from scipy import stats as scp
import win32api,win32process,win32con
from random import sample 
#from scipy.signal import savgol_filter
           
priorityclasses = [win32process.IDLE_PRIORITY_CLASS,
               win32process.BELOW_NORMAL_PRIORITY_CLASS,
               win32process.NORMAL_PRIORITY_CLASS,
               win32process.ABOVE_NORMAL_PRIORITY_CLASS,
               win32process.HIGH_PRIORITY_CLASS,
               win32process.REALTIME_PRIORITY_CLASS]  

DETERM=0.62
NNQRandm=512
QRandm_=0

def RandomQ(Nfx,NQRandm_=0):
    global NQRandm
    global QRandm_
    
    if not NQRandm_==0:
        NQRandm=NQRandm_
    
    KK=3e6
    liiX=np.zeros(Nfx,float)
    pp=0
    while pp<0.55:
        for ii in range(3):
            if NQRandm>=NNQRandm:
                QRandm_=np.asarray(range(NNQRandm),float)
                NQRandm=0
            try:                
                z=(QRandm_[NQRandm]+1)/KK           
                atim0=tm.time()        
                tm.sleep(z) 
                atim=tm.time()     
                dd=int(((atim-atim0)/z-1)/1000)
                zz=np.asarray(sample(list(range(Nfx)),Nfx),float)/KK
                lfib1340.LFib1340(dd).shuffle(zz)  
                lfib1340.LFib1340(int(2*dd/(1+np.sqrt(5)))).shuffle(QRandm_)
                
                if NQRandm>0:
                    liiX=liiX+zz
                NQRandm=NQRandm+1
            except:
                NQRandm=NNQRandm                

        k2, pp = scp.normaltest(liiX)
            
    rr2=[[],[]]
    rr2[0]= np.asarray(liiX[0:Nfx],float)
    rr2[1]= np.asarray(range(Nfx),int)
    m=[[rr2[j][l] for j in range(len(rr2))] for l in range(len(rr2[0]))]         
    m.sort(key=itemgetter(0))                  
    rr2=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]  
    liiXX=np.asarray(rr2[1],int)
    return liiXX


def filterFourierQ(arxx,arb,NNew,NChan,key=0):  
    Nfl=int(len(arb)/NChan)
    Nnl=NNew
    
    ar_=np.zeros(Nnl,float)
    farx=np.zeros(Nnl,float)
    
    az=int(np.floor(Nfl/Nnl))-1
    
    gg0=0    
    for l in range(NChan):        
        for i in range(az):
            for j in range(Nnl):
                ar_[j]=arb[Nfl-(az-i+1)*Nnl+j+Nfl*l]
                gg0=gg0+ar_[j]*ar_[j]
            ar_=abs(np.fft.fft(ar_))
            for j in range(Nnl):
                farx[j]=max(farx[j],ar_[j])
    gg0=np.sqrt(gg0)/(NChan*az*Nnl)
    
    gg=0
    arxr=arb.copy()
    for l in range(NChan):       
        farxx=np.fft.fft(arxx[Nfl-Nnl+Nfl*l:Nfl+Nfl*l])    
        mfarxx=np.abs(farxx)   
        amn=np.mean(mfarxx[1:])
        if not key<0: 
            srmfarxx=DETERM*amn
        else:
            srmfarxx=0.0001*amn
                    
        farxxx=np.zeros(Nnl,complex)     
        for j in range(1,Nnl):
            if mfarxx[j]>srmfarxx:
                farxxx[j]=farxx[j]/mfarxx[j]*farx[j] 
            
        if key>0:
            farxxx[1]=0*farxxx[1]
            farxxx[len(farxxx)-1]=0*farxxx[len(farxxx)-1]            
                   
        arxr[Nfl-Nnl+Nfl*l:Nfl+Nfl*l]=np.fft.ifft(farxxx).real[0:Nnl] 
#        arxr[Nfl-Nnl+Nfl*l:Nfl+Nfl*l]=arxr[Nfl-Nnl+Nfl*l:Nfl+Nfl*l]-arxr[Nfl-Nnl+Nfl*l]+arxr[Nfl-Nnl+Nfl*l-1]
        gg=gg+np.std(arxr[Nfl-Nnl+Nfl*l:Nfl+Nfl*l])
        
    return arxr

def RALF1FilterQ(dQ2):    
    Np=len(dQ2)
    Nf=len(dQ2[0])
       
    SdQ=np.mean(dQ2,0)  
    sSdQ=np.std(np.asarray(SdQ,float))
    for i in range(Np):
        SdQj_ = np.std(np.asarray(dQ2[i] - SdQ,float))
        SdQj__ = np.std(np.asarray(dQ2[i],float))            
        if SdQj__ >0. and sSdQ>0.:
            dQ2[i] = np.asarray(dQ2[i] +SdQ * ((SdQj_ - sSdQ)/ sSdQ ),np.float16)
        else:
            dQ2[i]=np.zeros(Nf,np.float16)        
    return dQ2
 
import warnings

def RALF1Calculation(arr_bx,arr_c,Nf,NNew,NNew0,NChan,D,Nhh,iProc):
    global NQRandm
    global QRandm_
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    Koe=1e-4 
    sz=Nf*NChan
    MM=2
    Nzz=int(Nhh*2)
    
    Ndel=MM
    NCh=int(np.ceil(sz/Ndel)) 
    Ndel0=MM
    NCh0=int(np.ceil(sz/Ndel0))   
          
    arr_b=np.asarray(arr_bx,float)
    
    arr_bZ=[]

    #arr_b[0]=0
    for l in range(NChan):
        #arr_b[l]=arr_bx[l]-arr_bx[l-1]        
        arr_bZ.append(arr_b[0+Nf*l:Nf-NNew+Nf*l])    
    arr_bZ=np.asarray(arr_bZ,np.float16)

    R4=np.zeros(sz,np.float16)  
    for l in range(NChan):
        R4[Nf-NNew+Nf*l:Nf+Nf*l]=D*Koe*2  
    
    hh=0          
    AMX=np.zeros((Nhh,sz),float)
    AMN=np.zeros((Nhh,sz),float) 
    arr_bbbxxx1=np.zeros((Nhh,sz),float)
    arr_bbbxxx2=np.zeros((Nhh,sz),float)
    max_dd1=np.zeros((Nhh+1,sz),np.float16)
    min_dd2=np.zeros((Nhh+1,sz),np.float16)
    r2=np.zeros((Nhh+1,sz),np.float16)
    NQRandm=512
    QRandm_=np.asarray(range(NQRandm),float)
    while hh<Nhh:  
        if hh==0:
            mn=np.mean(arr_bZ)         
            r2[hh]=np.asarray(arr_b,np.float16)
            for l in range(NChan):                
                r2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=mn   
            r2[hh]=r2[hh]-mn               
       
        liix=np.zeros((sz,sz),int) 
        dQ3_0=np.zeros((sz,sz),np.float16)
        mDD=np.zeros((sz,sz),np.float16)  
        
        aa=RandomQ(sz)  
        liiB=np.concatenate((aa,aa,aa))  
        aa=RandomQ(sz)  
        liiC=np.concatenate((aa,aa,aa))   
        aa=RandomQ(sz) 
        liiD=np.concatenate((aa,aa,aa))           
        for i in range(sz):    
            liix[i]=liiB[liiD[i+liiC[hh]]:sz+liiD[i+liiC[hh]]].copy()
            dQ3_0[i]=r2[hh,liix[i]].copy()
            mDD[i]=R4[liix[i]].copy()     
         
        dQ3=dQ3_0.copy() 
               
        ##########################################       
        sseq_=dQ3_0.reshape(sz*sz)*(1/(mDD.reshape(sz*sz)<D*Koe))  
        sseq_=np.asarray(list(filter(lambda x: abs(x)!= np.Inf, sseq_)),float) 
        sseq_=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, sseq_)),float)  
        
        dQ3mx=np.zeros((sz,sz),np.float16)-np.Inf
        dQ3mn=np.zeros((sz,sz),np.float16)+np.Inf  
        AsrXMx_=0
        AsrXMn_=0               
        AsrXMx=dQ3mx=np.zeros((sz,sz),np.float16)-np.Inf
        AsrXMn=dQ3mn=np.zeros((sz,sz),np.float16)+np.Inf 

        aa=RandomQ(sz) 
        NumFri0=np.concatenate((aa, aa, aa))  
        aa=RandomQ(sz) 
        NumFri0_=np.concatenate((aa, aa, aa)) 
        aa=RandomQ(sz) 
        rR0=np.concatenate((aa, aa, aa))   
        aa=RandomQ(sz) 
        liiC=np.concatenate((aa, aa, aa)) 
        aa=RandomQ(sz)  
        r5=aa.copy()
        r5=D*((r5/np.std(r5))/2+Koe*2) 
        r5=np.concatenate((r5, r5))
        aa=RandomQ(sz) 
        ss4=np.concatenate((aa, aa, aa, aa))                         
        zz=0  
        WW=0    
        while zz<Nzz and WW>-Nhh: 
            NumFri=NumFri0_[NumFri0[ss4[zz]]:NumFri0[ss4[zz]]+2*sz].copy()
            NumFri_=NumFri0[NumFri0_[ss4[zz]]:NumFri0_[ss4[zz]]+2*sz].copy()
            rR=rR0[liiC[ss4[zz]]:liiC[ss4[zz]]+2*sz].copy()
            rR_=rR0[liiC[ss4[len(ss4)-2*sz+zz]]:liiC[ss4[len(ss4)-2*sz+zz]]+2*sz].copy()
            kk=-1
            xxx=0
            while kk <Ndel-1 and xxx==0:  
                kk=kk+1                      
                ii=int(kk*NCh)
                k=-1
                while k<Ndel0-1 and xxx==0:     
                    k=k+1                       
                    i=int(k*NCh0) 
                    dQ4=np.zeros((NCh,NCh0),float)
                    mDD4=np.zeros((NCh,NCh0),float)
                    mDD4_=np.zeros((NCh,NCh0),float)
                    mDD4_A=np.zeros((NCh,NCh0),float) 
                    mDD4_B=np.zeros((NCh,NCh0),float)                                    
                    for ll in range(NCh0):
                        dQ4[:,ll]=(dQ3[NumFri[ii:ii+NCh],NumFri_[i+ll]])*1.
                        mDD4[:,ll]=(1-(mDD[NumFri[ii:ii+NCh],NumFri_[i+ll]]<D*Koe))*1.
                        mDD4_[:,ll]=(r5[rR_[ss4[ll]+zz]:rR_[ss4[ll]+zz]+NCh])*1.
                        mDD4_A[:,ll]=(r5[rR[ss4[ll]+zz]:rR[ss4[ll]+zz]+NCh]*(dQ4[:,ll]< D*Koe))*1.
                        mDD4_B[:,ll]=(r5[rR[ss4[ll]+zz]:rR[ss4[ll]+zz]+NCh]*(dQ4[:,ll]>-D*Koe))*1.
                        
                    mDD4_=(mDD4_*(1-mDD4))*1.
                    mDD4_=(mDD4_-np.mean(mDD4_))*2                                
                    P=np.zeros(3,float)
                                  
                    nNxA=sum(sum(mDD4<D*Koe))
                    nNxA_=sum(sum(1-mDD4<D*Koe))                    
                    if nNxA>nNxA_ and nNxA_>0:  
                        seqA=dQ4.reshape(NCh*NCh0)*(1/(mDD4.reshape(NCh*NCh0)<D*Koe)) 
                        seqA=np.asarray(list(filter(lambda x: abs(x)!= np.Inf, seqA)),float) 
                        seqA=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqA)),float)

                        mNxA=sum(sum(dQ4*(mDD4<D*Koe)))/nNxA                        
                        amNxA=np.sqrt(sum(sum((dQ4-mNxA)*(dQ4-mNxA)*(mDD4<D*Koe))))/nNxA
                        dQ4_=mNxA
                        
                        dQ4=dQ4-dQ4_
                        dQ4_A= np.asarray(XFilter.RALF1FilterX(  dQ4*(1-(dQ4<0))-mDD4_A+mDD4_,len(dQ4),len(dQ4[0]),1,0)+
                                          XFilter.RALF1FilterX(  dQ4*(1-(dQ4<0))-mDD4_B+mDD4_,len(dQ4),len(dQ4[0]),1,0),np.float16)
                        dQ4_B=-(np.asarray(XFilter.RALF1FilterX(-dQ4*(1-(dQ4>0))-mDD4_B+mDD4_,len(dQ4),len(dQ4[0]),1,0)+
                                           XFilter.RALF1FilterX(-dQ4*(1-(dQ4>0))-mDD4_A+mDD4_,len(dQ4),len(dQ4[0]),1,0),np.float16))
                        
                        dQ4=(dQ4_A+dQ4_B)/2                        
                        mNxB=sum(sum(dQ4*(mDD4<D*Koe)))/nNxA 
                        amNxB=np.sqrt(sum(sum((dQ4-mNxB)*(dQ4-mNxB)*(mDD4<D*Koe))))/nNxA                           
                        P[2]=mNxA
                        P[1]=mNxB
                        P[0]=amNxB/amNxA
                                                
                        dQ4_A= dQ4_+(dQ4_A-dQ4_A*(dQ4_A>0))
                        dQ4_B= dQ4_+(dQ4_B-dQ4_B*(dQ4_B<0))
                        dQ4=(dQ4_A+dQ4_B)/2                        
                        dQ4_A=dQ4.copy()
                        dQ4_B=dQ4.copy() 
                        dQ4_A=(dQ4_A-P[1])/P[0] +P[2]
                        dQ4_B=(dQ4_B-P[1])/P[0] +P[2]  
                    
                        seqB=dQ4.reshape(NCh*NCh0)*(1/(mDD4.reshape(NCh*NCh0)<D*Koe)) 
                        seqB=np.asarray(list(filter(lambda x: abs(x)!= np.Inf, seqB)),float) 
                        seqB=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqB)),float)
                        
                        if 100*scp.pearsonr(seqA,seqB)[0]>0:                                      
                            for ll in range(NCh0):
                                dQ3mx[NumFri[ii:ii+NCh],NumFri_[i+ll]]=np.maximum(dQ3mx[NumFri[ii:ii+NCh],NumFri_[i+ll]],dQ4_A[:,ll])
                                dQ3mn[NumFri[ii:ii+NCh],NumFri_[i+ll]]=np.minimum(dQ3mn[NumFri[ii:ii+NCh],NumFri_[i+ll]],dQ4_B[:,ll])
                       
                        else:
                            xxx=1
                            
                    else:     
                        xxx=1
                        
            if xxx==0:     
                AsrXMx=np.maximum(AsrXMx,dQ3mx)
                AsrXMn=np.minimum(AsrXMn,dQ3mn)
                AsrXMx_=(AsrXMx_*zz+AsrXMx)/(zz+1)
                AsrXMn_=(AsrXMn_*zz+AsrXMn)/(zz+1)
                    
                WW=0                                    
                zz=zz+1
            else:
                ss4[0:len(ss4)-1]=ss4[1:]
                ss4[len(ss4)-1]=ss4[0]
                WW=WW-1
                
        if WW==0:
            dQ4=(AsrXMx_+AsrXMn_)/2
            
            sseq=dQ4.reshape(sz*sz)*(1/(mDD.reshape(sz*sz)<D*Koe))  
            sseq=np.asarray(list(filter(lambda x: abs(x)!= np.Inf, sseq)),float) 
            sseq=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, sseq)),float)  
            WW=WW-1               
            dQ3=dQ3_0*(mDD<D*Koe)+(dQ4)*(np.asarray(1,np.float16)-(mDD<D*Koe))
            if not sum(sum(np.isnan(dQ3)))>0:
                try:
                    if 100*scp.pearsonr(sseq,sseq_)[0]>20:
                        WW=WW+1
                except:
                    WW=WW
        
        hh0=hh
            
        if not WW<0:
            aMx_=0
            aMn_=0
            aMx=np.zeros(sz,float)-np.Inf
            aMn=np.zeros(sz,float)+np.Inf
            for i in  range(sz):
                aMx[liix[i]]=np.maximum(aMx[liix[i]],dQ3[i])
                aMn[liix[i]]=np.minimum(aMn[liix[i]],dQ3[i])
                aMx_=aMx.copy()#(aMx_*i+aMx)/(i+1)
                aMn_=aMn.copy()#(aMn_*i+aMn)/(i+1)
                # aMx_=aMx
                # aMn_=aMn
                    
            ann=sum(np.isnan(aMx_ + aMn_))
            if ann==0: 
                if hh==0: 
                    AMX[hh]=aMx_.copy()
                    AMN[hh]=aMn_.copy()  
                    arr_bbbxxx1[hh]=AMX[hh].copy()
                    arr_bbbxxx2[hh]=AMN[hh].copy()
                else:
                    AMX[hh]=np.maximum(AMX[hh-1],aMx)
                    AMN[hh]=np.minimum(AMN[hh-1],aMn)
                    arr_bbbxxx1[hh]=AMX[hh].copy()#(arr_bbbxxx1[hh-1]*hh+AMX[hh])/(hh+1)
                    arr_bbbxxx2[hh]=AMN[hh].copy()#(arr_bbbxxx2[hh-1]*hh+AMN[hh])/(hh+1)                  

                ann=1                 
                dd1=filterFourierQ(arr_bbbxxx1[hh],arr_b,NNew,NChan)
                dd2=filterFourierQ(arr_bbbxxx2[hh],arr_b,NNew,NChan)                 
                     
                if sum(np.abs(dd1+dd2)==np.Inf)==0:  
                    for l in range(NChan):
                        max_dd1[hh,Nf*l:Nf-NNew+Nf*l]=r2[hh,Nf*l:Nf-NNew+Nf*l]+mn
                        min_dd2[hh,Nf*l:Nf-NNew+Nf*l]=r2[hh,Nf*l:Nf-NNew+Nf*l]+mn
                        dd1[Nf-NNew+Nf*l:Nf+Nf*l]=dd1[Nf-NNew+Nf*l:Nf+Nf*l]/np.std(dd1[Nf-NNew+Nf*l:Nf-NNew0+Nf*l])*np.std(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])
                        dd1[Nf-NNew+Nf*l:Nf+Nf*l]=dd1[Nf-NNew+Nf*l:Nf+Nf*l]-np.mean(dd1[Nf-NNew+Nf*l:Nf-NNew0+Nf*l])+np.mean(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])
                        dd2[Nf-NNew+Nf*l:Nf+Nf*l]=dd2[Nf-NNew+Nf*l:Nf+Nf*l]/np.std(dd2[Nf-NNew+Nf*l:Nf-NNew0+Nf*l])*np.std(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])
                        dd2[Nf-NNew+Nf*l:Nf+Nf*l]=dd2[Nf-NNew+Nf*l:Nf+Nf*l]-np.mean(dd2[Nf-NNew+Nf*l:Nf-NNew0+Nf*l])+np.mean(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])
                    
                    for l in range(NChan):
                        if hh==0:
                             max_dd1[hh,Nf-NNew+Nf*l:Nf+Nf*l]=dd1[Nf-NNew+Nf*l:Nf+Nf*l].copy()
                             min_dd2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=dd2[Nf-NNew+Nf*l:Nf+Nf*l].copy()
                        else:
                             max_dd1[hh,Nf-NNew+Nf*l:Nf+Nf*l]=(max_dd1[hh,Nf-NNew+Nf*l:Nf+Nf*l]*hh+np.maximum(max_dd1[hh-1,Nf-NNew+Nf*l:Nf+Nf*l],dd1[Nf-NNew+Nf*l:Nf+Nf*l]))/(hh+1)
                             min_dd2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=(min_dd2[hh,Nf-NNew+Nf*l:Nf+Nf*l]*hh+np.minimum(min_dd2[hh-1,Nf-NNew+Nf*l:Nf+Nf*l],dd2[Nf-NNew+Nf*l:Nf+Nf*l]))/(hh+1)
                                                                       
                    hh=hh+1
                    ann=0         
                    r2[hh]=(max_dd1[hh-1]+min_dd2[hh-1])/2                         
                    mn=np.mean(r2[hh])
                    r2[hh]=r2[hh]-mn

                    if hh==Nhh:
                        dd=filterFourierQ(r2[hh],arr_b,NNew,NChan)                    
                        if sum(abs(dd)==np.Inf)==0:
                            r2[hh]=dd+mn
                            for l in range(NChan):  
                                r2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=r2[hh,Nf-NNew+Nf*l:Nf+Nf*l]/np.std(r2[hh,Nf-NNew+Nf*l:Nf-NNew0+Nf*l])*np.std(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])
                                r2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=r2[hh,Nf-NNew+Nf*l:Nf+Nf*l]-np.mean(r2[hh,Nf-NNew+Nf*l:Nf-NNew0+Nf*l])+np.mean(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])
                            anamef="fralf.tmp"
                            fo = open(anamef, "w")
                            fo.write(str(iProc)+'\n')
                            fo.close()
                            return r2[hh]
                        else:
                            return r2[hh]/0  
                     

        if hh0==hh:
            if hh>1:
                hh=hh-2
            else:
                return r2[hh]/0                            

def RALf1FiltrQ(args):    
    pid = win32api.GetCurrentProcessId()
    handle = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS, True, pid)
    win32process.SetPriorityClass(handle, priorityclasses[1])
    NChan=int(args[1])
    NNew=int(args[2])
    Nhh=int(args[3])
    Nf=int(len(args)-4) 
       
    arr_bb=[]    
    for i in range(Nf):
        arr_bb.append(args[4+i])
    arr_bb=np.asarray(arr_bb,float)
            
    arr_b=arr_bb.copy() 
    Nf=int(arr_b.size/NChan)
    
    NNew0=int(NNew*1.1) 
    arr_c=[]
    for l in range(NChan):
        arr_c.append(arr_b[Nf-NNew0+Nf*l:Nf-NNew+Nf*l]) 
    arr_c=np.asarray(arr_c,np.float16)    

    arr_bZ=[]
    for l in range(NChan):
        arr_bZ.append(arr_b[0+Nf*l:Nf-NNew0+Nf*l])   
    arr_bZ=np.asarray(arr_bZ,float)
    D=np.std(arr_bZ)
    arr_b=np.asarray(arr_bb,np.float16)
    
    while 1==1: 
        hh=0
        ann=0
        arr_bbx=[]
        Nch=0
        Koef=np.zeros(Nhh,float)
        KoefA=np.zeros(Nhh,float)
        while hh<Nhh:
            if hh<Nhh:       
                arr_bbbxxx=RALF1Calculation(arr_b,arr_c,Nf,NNew0,NNew,NChan,D,2*Nhh,args[0])
                if (sum(np.abs(arr_bbbxxx)==np.Inf)==0 and sum(np.isnan(arr_bbbxxx))==0):                
                    Nf_=int(NNew0*1.8)
                    NNew_=Nf_-NNew0
                    arr_bbbxxx_=np.zeros(Nf_*NChan,np.float16)
                    for l in range(NChan):
                        dd_=arr_bbbxxx[Nf-1+Nf*l:Nf-NNew0+Nf*l:-1].copy()
                        arr_bbbxxx_[0+Nf_*l:Nf_+Nf_*l]=(np.concatenate((dd_,np.ones(Nf_-len(dd_),float)*dd_[len(dd_)-1])))  

                    NNew0_=int(NNew_*1.1) 
                    arr_c_=[]
                    for l in range(NChan):
                        dd_=arr_bbbxxx[Nf-NNew0_+(NNew0_-NNew_)+Nf*l:Nf-NNew0_+Nf*l:-1].copy()
                        arr_c_.append(dd_)  
                    arr_c_=np.asarray(arr_c_,np.float16)    
                    
                    arr_bbbxxx_y=RALF1Calculation(arr_bbbxxx_,arr_c_,Nf_,NNew0_,NNew_,NChan,D,2*Nhh,args[0])
                    if (sum(np.abs(arr_bbbxxx_y)==np.Inf)==0 and sum(np.isnan(arr_bbbxxx_y))==0): 
                        arr_bbbxxx_yy=[]
                        
                        for l in range(NChan):
                            dd_=arr_bbbxxx_y[Nf_-1+Nf_*l:Nf_-NNew0_+Nf_*l:-1].copy()
                            arr_bbbxxx_yy.append(dd_) 
                            if l==0:
                                mm1=arr_b[Nf-NNew0-len(dd_):Nf-NNew0].copy()
                                mm2=arr_bbbxxx_yy[l].copy()
                            else:
                                mm1=np.concatenate((mm1,arr_b[Nf-NNew0-len(dd_):Nf-NNew0]))
                                mm2=np.concatenate((mm2,arr_bbbxxx_yy[l]))                
                        
                        ann=(sum(np.abs(mm1)==np.Inf)>0 + sum(np.isnan(mm1))>0+
                             sum(np.abs(mm2)==np.Inf)>0 + sum(np.isnan(mm2))>0)
                        
                        if ann==0 and len(mm1)>1 and len(mm1)==len(mm2): 
                            mm1=mm1-sum(mm1)/len(mm1)
                            mm2=mm2-sum(mm2)/len(mm1)
                       
                            if np.std(mm1)>0 and np.std(mm2)>0:
                                anamef="fralf_.tmp"
                                fo = open(anamef, "w")
                                fo.write(str(args[0])+'\n')
                                fo.close() 
                                KoefA[hh]=100*(scp.pearsonr(mm1,mm2)[0])
                                #mm1=mm1*np.std(mm2)/np.std(mm1)                       
                                Koef[hh]=-np.std(mm1-mm2)
                                arr_bbx.append(arr_bbbxxx) 
                                hh=hh+1
                
            if hh==Nhh:            
                arr_bbx=np.asarray(arr_bbx,np.float16)
                r2s=np.zeros((2,Nhh),float)
                r2s[0]= np.asarray(Koef,float)
                r2s[1]= np.asarray(range(Nhh),float)
                m=[[r2s[j][l] for j in range(len(r2s))] for l in range(len(r2s[0]))]         
                m.sort(key=itemgetter(0))                  
                r2s=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]  
                Nch=int(r2s[1][Nhh-1])
                print(KoefA)
                if np.isnan(KoefA[Nch]):
                    KoefA[Nch]=0            
                if KoefA[Nch]>0:
                    for l in range(NChan):
                        arr_b[Nf-NNew+Nf*l:Nf+Nf*l]=arr_bbx[Nch][Nf-NNew+Nf*l:Nf+Nf*l].copy()    
                    #arr_b=filterFourierQ(arr_b,arr_b,NNew,NChan,0)
                    return arr_b

if __name__ == '__main__':
    RALf1FiltrQ(sys.argv)
