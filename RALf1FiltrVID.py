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

DETERM=0.5

def RandomQ(Nfx,NQRandm_=0):
    QRandm_=np.asarray(range(512),float)
    NQRandm=0  
    KK=3e6
    liiX=np.zeros(Nfx,float)
    pp=0
    while pp<0.55:
        for ii in range(3):
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
                NQRandm=0                

        k2, pp = scp.skewtest(liiX)
            
    rr2=[[],[]]
    rr2[0]= liiX.copy()
    rr2[1]= np.asarray(np.linspace(0,Nfx-1,Nfx),int)
    m=[[rr2[j][l] for j in range(len(rr2))] for l in range(len(rr2[0]))]         
    m.sort(key=itemgetter(0))                  
    rr2=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]  
    liiXX=np.asarray(rr2[1],int)
    return liiXX


def filterFourierQ(arxx,arb,NNew,NChan,key=0):  
    Nfl=int(len(arb)/NChan)
    Nnl=NNew
    
    ar_=np.zeros(Nnl,float)
    farx=np.zeros(2*Nnl,float)-np.Inf
    
    az=int(np.floor(Nfl/Nnl))-1
    
    gg0=0    
    for l in range(NChan):        
        for i in range(az):
            ar_=arb[Nfl-(az-i+1)*Nnl+Nfl*l:Nfl+Nnl-(az-i+1)*Nnl+Nfl*l].copy()
            gg0=gg0+np.sum(ar_*ar_)
            ar__=abs(np.fft.fft(np.concatenate((ar_,ar_)))) 
            farx=np.maximum(farx,ar__)
    gg0=np.sqrt(gg0)/(NChan*az*Nnl)
    
    gg=0
    arxr=arb.copy()
    for l in range(NChan):      
        farxx=np.fft.fft(np.concatenate((arxx[Nfl*(l+1)-Nnl:Nfl*(l+1)],
                                         arxx[Nfl*(l+1)-Nnl:Nfl*(l+1)])))    
        mfarxx=np.abs(farxx)   
        A1=np.mean(mfarxx[1:])+DETERM*np.std(mfarxx[1:])
        A2=np.mean(mfarxx[1:])-DETERM*np.std(mfarxx[1:])     
        farxxx=np.zeros(2*Nnl,complex)     
        for j in range(1,2*Nnl):
            if mfarxx[j]>A1:
                farxxx[j]=farxx[j]/mfarxx[j]*farx[j] 
            if mfarxx[j]<A2:
                farxxx[j]=0*farxx[j]
            
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
        asta=np.Inf
        # dQ3_0=dQ3_0.reshape((sz*sz))  
        # asta=dQ3_0[0]
        # dQ3[0]=0
        # dQ3[1:]=np.diff(dQ3_0)
        # dQ3=dQ3.reshape((sz,sz))
        # dQ3_0=dQ3.copy()
                 
        ##########################################       
        sseq_=dQ3_0.reshape(sz*sz)*(1/(mDD.reshape(sz*sz)<D*Koe))  
        sseq_=np.asarray(list(filter(lambda x: abs(x)!= np.Inf, sseq_)),float) 
        sseq_=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, sseq_)),float)  
        
        AsrXMx_=np.zeros((Nzz,sz,sz),np.float16)-np.Inf
        AsrXMn_=np.zeros((Nzz,sz,sz),np.float16)+np.Inf 
        dQ3mx=np.zeros((sz,sz),np.float16)-np.Inf
        dQ3mn=np.zeros((sz,sz),np.float16)+np.Inf 

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
        r5=D*((r5/np.std(r5))/2+Koe*2)*2
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
                        try:                        
                            seqB=dQ4.reshape(NCh*NCh0)*(1/(mDD4.reshape(NCh*NCh0)<D*Koe)) 
                            seqB=np.asarray(list(filter(lambda x: abs(x)!= np.Inf, seqB)),float) 
                            seqB=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqB)),float)
                            
                            if 100*scp.pearsonr(seqA,seqB)[0]>0:                                      
                                for ll in range(NCh0):
                                    dQ3mx[NumFri[ii:ii+NCh],NumFri_[i+ll]]=np.maximum(dQ3mx[NumFri[ii:ii+NCh],NumFri_[i+ll]],dQ4_A[:,ll])
                                    dQ3mn[NumFri[ii:ii+NCh],NumFri_[i+ll]]=np.minimum(dQ3mn[NumFri[ii:ii+NCh],NumFri_[i+ll]],dQ4_B[:,ll])
                           
                            else:
                                xxx=1
                        except:
                            xxx=1
                    else:     
                        xxx=1
                        
            if xxx==0:     
                if zz>0:
                    AsrXMx_[zz]=np.maximum(AsrXMx_[zz-1],dQ3mx)
                    AsrXMn_[zz]=np.minimum(AsrXMn_[zz-1],dQ3mn)
                else:
                    AsrXMx_[zz]=dQ3mx.copy()
                    AsrXMn_[zz]=dQ3mn.copy()
                    
                WW=0                                    
                zz=zz+1
            else:
                ss4[0:len(ss4)-1]=ss4[1:]
                ss4[len(ss4)-1]=ss4[0]
                WW=WW-1
                
        if WW==0:
            dd1=np.median(AsrXMx_[0:Nzz],axis=0)
            dd2=np.median(AsrXMn_[0:Nzz],axis=0)
            dQ3=(dd1+dd2)/2
            
            sseq=dQ3.reshape(sz*sz)*(1/(mDD.reshape(sz*sz)<D*Koe))  
            sseq=np.asarray(list(filter(lambda x: abs(x)!= np.Inf, sseq)),float) 
            sseq=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, sseq)),float)  
            WW=WW-1               
            dQ3=dQ3_0*(mDD<D*Koe)+(dQ3)*(np.asarray(1,np.float16)-(mDD<D*Koe))
            if not sum(sum(np.isnan(dQ3)))>0:
                try:
                    if 100*scp.pearsonr(sseq,sseq_)[0]>20:
                        WW=WW+1
                except:
                    WW=WW        
        hh0=hh
            
        if not WW<0: 
            if not asta==np.Inf:
                dd1=dd1.reshape((sz*sz))
                dd1=asta+np.cumsum(dd1)
                dd1=dd1.reshape((sz,sz))
                dd2=dd2.reshape((sz*sz))
                dd2=asta+np.cumsum(dd2)
                dd2=dd2.reshape((sz,sz))
            
            aMx=np.zeros(sz,float)-np.Inf
            aMn=np.zeros(sz,float)+np.Inf
            for i in  range(sz):
                aMx[liix[i]]=np.maximum(aMx[liix[i]],dd1[i]+mn)
                aMn[liix[i]]=np.minimum(aMn[liix[i]],dd2[i]+mn)
                    
            ann=sum(np.isnan(aMx + aMn))
            if ann==0: 
                if hh==0: 
                    AMX[hh]=aMx.copy()
                    AMN[hh]=aMn.copy()  

                else:
                    AMX[hh]=np.maximum(AMX[hh-1],aMx)
                    AMN[hh]=np.minimum(AMN[hh-1],aMn)      
                    
                ann=1                 
                dd1=filterFourierQ(AMX[hh],arr_b,NNew,NChan)
                dd2=filterFourierQ(AMN[hh],arr_b,NNew,NChan)                 
                     
                if sum(np.abs(dd1+dd2)==np.Inf)==0:  
                    for l in range(NChan):
                        dd1[Nf-NNew+Nf*l:Nf+Nf*l]=dd1[Nf-NNew+Nf*l:Nf+Nf*l]/np.std(dd1[Nf-NNew+Nf*l:Nf-NNew0+Nf*l])*np.std(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])
                        dd1[Nf-NNew+Nf*l:Nf+Nf*l]=dd1[Nf-NNew+Nf*l:Nf+Nf*l]-np.mean(dd1[Nf-NNew+Nf*l:Nf-NNew0+Nf*l])+np.mean(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])
                        dd2[Nf-NNew+Nf*l:Nf+Nf*l]=dd2[Nf-NNew+Nf*l:Nf+Nf*l]/np.std(dd2[Nf-NNew+Nf*l:Nf-NNew0+Nf*l])*np.std(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])
                        dd2[Nf-NNew+Nf*l:Nf+Nf*l]=dd2[Nf-NNew+Nf*l:Nf+Nf*l]-np.mean(dd2[Nf-NNew+Nf*l:Nf-NNew0+Nf*l])+np.mean(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])
                    
                    mn=0
                    for l in range(NChan):
                        if hh==0:
                             max_dd1[hh,Nf-NNew+Nf*l:Nf+Nf*l]=dd1[Nf-NNew+Nf*l:Nf+Nf*l].copy()
                             min_dd2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=dd2[Nf-NNew+Nf*l:Nf+Nf*l].copy()
                        else:
                             max_dd1[hh,Nf-NNew+Nf*l:Nf+Nf*l]=np.maximum(max_dd1[hh-1,Nf-NNew+Nf*l:Nf+Nf*l],dd1[Nf-NNew+Nf*l:Nf+Nf*l])
                             min_dd2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=np.minimum(min_dd2[hh-1,Nf-NNew+Nf*l:Nf+Nf*l],dd2[Nf-NNew+Nf*l:Nf+Nf*l])
                                               
                    hh=hh+1
                    ann=0  
                    dd1=np.mean(max_dd1[0:hh,Nf-NNew+Nf*l:Nf+Nf*l],axis=0)
                    dd2=np.mean(min_dd2[0:hh,Nf-NNew+Nf*l:Nf+Nf*l],axis=0)
                    asr1=np.abs(dd1-mn)>np.abs(dd2-mn)
                    asr2=np.abs(dd1-mn)<np.abs(dd2-mn)
                    r2[hh]=r2[hh-1].copy()
                    for l in range(NChan):  
                        r2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=(dd1*asr1+dd2*asr2)
                    r2[hh]=filterFourierQ(r2[hh],arr_b,NNew,NChan) 
                    r2[hh]=(r2[hh-1]*(hh-1)+r2[hh])/hh
                    for l in range(NChan):  
                        r2[hh,Nf*l:Nf-NNew+Nf*l]=r2[hh-1,Nf*l:Nf-NNew+Nf*l].copy()
                        r2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=r2[hh,Nf-NNew+Nf*l:Nf+Nf*l]/np.std(r2[hh,Nf-NNew+Nf*l:Nf-NNew0+Nf*l])*np.std(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])
                        r2[hh,Nf-NNew+Nf*l:Nf+Nf*l]=r2[hh,Nf-NNew+Nf*l:Nf+Nf*l]-np.mean(r2[hh,Nf-NNew+Nf*l:Nf-NNew0+Nf*l])+np.mean(arr_c[(NNew-NNew0)*l:NNew-NNew0+(NNew-NNew0)*l])

                    mn=0
                    for l in range(NChan):                      
                        mn=mn+np.mean(r2[hh,Nf*l:Nf-NNew+Nf*l])
                    mn=mn/NChan
                    r2[hh]=r2[hh]-mn

                    if hh==Nhh:             
                        if sum(abs(r2[hh])==np.Inf)==0:
                            r2[hh]=r2[hh]+mn
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
    IIIX=RandomQ(len(arr_bb))
    arr_b=arr_bb.copy()
    astar=arr_b[IIIX[0]]
    arr_b[0]=0
    arr_b[1:]=np.diff(arr_bb[IIIX])
    arr_b[IIIX]=arr_b.copy()
    
    Nf=int(len(arr_b)/NChan)    
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
    NumIt=2*Nhh
    Nhh0=Nhh
    while 1==1: 
        hh=0
        ann=0
        arr_bbx=np.zeros((Nhh*3,Nf),float)
        Nch=0
        Koef=np.zeros(Nhh*3,float)                
        KoefA=np.zeros(Nhh*3,float)

        while hh<Nhh:
            if hh<Nhh:       
                arr_bbbxxx=RALF1Calculation(arr_b,arr_c,Nf,NNew0,NNew,NChan,D,NumIt,args[0])
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
                    
                    arr_bbbxxx_y=RALF1Calculation(arr_bbbxxx_,arr_c_,Nf_,NNew0_,NNew_,NChan,D,NumIt,args[0])
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
                                coef=100*(scp.pearsonr(mm1,mm2)[0])
                                print('%s'%(coef))
                                KoefA[hh]=coef
         
                                #mm1=mm1*np.std(mm2)/np.std(mm1)                       
                                Koef[hh]=-np.std(mm1-mm2)
                                arr_bbx[hh]=arr_bbbxxx.copy() 
                                hh=hh+1
                
            if hh==Nhh:            
                arr_bbx_=np.asarray(arr_bbx[0:Nhh],np.float16)
                r2s=np.zeros((2,Nhh),float)
                r2s[0]= np.asarray(Koef[0:Nhh],float)
                r2s[1]= np.asarray(range(Nhh),float)
                m=[[r2s[j][l] for j in range(len(r2s))] for l in range(len(r2s[0]))]         
                m.sort(key=itemgetter(0))                  
                r2s=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]  
                Nch=int(r2s[1][Nhh-1])
                if np.isnan(KoefA[Nch]):
                    KoefA[Nch]=0            
                if KoefA[Nch]>10:
                    print(KoefA[0:Nhh])
                    REZ=arr_bbx_[Nch].copy()
                    REZ[IIIX]=astar+np.cumsum(arr_bbx_[Nch])
                    return REZ
                else:
                    if (Nhh<len(KoefA)):
                        Nhh=Nhh+1
                    else:
                        hh=0
                        Nhh=Nhh0


if __name__ == '__main__':
    RALf1FiltrQ(sys.argv)