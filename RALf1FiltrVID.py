import numpy as np
from operator import itemgetter
import time as tm
import RALF1FilterX as XFilter
import sys
import lfib1340 
from scipy import stats as scp
import scipy.interpolate as scpyi 
#from scipy.signal import savgol_filter
import win32api,win32process,win32con
import mersenne_twister as gen

def randomQ(Nfx):
    KK=1e6
    liiX=np.zeros(Nfx,float)
    for ii in range(3):            
        atim0=tm.time() 
        z=np.random.randint(Nfx)/KK
        tm.sleep(z) 
        atim=tm.time()
        delt=(atim-atim0-z)*KK            
        genera = gen.mersenne_rng(int(delt))
        z=genera.get_random_number()/KK/KK
        i=0
        while i<Nfx:
            atim0=tm.time() 
            z=genera.get_random_number()/KK/KK#np.random.randint(Nf)/KK
            tm.sleep(z) 
            atim=tm.time()
            i=i+1
            liiX[i-1]=liiX[i-1]+atim-atim0-z
            
    r2=[[],[]]
    r2[0]= (liiX[0:Nfx]).copy()
    r2[1]= np.asarray(range(Nfx),int)
    m=[[r2[j][l] for j in range(len(r2))] for l in range(len(r2[0]))]         
    m.sort(key=itemgetter(0))                  
    r2=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]  
    liiX=np.asarray(r2[1],int)
    return liiX

priorityclasses = [win32process.IDLE_PRIORITY_CLASS,
               win32process.BELOW_NORMAL_PRIORITY_CLASS,
               win32process.NORMAL_PRIORITY_CLASS,
               win32process.ABOVE_NORMAL_PRIORITY_CLASS,
               win32process.HIGH_PRIORITY_CLASS,
               win32process.REALTIME_PRIORITY_CLASS]

def filterFourierQ(arxx,arb,NNew,NChan):  
    Nfl=int(len(arb)/NChan)
    Nnl=NNew
    
    ar_=np.zeros(Nnl,float)
    farx=np.zeros(Nnl,float)
    
    az=int(np.floor(Nfl/Nnl))-1
    
    for l in range(NChan):        
        for i in range(az):
            for j in range(Nnl):
                ar_[j]=arb[Nfl-(az-i+1)*Nnl+j+Nfl*l]
            ar_=abs(np.fft.fft(ar_))
            for j in range(Nnl):
                farx[j]=max(farx[j],ar_[j])
    
    farx[0]=1e-32
    arxr=np.zeros(Nfl*NChan,float)   
    for l in range(NChan):       
        farxx=np.fft.fft(arxx[Nfl-Nnl+Nfl*l:Nfl+Nfl*l])    
        mfarxx=abs(farxx) 
        mfarxx[0]=1e-32
        srmfarxx=.62*np.mean(mfarxx[1:])
        farxxx=np.zeros(Nnl,complex)     
        for j in range(1,Nnl):
            if mfarxx[j]>srmfarxx:
                farxxx[j]=farxx[j]/mfarxx[j]*farx[j]            
            else:
                farxxx[j]=0   
        farxxx[0]=farxx[0]
        farxxx2=np.zeros(2*Nnl,complex)
        #farxxx2[0:2*Nnl:2]=farxxx.copy()
        farxxx2=farxxx.copy()
        arxr[Nfl-Nnl+Nfl*l:Nfl+Nfl*l]=np.fft.ifft(farxxx2).real[0:Nnl] 
        #arxr[Nfl-Nnl+Nfl*l:Nfl+Nfl*l]=savgol_filter(arxr[Nfl-Nnl+Nfl*l:Nfl+Nfl*l], 15, 5)
        arxr[0+Nfl*l:Nfl-Nnl+Nfl*l]=arb[0+Nfl*l:Nfl-Nnl+Nfl*l].copy() 
        #arxr[Nfl-Nnl+Nfl*l:Nfl+Nfl*l]=arxr[Nfl-Nnl+Nfl*l:Nfl+Nfl*l]-arxr[Nfl-Nnl]+arb[Nfl-Nnl-1+Nfl*l]

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

def randomX(Nfx):
    KK=3e6
    liiX=np.zeros(Nfx,float)
    pp=0
    while pp<0.055:
        for ii in range(3): 
            z=np.random.randint(Nfx)/KK           
            atim0=tm.time()        
            tm.sleep(z) 
            atim=tm.time()     
            dd=int((atim-atim0-z)*KK)
            zz=np.asarray(range(Nfx),float)/KK
            lfib1340.LFib1340(dd).shuffle(zz)   
            liiX=liiX+zz
        
        k2, pp = scp.normaltest(liiX)
            
    r2=[[],[]]
    r2[0]= np.asarray(liiX[0:Nfx],float)
    r2[1]= np.asarray(range(Nfx),int)
    m=[[r2[j][l] for j in range(len(r2))] for l in range(len(r2[0]))]         
    m.sort(key=itemgetter(0))                  
    r2=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]  
    liiXX=np.asarray(r2[1],int)
    return liiXX
  
def RALF1Calculation(arr_bx,Nf,NNew,NChan,D,Nhh):
    Koe=1e-3   
    tSp=2
    arr_bZ=[]
    arr_b=np.asarray(arr_bx,float)
    #arr_b[0]=0
    for l in range(NChan):
        #arr_b[l]=arr_bx[l]-arr_bx[l-1]        
        arr_bZ.append(arr_b[0+Nf*l:Nf-NNew+Nf*l])    
    arr_bZ=np.asarray(arr_bZ,np.float16)
    mn=np.mean(arr_bZ)
    sz=Nf*NChan
    
    hh=0
    ann=0
     
    arr_bbx=[]
    while hh<Nhh:    
        aa=randomX(sz) 
        liiB=np.concatenate((aa,aa))        
        aa=randomX((Nf+1)*NChan) 
        liiC=np.concatenate((aa,aa))         
        liiD=randomX(sz)
        liiE=randomX(sz)
        
        r4=np.zeros(sz,np.float16)
        for l in range(NChan):            
            r4[Nf-NNew+Nf*l:Nf+Nf*l]=randomX(NNew)/NNew 
            r4[Nf-NNew+Nf*l:Nf+Nf*l]=D*((r4[Nf-NNew+Nf*l:Nf+Nf*l]/np.std(r4[Nf-NNew+Nf*l:Nf+Nf*l]))/2+Koe*10) 

        r2=np.asarray(arr_b,np.float16)
        for l in range(NChan):                
            r2[Nf-NNew+Nf*l:Nf+Nf*l]=mn
    
        R4=np.asarray(r4,np.float16)       
            
        liix=np.zeros((sz,sz*tSp),np.float16)
        dQ3=np.zeros((sz,sz*tSp),np.float16)   
        mDD=np.zeros((sz,sz*tSp),np.float16)   
             
        for i in range(sz):                                                     
            r1=liiB[int(liiD[i]):sz+int(liiD[i])]                                     
            ge=scpyi.interp1d(np.asarray(range(sz),np.float16),r1)                              
            liix[i]=np.asarray(ge(np.linspace(0,sz-1,sz*tSp)),np.float16)
            liix[i]=np.asarray(liix[i]-min(liix[i]),float)
            liix[i]=np.asarray(liix[i]*(sz-1.033)/max(liix[i]),np.float16)
                          
            ge=scpyi.interp1d(np.asarray(range(sz),np.float16) ,r2,kind='linear')
            dQ3[i]=np.asarray(ge(liix[i]),np.float16)            
            for l in range(NChan):
                bb=np.asarray(0.5*np.around(2*np.arange(l+int(liiE[i]),l+sz+int(liiE[i]),sz/NNew)),int)
                bb=np.asarray(0.5*np.around(2*liiC[bb]*NNew/(sz+1)),int)
                #bb_=np.asarray(liiC[np.asarray(np.arange(l+sz+int(liiE[i]),l+int(liiE[i]),-sz/NNew),int)]*K,int)
                R4[Nf-len(bb)+Nf*l:Nf+Nf*l]=(r4[Nf-len(bb)+Nf*l+bb])#+
                                          #r4[Nf-NNew+Nf*l+bb_[0:NNew]])/np.sqrt(2)        
            ge=scpyi.interp1d(np.asarray(range(sz),np.float16) ,R4,kind='linear') 
            mDD[i]=np.asarray(ge(liix[i]),np.float16)
            
        dQ3=dQ3-mn
        zz=8
        dQ3mx=np.zeros((sz,sz*tSp),np.float16)-np.Inf
        dQ3mn=np.zeros((sz,sz*tSp),np.float16)+np.Inf
        Ndel=3#int(np.ceil(np.sqrt(sz)))
        NCh=int(np.ceil(sz/Ndel)) 
        Ndel0=3
        NCh0=int(np.ceil(sz*tSp/Ndel0)) 
        while zz>=0:        
            w=1
            while w>0:    
                try:                   
                    NumFri=randomX(sz)
                    NumFri_=randomX(sz*tSp)                 
                    NumFri=np.concatenate((NumFri, NumFri))                  
                    NumFri_=np.concatenate((NumFri_, NumFri_))  
                    for kk in range(Ndel):
                        ii=int(kk*NCh)
                        for k in range(Ndel0):
                            i=int(k*NCh0) 
                            dQ4=[]
                            mDD4=[]
                            for ll in range(NCh0):
                                dQ4.append(dQ3[NumFri[ii+0:ii+NCh],NumFri_[i+ll]])
                                mDD4.append(mDD[NumFri[ii+0:ii+NCh],NumFri_[i+ll]])
                            dQ4=np.asarray(dQ4,np.float16).transpose()
                            mDD4=np.asarray(mDD4,np.float16).transpose()
                            dQ4mn=np.mean(dQ4*(1-(np.abs(mDD4)<D*Koe)))
                            dQ4=dQ4-dQ4mn                 
                            dQ4=( XFilter.RALF1FilterX(  dQ4*(1-(dQ4<0))+mDD4,len(dQ4),len(dQ4[0]),1,0)-                    
                                  XFilter.RALF1FilterX( -dQ4*(1-(dQ4>0))+mDD4,len(dQ4),len(dQ4[0]),1,0))            
                            dQ4=dQ4+dQ4mn
                            for ll in range(NCh0):
                                dQ3mx[NumFri[ii+0:ii+NCh],NumFri_[i+ll]]=np.maximum(dQ3mx[NumFri[ii+0:ii+NCh],NumFri_[i+ll]],dQ4[:,ll])
                                dQ3mn[NumFri[ii+0:ii+NCh],NumFri_[i+ll]]=np.minimum(dQ3mn[NumFri[ii+0:ii+NCh],NumFri_[i+ll]],dQ4[:,ll])
                        
                    w=0
                except:
                    w=1                     
                zz=zz-1
                
        dQ3=(dQ3mx+dQ3mn)/2
        del(dQ3mx)
        del(dQ3mn)  
        dQ2X=dQ3.copy()
        dQ2Y=dQ2X.copy()     
        # dQ3A=dQ3.copy()        
        # dQ3B=dQ3A-dQ3A*np.asarray(dQ3A<0,int)   
        # dQ2X=XFilter.RALF1FilterX(dQ3B+mDD,sz,sz*tSp,1,0)
        # dQ3C=-(dQ3A-dQ3A*np.asarray(dQ3A>0,int))   
        # dQ2Y=-XFilter.RALF1FilterX(dQ3C+mDD,sz,sz*tSp,1,0)
        # dQ2X=(dQ2X+dQ2Y)
        # dQ2Y=dQ2X.copy()
        dQ5mx=np.zeros((sz,sz),np.float16)
        dQ5mn=np.zeros((sz,sz),np.float16)
        r4=np.zeros((3,sz*tSp),np.float16)
        for i in range(sz):
            r4[0]= (liix[i]).copy()
            r4[1]= (dQ2X[i]).copy()
            r4[2]= (dQ2Y[i]).copy()
            m=[[r4[j][l] for j in range(len(r4))] for l in range(len(r4[0]))]         
            m.sort(key=itemgetter(0))                  
            r4=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]         
        
            anum1=-1
            for j in range(sz): 
                jjj0=anum1
                anum0=max(0,min(jjj0,sz*tSp-1))
                while jjj0<sz*tSp-1 and round(r4[0][anum0])<j:
                    jjj0=jjj0+1
                    anum0=max(0,min(jjj0,sz*tSp-1))
                jjj1=anum0+1
                anum1=max(0,min(jjj1,sz*tSp))
                while jjj1<sz*tSp and round(r4[0][anum1])<j+1:
                    jjj1=jjj1+1
                    anum1=max(0,min(jjj1,sz*tSp))
                dQ5mx[i][j]=max(r4[1][anum0:anum1])
                dQ5mn[i][j]=min(r4[2][anum0:anum1])
                
        dQ5mx=dQ5mx.transpose()
        dQ5mn=dQ5mn.transpose()
        
        dQ5mx=dQ5mx-dQ5mx*(dQ5mx<0)
        dQ5mn=dQ5mn-dQ5mn*(dQ5mn>0)
        
        aMx=np.zeros(sz,np.float16)
        aMn=np.zeros(sz,np.float16)
        
        for i in range(sz):
            aMx[i]=max(dQ5mx[i])  
            aMn[i]=min(dQ5mn[i])          
            
        Nfl=int(len(arr_bx)/NChan)
        # for l in range(NChan):
        #     aMx[Nfl-NNew+Nfl*l:Nfl+Nfl*l]= savgol_filter(aMx[Nfl-NNew+Nfl*l:Nfl+Nfl*l], 15, 5)
        #     aMn[Nfl-NNew+Nfl*l:Nfl+Nfl*l]= savgol_filter(aMn[Nfl-NNew+Nfl*l:Nfl+Nfl*l], 15, 5)
        
        arr_bbbxxx=aMx + aMn  
            
        arr_bbbxxx=filterFourierQ(arr_bbbxxx,arr_bx,NNew,NChan)
        ann=sum(np.isnan(arr_bbbxxx))
        if ann==0: 
            arr_bbx.append(arr_bbbxxx)           
        hh=hh+1
    
    arr_bbx=np.asarray(arr_bbx,float).transpose()
    for l in range(NChan):
        for ii in range(NNew):  
            arr_b[Nfl-NNew+Nfl*l+ii]=(max(arr_bbx[Nfl-NNew+Nfl*l+ii])+min(arr_bbx[Nfl-NNew+Nfl*l+ii]))/2

    arr_b=filterFourierQ(arr_b,arr_b,NNew,NChan)+mn
    return arr_bbbxxx

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
    arr_bZ=[]
    for l in range(NChan):
        arr_bZ.append(arr_b[0+Nf*l:Nf-NNew+Nf*l])    
    arr_bZ=np.asarray(arr_bZ,float)
    D=np.std(arr_bZ)
    arr_b=np.asarray(arr_bb,np.float16)    
    NNew=int(NNew*1.1)   
    while 1==1: 
        hh=0
        ann=0
        arr_bbx=[]
        Nch=0
        Koef=np.zeros(Nhh,float)
        KoefA=np.zeros(Nhh,float)
        while hh<Nhh:
            arr_bbbxxx=RALF1Calculation(arr_b,Nf,NNew,NChan,D,Nhh)
            Nf_=int(NNew*1.8)
            NNew_=Nf_-NNew
            arr_bbbxxx_=np.zeros(Nf_*NChan,np.float16)
            for l in range(NChan):
                dd_=arr_bbbxxx[Nf-1+Nf*l:Nf-NNew+Nf*l:-1].copy()
                arr_bbbxxx_[0+Nf_*l:Nf_+Nf_*l]=(np.concatenate((dd_,np.ones(Nf_-len(dd_),float)*dd_[len(dd_)-1])))  
            if not (np.abs(np.mean(arr_bbbxxx_))==np.Inf or np.abs(np.mean(arr_bbbxxx_))==np.NaN):
                arr_bbbxxx_y=RALF1Calculation(arr_bbbxxx_,Nf_,NNew_,NChan,D,Nhh)
                
                arr_bbbxxx_yy=[]
                
                for l in range(NChan):
                    dd_=arr_bbbxxx_y[Nf_-1+Nf_*l:Nf_-NNew_+Nf_*l:-1].copy()
                    arr_bbbxxx_yy.append(dd_) 
                    if l==0:
                        mm1=arr_b[Nf-NNew-len(dd_):Nf-NNew].copy()
                        mm2=arr_bbbxxx_yy[l].copy()
                    else:
                        mm1=np.concatenate((mm1,arr_b[Nf-NNew-len(dd_):Nf-NNew]))
                        mm2=np.concatenate((mm2,arr_bbbxxx_yy[l]))
                        
                ann=sum(np.isnan(arr_bbbxxx))
                if ann==0: 
                    Koef[hh]=-np.std(mm1-mm2)
                    if np.mean(np.abs(mm1-np.mean(mm1)))>0:
                        KoefA[hh]=100*scp.spearmanr(mm1,mm2)[0]
                    else:
                        KoefA[hh]=0
                    arr_bbx.append(arr_bbbxxx)           
                    hh=hh+1
            else:
                hh=Nhh+2
        if hh<Nhh+2:
            arr_bbx=np.asarray(arr_bbx,np.float16).transpose()
            
            r2=np.zeros((2,Nhh),float)
            r2[0]= np.asarray(KoefA,float)
            r2[1]= np.asarray(range(Nhh),float)
            m=[[r2[j][l] for j in range(len(r2))] for l in range(len(r2[0]))]         
            m.sort(key=itemgetter(0))                  
            r2=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]  
            Nch=int(r2[1][Nhh-1])
            if np.isnan(KoefA[Nch]):
                KoefA[Nch]=0            
            if KoefA[Nch]>0:
                for l in range(NChan):
                    arr_b[Nf-NNew+Nf*l:Nf+Nf*l]=arr_bbx[Nf-NNew+Nf*l:Nf+Nf*l,Nch].copy()    
                arr_b=filterFourierQ(arr_b,arr_b,NNew,NChan)
                return arr_b

if __name__ == '__main__':
    RALf1FiltrQ(sys.argv)
