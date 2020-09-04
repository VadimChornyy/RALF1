import numpy as np
from operator import itemgetter
import time as tm
import RALF1FilterX as XFilter
import sys
import lfib1340 
#from scipy.signal import savgol_filter
from scipy import stats as scp
import scipy.interpolate as scpyi 

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

def RandomQ(Nfx):
    KK=3e6
    liiX=np.zeros(Nfx,float)
    for ii in range(3): 
        z=np.random.randint(Nfx)/KK           
        atim0=tm.time()        
        tm.sleep(z) 
        atim=tm.time()     
        dd=int((atim-atim0-z)*KK)
        zz=np.asarray(range(Nfx),float)/KK
        lfib1340.LFib1340(dd).shuffle(zz)   
        liiX=liiX+zz
            
    r2=np.zeros((2,Nfx),float)
    r2[0]= np.asarray(liiX[0:Nfx],float)
    r2[1]= np.asarray(range(Nfx),int)
    m=[[r2[j][l] for j in range(len(r2))] for l in range(len(r2[0]))]         
    m.sort(key=itemgetter(0))                  
    r2=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]  
    liiXX=np.asarray(r2[1],int)
    return liiXX
  
def RALF1Calculation(arr_bx,Nf,NNew,NChan,D):
    Koe=1e-6 
    tSp=1
    arr_bZ=[]
    arr_b=np.asarray(arr_bx,float)
    #arr_b[0]=0
    for l in range(NChan):
        #arr_b[l]=arr_bx[l]-arr_bx[l-1]        
        arr_bZ.append(arr_b[0+Nf*l:Nf-NNew+Nf*l])    
    arr_bZ=np.asarray(arr_bZ,float)
    mn=np.mean(arr_bZ)
    
    aa=RandomQ(Nf*NChan) 
    liiB=np.concatenate((aa,aa))        
    aa=RandomQ((Nf+1)*NChan) 
    liiC=np.concatenate((aa,aa))         
    liiD=RandomQ(Nf*NChan)
    liiE=RandomQ(Nf*NChan)
    
    r4=np.zeros(Nf*NChan,float)
    for l in range(NChan):            
        r4[Nf-NNew+Nf*l:Nf+Nf*l]=RandomQ(NNew)/NNew 
        r4[Nf-NNew+Nf*l:Nf+Nf*l]=D*((r4[Nf-NNew+Nf*l:Nf+Nf*l]/np.std(r4[Nf-NNew+Nf*l:Nf+Nf*l]))/2+Koe*10) 
                        
    r2=np.asarray(arr_b,np.float16)
    for l in range(NChan):                
        r2[Nf-NNew+Nf*l:Nf+Nf*l]=mn
    r2=r2-mn
    R4=np.asarray(r4,np.float16)
    K=NNew/(Nf+1)/NChan
    sz=int(NChan*Nf)
    liix=[[] for kk in range(sz)]  
    
    w=1
    while w>0:
        try:
            dQ3=np.zeros((sz,sz*tSp),np.float16)
            mDD=np.zeros((sz,sz*tSp),np.float16)
            w=0
        except:
            w=1
    
    for i in range(sz):    
        r1=liiB[int(liiD[i]):sz+int(liiD[i])]                                     
        ge=scpyi.interp1d(np.asarray(range(sz),float),r1)                              
        liix[i]=np.asarray(ge(np.linspace(0,sz-1,sz*tSp)),int)
            
        ge=scpyi.interp1d(np.asarray(range(sz),float),r2) 
        dQ3[i]=np.asarray(ge(liix[i]),np.float16)
        
        for l in range(NChan):
            bb=np.asarray(liiC[np.asarray(np.arange(l+int(liiE[i]),l+sz+int(liiE[i]),sz/NNew),int)]*K,int)
            #bb_=np.asarray(liiC[np.asarray(np.arange(l+sz+int(liiE[i]),l+int(liiE[i]),-sz/NNew),int)]*K,int)
            R4[Nf-len(bb)+Nf*l:Nf+Nf*l]=(r4[Nf-len(bb)+Nf*l+bb])#+
                                      #r4[Nf-NNew+Nf*l+bb_[0:NNew]])/np.sqrt(2) 
        
        ge=scpyi.interp1d(np.asarray(range(sz),float),R4)
        mDD[i]=np.asarray(ge(liix[i]),np.float16)
        tm.sleep(0.002)
    w=1
    while w>0:
        try: 
            dQ3=( XFilter.RALF1FilterX(  dQ3*(1-(dQ3<0))+mDD,len(dQ3),len(dQ3[0]),1,0)-                    
                  XFilter.RALF1FilterX( -dQ3*(1-(dQ3>0))+mDD,len(dQ3),len(dQ3[0]),1,0))            
            w=0
        except:
            w=1
       
    dQ5mx=np.zeros((sz,sz),float)
    dQ5mn=np.zeros((sz,sz),float)
    r4=np.zeros((2,sz*tSp),float)
    for i in range(sz):
        r4[0]= np.asarray(liix[i],int).copy()
        r4[1]= (dQ3[i]).copy()
        m=[[r4[j][l] for j in range(len(r4))] for l in range(len(r4[0]))]         
        m.sort(key=itemgetter(0))                  
        r4=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]         

        anum1=-1;
        for jj in range(sz): 
            jjj0=anum1;
            anum0=max(0,min(jjj0,sz*tSp-1));
            while jjj0<sz*tSp-1 and int(r4[0][anum0])<jj:
                jjj0=jjj0+1;
                anum0=max(0,min(jjj0,sz*tSp-1));
            jjj1=anum0+1;
            anum1=max(0,min(jjj1,sz*tSp));
            while jjj1<sz*tSp and int(r4[0][anum1])<jj+1:
                jjj1=jjj1+1;
                anum1=max(0,min(jjj1,sz*tSp));
            dQ5mx[i][jj]=max(r4[1][anum0:anum1])
            dQ5mn[i][jj]=min(r4[1][anum0:anum1])
            
    dQ5mx=dQ5mx.transpose()
    dQ5mn=dQ5mn.transpose()
    
    aMx=np.zeros(sz,float)
    aMn=np.zeros(sz,float)
    
    for i in range(sz):
        aMx[i]=max(dQ5mx[i])  
        aMn[i]=min(dQ5mn[i])          
    
    arr_bbbxxx=(aMx + aMn)/2 
    arr_bbbxxx=filterFourierQ(arr_bbbxxx,arr_bx,NNew,NChan)+mn
    return arr_bbbxxx

def RALf1FiltrQ(args):
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
            arr_bbbxxx=RALF1Calculation(arr_b,Nf,NNew,NChan,D)
            Nf_=int(NNew*1.8)
            NNew_=Nf_-NNew
            arr_bbbxxx_=np.zeros(Nf_*NChan,np.float16)
            for l in range(NChan):
                dd_=arr_bbbxxx[Nf-1+Nf*l:Nf-NNew+Nf*l:-1].copy()
                arr_bbbxxx_[0+Nf_*l:Nf_+Nf_*l]=(np.concatenate((dd_,np.ones(Nf_-len(dd_),float)*dd_[len(dd_)-1])))  
            if not (np.abs(np.mean(arr_bbbxxx_))==np.Inf or np.abs(np.mean(arr_bbbxxx_))==np.NaN):
                arr_bbbxxx_y=RALF1Calculation(arr_bbbxxx_,Nf_,NNew_,NChan,D)
                
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
            if KoefA[Nch]>10:
                for l in range(NChan):
                    arr_b[Nf-NNew+Nf*l:Nf+Nf*l]=arr_bbx[Nf-NNew+Nf*l:Nf+Nf*l,Nch].copy()    
                arr_b=filterFourierQ(arr_b,arr_b,NNew,NChan)
                return arr_b

if __name__ == '__main__':
    RALf1FiltrQ(sys.argv)