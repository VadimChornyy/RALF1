import numpy as np
from operator import itemgetter
import time as tm
import RALF1FilterX as XFilter
import sys
import lfib1340 
from scipy import stats as scp
import win32api,win32process,win32con
#from scipy.signal import savgol_filter
           
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
        farxxx2=farxxx.copy()
        arxr[Nfl-Nnl+Nfl*l:Nfl+Nfl*l]=np.fft.ifft(farxxx2).real[0:Nnl] 
        arxr[0+Nfl*l:Nfl-Nnl+Nfl*l]=arb[0+Nfl*l:Nfl-Nnl+Nfl*l].copy() 
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
    Koe=1e-4 
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
    aa=RandomQ(sz) 
    liiC=np.concatenate((aa,aa))
    
    aa=RandomQ(sz) 
    liiB=np.concatenate((aa,aa))   
    aa=RandomQ(sz) 
    liiD=np.concatenate((aa,aa))
    aa=RandomQ(sz) 
    liiE=np.concatenate((aa,aa))
    
    r4=np.zeros(sz,np.float16)
    for l in range(NChan):            
        r4[Nf-NNew+Nf*l:Nf+Nf*l]=RandomQ(NNew)/NNew 
        r4[Nf-NNew+Nf*l:Nf+Nf*l]=D*((r4[Nf-NNew+Nf*l:Nf+Nf*l]/np.std(r4[Nf-NNew+Nf*l:Nf+Nf*l]))/2+Koe*10) 

    r2=np.asarray(arr_b,np.float16)
    for l in range(NChan):                
        r2[Nf-NNew+Nf*l:Nf+Nf*l]=mn

    R4=np.asarray(r4,np.float16)       
        
    K=NNew/(Nf+1)/NChan

    liix=np.zeros((sz,sz),int) 
    dQ3_0=np.zeros((sz,sz),np.float16)
    mDD=np.zeros((sz,sz),np.float16)  

    MM=int(np.ceil(sz/300)+1)

    Ndel=MM#int(np.ceil(np.sqrt(sz)))
    NCh=int(np.ceil(sz/Ndel)) 
    Ndel0=MM
    Nzz=Nhh#*10#int(np.ceil(np.sqrt(Ndel)))
    NCh0=int(np.ceil(sz/Ndel0))    
    dQ4=np.zeros((NCh,NCh0),np.float16)
    mDD4=np.zeros((NCh,NCh0),np.float16)      
    while hh<Nhh:              
        for i in range(sz):    
            r1=liiB[int(liiD[i+hh]):sz+int(liiD[i+hh])]                                     
            liix[i]=r1.copy()                    
            dQ3_0[i]=r2[r1].copy()   
            for l in range(NChan):
                bb=np.asarray(liiC[np.asarray(np.arange(l+int(liiE[i+hh]),l+sz+int(liiE[i+hh]),sz/NNew),int)]*K,int)
                R4[Nf-len(bb)+Nf*l:Nf+Nf*l]=r4[Nf-len(bb)+Nf*l+bb].copy()                    
            mDD[i]=R4[r1].copy()     
            tm.sleep(0.002) 
        
        dQ3_0=dQ3_0-mn   

        ##########################################
        w=1
        ss4=0
        while w>0:
            try:   
                dQ3=dQ3_0.copy()      
                NumFri0=RandomQ(sz)
                NumFri0_=RandomQ(sz)                 
                NumFri0=np.concatenate((NumFri0, NumFri0, NumFri0))                  
                NumFri0_=np.concatenate((NumFri0_, NumFri0_, NumFri0_))  
                r5=RandomQ(sz) 
                r5=D*((r5/np.std(r5))/2+Koe*10) 
                r5=np.concatenate((r5, r5))
                zz=0
                while zz<Nzz:
                    dQ3mx=np.zeros((sz,sz),np.float16)-np.Inf
                    dQ3mn=np.zeros((sz,sz),np.float16)+np.Inf         
                    ss4=ss4+1
                    ss4_=ss4-int(ss4/sz)*sz
                    NumFri=NumFri0[ss4_:].copy()
                    NumFri_=NumFri0_[ss4_:].copy()                        
                    for kk in range(Ndel):
                        ii=int(kk*NCh)
                        for k in range(Ndel0):
                            i=int(k*NCh0) 
                            dQ4=np.zeros((NCh,NCh0),float)
                            mDD4=np.zeros((NCh,NCh0),float)                                    
                            for ll in range(NCh0):
                                dQ4[:,ll]=(dQ3[NumFri[ii:ii+NCh],NumFri_[i+ll]])*1.
                                mDD4[:,ll]=(r5[ll+k+kk:ll+k+kk+NCh]*(1-(mDD[NumFri[ii:ii+NCh],NumFri_[i+ll]]<D*Koe)))*1.
                                #mDD4[:,ll]=(mDD[NumFri[ii:ii+NCh],NumFri_[i+ll]])
                            dQ4_=np.mean(dQ4*(mDD4<D*Koe))
                            dQ4=dQ4-dQ4_
                            dQ4_A= np.asarray(XFilter.RALF1FilterX(  dQ4*(1-(dQ4<0))+mDD4,len(dQ4),len(dQ4[0]),1,0),np.float16)
                            dQ4_B=-np.asarray(XFilter.RALF1FilterX( -dQ4*(1-(dQ4>0))+mDD4,len(dQ4),len(dQ4[0]),1,0),np.float16)
                            dQ4A=dQ4_A+dQ4_B+dQ4_
                            dQ4_A=dQ4A.copy()
                            dQ4_B=dQ4A.copy()
                            for ll in range(NCh0):
                                dQ3mx[NumFri[ii:ii+NCh],NumFri_[i+ll]]=np.maximum(dQ3mx[NumFri[ii:ii+NCh],NumFri_[i+ll]],dQ4_A[:,ll])
                                dQ3mn[NumFri[ii:ii+NCh],NumFri_[i+ll]]=np.minimum(dQ3mn[NumFri[ii:ii+NCh],NumFri_[i+ll]],dQ4_B[:,ll])

                    zz=zz+1
                    if zz==1:
                        AsrXMx=dQ3mx.copy()
                        AsrXMn=dQ3mn.copy()
                    else:
                        AsrXMx=(AsrXMx*(zz-1)+(dQ3mx))/zz
                        AsrXMn=(AsrXMn*(zz-1)+(dQ3mn))/zz
                    
                Asr=AsrXMx+AsrXMn
                Asr=Asr*np.std(mDD*(1-(mDD<D*Koe)))/np.std(Asr*((mDD<D*Koe)*1))
                Asr_=np.mean(Asr*(mDD<D*Koe))
                Asr=Asr-Asr_
                Asr=dQ3_0*(mDD<D*Koe)+Asr*(np.asarray(1,np.float16)-(mDD<D*Koe))
               
                dQ3=np.asarray( XFilter.RALF1FilterX(Asr*(1-(Asr<0))+mDD,sz,sz,1,0)-
                     XFilter.RALF1FilterX(-Asr*(1-(Asr>0))+mDD,sz,sz,1,0),np.float16) #+Asr_ 
                #dQ3=dQ3+Asr_
                
                # Asr=(mDD+mDD.transpose()+D*Koe*10)*(np.asarray(1,np.float16)-(mDD<D*Koe))    
                # ddd=np.std(np.asarray(mDD*(1-(mDD<D*Koe)),float))/np.std(np.asarray(Asr*(1-(mDD<D*Koe)),float))
                # mDD=np.asarray(Asr*ddd,np.float16)
                
                # dQ3=dQ3*np.std(mDD*(1-(mDD<D*Koe)))/np.std(dQ3*((mDD<D*Koe)*1))
                # dQ3=dQ3_0*(mDD<D*Koe)+dQ3*(np.asarray(1,np.float16)-(mDD<D*Koe))
                                          
                w=w-1
            except:
                w=w
            ##########################################
        #dQ3_0[:][liix]=dQ3
        for i in  range(sz):
            dQ3_0[i][liix[i]]=dQ3[i].copy()
                   
        aMx=np.max(dQ3_0,0)
        aMn=np.min(dQ3_0,0)          
        
        # Nfl=int(len(arr_bx)/NChan)
        # for l in range(NChan):      
        #     aMx[0+Nfl*l:Nfl+Nfl*l]= savgol_filter(aMx[0+Nfl*l:Nfl+Nfl*l], 11, 5)
        #     aMn[0+Nfl*l:Nfl+Nfl*l]= savgol_filter(aMn[0+Nfl*l:Nfl+Nfl*l], 11, 5)
        
        ann=sum(np.isnan(aMx + aMn))
        if ann==0: 
            if hh==0:                
                AMX=aMx.copy()
                AMN=aMn.copy()                
            else:
                AMX=np.maximum(AMX,aMx)
                AMN=np.minimum(AMN,aMn)
                
            arr_bbbxxx=AMX+AMN
            arr_bbbxxx=filterFourierQ(arr_bbbxxx,arr_b,NNew,NChan)
            if hh==0:
                arr_bbx=arr_bbbxxx.copy()            
            else:               
                arr_bbx=(arr_bbx*hh+arr_bbbxxx)/(hh+1)                  
            hh=hh+1

    return arr_bbx+mn

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
            if (sum(np.abs(arr_bbbxxx_)==np.Inf)==0 and sum(np.isnan(arr_bbbxxx_))==0):
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
                
                ann=(sum(np.abs(mm1)==np.Inf)>0 + sum(np.isnan(mm1))>0+
                     sum(np.abs(mm2)==np.Inf)>0 + sum(np.isnan(mm2))>0)
                
                if ann==0 and len(mm1)>1 and len(mm1)==len(mm2):
                    mm1=mm1-sum(mm1)/len(mm1)
                    mm2=mm2-sum(mm2)/len(mm1)
               
                    if np.std(mm1)>0 and np.std(mm2)>0:
                        KoefA[hh]=100*scp.spearmanr(mm1,mm2)[0]
                        mm1=mm1*np.std(mm2)/np.std(mm1)                       
                        Koef[hh]=-np.std(mm1-mm2)
                        arr_bbx.append(arr_bbbxxx)           
                        hh=hh+1
            else:
                hh=Nhh+2
        if hh<Nhh+2:            
            arr_bbx=np.asarray(arr_bbx,np.float16).transpose()  
            Koef=Koef/np.std(Koef)+KoefA/np.std(KoefA)
            r2=np.zeros((2,Nhh),float)
            r2[0]= np.asarray(Koef,float)
            r2[1]= np.asarray(range(Nhh),float)
            m=[[r2[j][l] for j in range(len(r2))] for l in range(len(r2[0]))]         
            m.sort(key=itemgetter(0))                  
            r2=[[m[j][l] for j in range(len(m))] for l in range(len(m[0]))]  
            Nch=int(r2[1][Nhh-1])
            print(KoefA)
            if np.isnan(KoefA[Nch]):
                KoefA[Nch]=0            
            if KoefA[Nch]>20:
                for l in range(NChan):
                    arr_b[Nf-NNew+Nf*l:Nf+Nf*l]=arr_bbx[Nf-NNew+Nf*l:Nf+Nf*l,Nch].copy()    
                arr_b=filterFourierQ(arr_b,arr_b,NNew,NChan)
                return arr_b

if __name__ == '__main__':
    RALf1FiltrQ(sys.argv)
