import numpy as np
from operator import itemgetter
import time as tm
import sys
import lfib1340 
#from scipy.signal import savgol_filter
  
def RALF1FilterQ(dQ2):
    Np=len(dQ2)  
    SdQ=np.mean(dQ2,0)
    stddQ=np.std(np.asarray(dQ2,float))
        
    sSdQ=np.std(np.asarray(SdQ,float))
    for i in range(Np):
        SdQj_ = np.std(np.asarray(dQ2[i] - SdQ,float))
        SdQj__ = np.std(np.asarray(dQ2[i],float))              
        if SdQj__ >0 and sSdQ>0:
            #dd = dd * (1+sSdQ / SdQj__) / 2 
            dQ2[i]=np.asarray(dQ2[i] +SdQ * (SdQj_ / sSdQ - 1),np.float16)
        else:
            dQ2[i]=dQ2[i]*0  
    stddQ_=np.std(np.asarray(dQ2,float))
    if stddQ_>0:
        dQ2=dQ2*stddQ/stddQ_
    return dQ2

def RandomV(Nfx):
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

def filterFourierV(arxx,arb,NNew,NChan):  
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
    srfarx=.63*np.mean(farx[1:])/2
    arxr=np.zeros(Nfl*NChan,float)   
    for l in range(NChan):       
        farxx=np.fft.fft(arxx[Nfl-Nnl+Nfl*l:Nfl+Nfl*l])    
        mfarxx=abs(farxx) 
        mfarxx[0]=1e-32
        srmfarxx=.63*np.mean(mfarxx[1:])
        farxxx=np.zeros(Nnl,complex)     
        for j in range(Nnl):
            if mfarxx[j]>srmfarxx:# and farx[j]>srfarx:
                farxxx[j]=farxx[j]/mfarxx[j]*farx[j]            
            else:
                farxxx[j]=0        
        arxr[Nfl-Nnl+Nfl*l:Nfl+Nfl*l]=np.fft.ifft(farxxx).real  
        arxr[0+Nfl*l:Nfl-Nnl+Nfl*l]=arb[0+Nfl*l:Nfl-Nnl+Nfl*l].copy() 

    return arxr

def RALf1FiltrV(args):
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
    mn=np.mean(arr_bZ)
    arr_b=np.asarray(arr_bb,np.float16)    
        
    hh=0
    ann=0
     
    NNew=int(NNew*1.1)
    arr_bbx=[]
    while hh<Nhh:
        liiB=np.zeros(2*Nf*NChan,int)
        aa=RandomV(Nf*NChan) 
        liiB[0:Nf*NChan]=aa
        liiB[Nf*NChan:2*Nf*NChan]=aa        
        
        liiC=np.zeros(2*(Nf+1)*NChan,int)
        aa=RandomV((Nf+1)*NChan) 
        liiC[0:(Nf+1)*NChan]=aa
        liiC[(Nf+1)*NChan:2*(Nf+1)*NChan]=aa   
        
        liiD=RandomV(Nf*NChan)
        liiE=RandomV(Nf*NChan)
        
        r4=np.zeros(Nf*NChan,float)
        for l in range(NChan):            
            r4[Nf-NNew+Nf*l:Nf+Nf*l]=RandomV(NNew)/NNew 
            r4[Nf-NNew+Nf*l:Nf+Nf*l]=D*(r4[Nf-NNew+Nf*l:Nf+Nf*l]/np.std(r4[Nf-NNew+Nf*l:Nf+Nf*l])/2+1e-6) 
                            
        r2=np.asarray(arr_b,np.float16)
        for l in range(NChan):                
            r2[Nf-NNew+Nf*l:Nf+Nf*l]=mn
        r2=r2-mn
        R4=np.asarray(r4,np.float16)
        K=NNew/(Nf+1)/NChan
        sz=int(NChan*Nf)
        liix=[[] for kk in range(Nf*NChan)]  
        dQ3=np.zeros((sz,sz),np.float16)
        mDD=np.zeros((sz,sz),np.float16)        
        for i in range(sz):    
            r1=liiB[int(liiD[i]):sz+int(liiD[i])]                                     
            liix[i].append(np.asarray(r1,int))                     
            dQ3[i]=r2[r1]
            for l in range(NChan):
                bb=np.asarray(liiC[np.asarray(np.arange(l+int(liiE[i]),l+sz+int(liiE[i]),sz/NNew),int)]*K,int)
                R4[Nf-NNew+Nf*l:Nf+Nf*l]=r4[Nf-NNew+Nf*l+bb[0:NNew]]                 
            mDD[i]=R4[r1]  
            tm.sleep(0.002)
             
        line=1
        while line==1: 
            tm.sleep(0.2)
            #anamef="fralf.tmp"
            try:
            #     fo = open(anamef, "r")
            #     line = int(fo.readline()) 
            # except:
            #     fo = open(anamef, "w")
            #     fo.write(str(0)+'\n') 
                line=0                       
            #fo.close()
            #if line==0:
                # fo = open(anamef, "w")
                # fo.write(str(1)+'\n')
                # fo.close()               
             
                Ndel=1#int(np.ceil(np.sqrt(sz)))
                NCh=int(np.ceil(sz/Ndel)) 
                Ndel0=1#2*Ndel
                NCh0=int(np.ceil(sz/Ndel0))                    
                dQ3mx=np.zeros((sz,sz),np.float16)-np.Inf
                dQ3mn=np.zeros((sz,sz),np.float16)+np.Inf
                dQ4=np.zeros((NCh,NCh0),np.float16)
                mDD4=np.zeros((NCh,NCh0),np.float16)
                NumFri=RandomV(sz)
                NumFri_=RandomV(sz)                 
                NumFri=np.concatenate((NumFri, NumFri, NumFri))                  
                NumFri_=np.concatenate((NumFri_, NumFri_, NumFri_))  
                zz=Nhh-1#int(np.ceil(np.sqrt(Ndel)))
                while zz>=0:
                    for kk in range(Ndel):
                        ii=int(kk*NCh)
                        for k in range(Ndel0):
                            i=int(k*NCh0) 
                            dQ4=[]
                            mDD4=[]
                            for ll in range(NCh0):
                                dQ4.append(dQ3[NumFri[zz+ii+0:zz+ii+NCh],NumFri_[i+ll]])
                                mDD4.append(mDD[NumFri[zz+ii+0:zz+ii+NCh],NumFri_[i+ll]])
                            dQ4=np.asarray(dQ4,np.float16).transpose()
                            mDD4=np.asarray(mDD4,np.float16).transpose()
                            dQ4mn=np.mean(dQ4*(1-(np.abs(mDD4)<D*1e-6)))
                            dQ4=dQ4-dQ4mn                 
                            dQ4=(RALF1FilterQ(  dQ4-dQ4*(dQ4<0) +mDD4)-
                                 RALF1FilterQ(-(dQ4-dQ4*(dQ4>0))+mDD4))
                            dQ4=dQ4+dQ4mn
                            for ll in range(NCh0):
                                dQ3mx[NumFri[zz+ii+0:zz+ii+NCh],NumFri_[i+ll]]=np.maximum(dQ3mx[NumFri[zz+ii+0:zz+ii+NCh],NumFri_[i+ll]],dQ4[:,ll])
                                dQ3mn[NumFri[zz+ii+0:zz+ii+NCh],NumFri_[i+ll]]=np.minimum(dQ3mn[NumFri[zz+ii+0:zz+ii+NCh],NumFri_[i+ll]],dQ4[:,ll])
                    zz=zz-1   
                dQ3=(dQ3mx+dQ3mn)/2
                del(dQ3mx)
                del(dQ3mn)
                # fo = open(anamef, "w")
                # fo.write(str(0)+'\n')
                # fo.close()   
            except:
                line=0
        
        del(mDD)       
        for i in range(sz):
            r1=np.asarray(liix[i],int)
            dQ3[i][r1]=dQ3[i].copy()     
        aMx=np.max(dQ3,0)
        aMn=np.min(dQ3,0)        
        del(liix)
        del(dQ3)
        
        # for l in range(NChan):                
        #     aMx[0+Nf*l:Nf+Nf*l]= savgol_filter(aMx[0+Nf*l:Nf+Nf*l], 11, 5)
        #     aMn[0+Nf*l:Nf+Nf*l]= savgol_filter(aMn[0+Nf*l:Nf+Nf*l], 11, 5)
        arr_bbbxxx=aMx + aMn  
        
        arr_bbbxxx=filterFourierV(arr_bbbxxx,arr_b,NNew,NChan)
        
        ann=sum(np.isnan(arr_bbbxxx))
        if ann==0: 
            arr_bbx.append(arr_bbbxxx)           
            hh=hh+1
    
    arr_bbx=np.asarray(arr_bbx,np.float16).transpose()
    for l in range(NChan): 
        for ii in range(NNew):  
            arr_b[ii+Nf-NNew+Nf*l]=(max(arr_bbx[ii+Nf-NNew+Nf*l])+min(arr_bbx[ii+Nf-NNew+Nf*l]))/2        
    #arr_b=filterFourierV(arr_b,arr_b,NNew,NChan)
         
    return arr_b+mn

if __name__ == '__main__':
    RALf1FiltrV(sys.argv)