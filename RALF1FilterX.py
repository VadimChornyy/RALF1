import numpy as np
def meani(bb):
    aa=0
    sz=len(bb)
    for i in range(sz):
        aa=np.float16((aa*i+bb[i])/(i+1))  
    return aa    
    
def stdeviat(bb):
    aa=0
    sz=len(bb)
    bb=bb-meani(bb)
    for i in range(sz):
        aa=np.float16((aa*i*i+bb[i]*bb[i])/(i+1)/(i+1))
    aa=np.sqrt(aa)
    return aa

def RALF1FilterX(dQ1,Np,Nf,key,key2):
    dQ2=np.asarray(dQ1,np.float16)
    if key>0: 
        SdQj=np.ones(Np,np.float16)   
        dSdQj=np.zeros(Np,np.float16)
        znakSdQj=np.ones(Np,np.float16) 
        if key2>0:
            for i in range(Np):
                znakSdQj[i]=meani(dQ2[i])        
            znakSdQj=((znakSdQj<0)-0.5)*2*(-1) 
 
        for i in range(Np): 
            #SdQj[i]=stdeviat(dQ2[i])
            dQ2[i] = dQ2[i] * znakSdQj[i]            
            if SdQj[i]==0:
                dQ2[i]=np.zeros(Nf,np.float16)
            else:
                dQ2[i]=dQ2[i] / SdQj[i]
        
        dQ2_=dQ2.transpose()
        SdQ=np.zeros(Nf,np.float16)         
        for j in range(Nf):
            SdQ[j]=meani(dQ2_[j])  
        sSdQ=stdeviat(np.asarray(SdQ,np.float16))
        for i in range(Np):
            SdQj_ = stdeviat(np.asarray(dQ2[i] - SdQ,np.float16))
            SdQj__ = stdeviat(np.asarray(dQ2[i],np.float16))            
            if SdQj__ >0 and sSdQ>0:
                dQ2[i] = dQ2[i] +SdQ * (SdQj_ / sSdQ - 1)
                dQ2[i] = dQ2[i] *znakSdQj[i] * SdQj[i]#* (1+sSdQ / SdQj__) / 2 
            else:
                dQ2[i]=np.zeros(Nf,np.float16)        
            dQ2[i]=dQ2[i]+dSdQj[i]
        #dQ2=dQ2-np.mean(dQ2)+mQ
    return np.asarray(dQ2,np.float16)
