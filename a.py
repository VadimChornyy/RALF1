import hickle as hkl
import numpy as np
from RALf1FiltrVID import filterFourierQ
from scipy.signal import savgol_filter
import pylab as plt
from PIL import Image  
import cv2 as cv
import os
def fig2dsk ( fig ):
    fig.savefig(wrkdir +'_tm_dynamic.png',dpi=150,transparent=False,bbox_inches = 'tight')
    frame=Image.open(wrkdir +'_tm_dynamic.png')
    return frame
lSrez=0.99
arrrxx=hkl.load("ralfrez.rlf2")
# i=0 
# arrrxx_=[]
# for ii in arrrxx:
#     if not i==4:
#         arrrxx_.append(ii)
#     i=i+1
# hkl.dump(arrrxx_,"ralfrez.rlf2")
try:
    ZDat=np.asarray(arrrxx,float)#.transpose()
    anI=len(ZDat)
    # for i in range(anI):  
    #     ZDat[i]= savgol_filter(ZDat[i],65, 5)
    for ii in range(3):       
        ar0x=np.mean(ZDat,axis=0)
        ar0x_=.4*(np.mean(abs((ZDat)-(ar0x)),axis=0))
            
        lnn=len(ZDat[0])
        NNew=int(.35*lnn)
        for i in range(anI):    
            for j in range(lnn):    
                if not abs(ZDat[i,j]-(ar0x[j]))<=ar0x_[j]:     
                    if ZDat[i,j]<((ar0x[j])-ar0x_[j]):            
                        ZDat[i,j]=(ar0x[j])-ar0x_[j]
                    else:
                        if ZDat[i,j]>((ar0x[j])+ar0x_[j]):
                            ZDat[i,j]=(ar0x[j])+ar0x_[j]
        
        P=np.zeros(3,float)
        for i in range(anI):
            dd=ZDat[i][lnn-NNew:].copy()                         
            x=ar0x[lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))].copy()
            ZDat[i][lnn-NNew:]=filterFourierQ(ZDat[i],(ar0x),NNew,1)[lnn-NNew:]
            # ZDat[i]= savgol_filter(ZDat[i],65, 5)
            P[0:2]=np.polyfit(x,ZDat[i][lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))],1)
            if not P[0]>0:
                P[0:2]=np.polyfit(dd,ZDat[i][lnn-NNew:],1)
            ZDat[i][lnn-NNew:]=(ZDat[i][lnn-NNew:]-P[1])/P[0]
            ZDat[i]= savgol_filter(ZDat[i], 14, 5)  
                                
    bbbbb=ZDat[:,:].transpose().copy()
    aaaaa=np.mean(bbbbb.transpose(),axis=0)
        
    %varexp --plot bbbbb 
    len(arrrxx)
except:
    len(arrrxx)
len(arrrxx)
aname="MANA"

fourcc = 0
ImApp=[]
gray_sz1=0

dNIt=4
Ngroup=3
Lo=1
Nproc=Ngroup*(os.cpu_count())
wrkdir = r""
[hhhao,Arr_AAA]=(hkl.load(wrkdir + aname+".rlf1"))
NIter=100

for hhhai in range(hhhao):  
    hhha=hhhai+1
    for iGr in range(Ngroup):  
        ZDat=(Arr_AAA[iGr*NIter*int(Nproc/Ngroup)+max(0,(hhhai+1)-dNIt)*int(Nproc/Ngroup):iGr*NIter*int(Nproc/Ngroup)+(hhhai+1)*int(Nproc/Ngroup)]).copy()
        if iGr==0:
            xxxx=ZDat.copy()
        else:
            xxxx=np.concatenate((xxxx, ZDat))
    ZDat=xxxx.copy()   
    anI=len(ZDat)
    # for i in range(anI):  
    #     ZDat[i]= savgol_filter(ZDat[i], 65, 5)
    for ii in range(3):       
        ar0x=np.mean(ZDat,axis=0)
        ar0x_=0.4*(np.mean(abs((ZDat)-(ar0x)),axis=0))
            
        lnn=len(ZDat[0])
        NNew=int(.35*lnn)
        for i in range(anI):    
            for j in range(lnn):    
                if not abs(ZDat[i,j]-(ar0x[j]))<=ar0x_[j]:     
                    if ZDat[i,j]<((ar0x[j])-ar0x_[j]):            
                        ZDat[i,j]=(ar0x[j])-ar0x_[j]
                    else:
                        if ZDat[i,j]>((ar0x[j])+ar0x_[j]):
                            ZDat[i,j]=(ar0x[j])+ar0x_[j]

        P=np.zeros(3,float)
        for i in range(anI):
            dd=ZDat[i][lnn-NNew:].copy()                         
            x=ar0x[lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))].copy()
            ZDat[i][lnn-NNew:]=filterFourierQ(ZDat[i],(ar0x),NNew,1)[lnn-NNew:]
            # ZDat[i]= savgol_filter(ZDat[i], 65, 5)
            P[0:2]=np.polyfit(x,ZDat[i][lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))],1)
            if not P[0]>0:
                P[0:2]=np.polyfit(dd,ZDat[i][lnn-NNew:],1)
            ZDat[i][lnn-NNew:]=(ZDat[i][lnn-NNew:]-P[1])/P[0]    
            ZDat[i]= savgol_filter(ZDat[i], 14, 5)                  
    #ZDat=Arr_AAA[iGr*NIter*int(Nproc/Ngroup)+max(0,(hhh+1)-dNIt)*int(Nproc/Ngroup):iGr*NIter*int(Nproc/Ngroup)+(hhh+1)*int(Nproc/Ngroup)].copy()
    if hhhai==0:
        bbbbb=np.mean(ZDat[:,:],axis=0).transpose().copy()
    else:
        bbbbb=(bbbbb*hhhai+np.mean(ZDat[:,:],axis=0).transpose())/(hhhai+1)
    aaaaa=(np.amax(bbbbb.transpose(),axis=0)+np.amin(bbbbb.transpose(),axis=0))/2
    fig = plt.figure()
    axes = fig.add_axes([0.2, 0.2, 0.8, 0.8])    
    axes.plot(bbbbb,'o')
    frame=fig2dsk(fig) 
    cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)   
    if gray_sz1==0:       
        gray_sz1=len(cimg[0])
        gray_sz2=len(cimg)
    else:
        gray_sz1=min(gray_sz1,len(cimg[0]))
        gray_sz2=min(gray_sz2,len(cimg))
    ImApp.append(frame)
    if fourcc==0:
        aDur=3
        fourcc = cv.VideoWriter_fourcc(*'MP4V')
    out = cv.VideoWriter(wrkdir + aname+'_tm_'+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
    for icl in range(len(ImApp)):
        cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
        out.write(cimgx[0:gray_sz2,0:gray_sz1,:]) 
    out.release()
    del(out)
    plt.show()

