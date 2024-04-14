import hickle as hkl
import numpy as np
from RALf1FiltrVID import filterFourierQ
from scipy.signal import savgol_filter

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
    for ii in range(3):       
        ar0x=np.median(ZDat,axis=0)
        ar0x_=.4*(np.median(abs((ZDat)-(ar0x)),axis=0))
            
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
        for i in range(anI):  
            ZDat[i]= savgol_filter(ZDat[i], 14, 5)
        P=np.zeros(3,float)
        for i in range(anI):
            dd=ZDat[i][lnn-NNew:].copy()                         
            x=ar0x[lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))].copy()
            ZDat[i][lnn-NNew:]=filterFourierQ(ZDat[i],(ar0x),NNew,1)[lnn-NNew:]
            P[0:2]=np.polyfit(x,ZDat[i][lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))],1)
            if not P[0]>0:
                P[0:2]=np.polyfit(dd,ZDat[i][lnn-NNew:],1)
            ZDat[i][lnn-NNew:]=(ZDat[i][lnn-NNew:]-P[1])/P[0]                      
    bbbbb=ZDat[:,:].transpose().copy()
    aaaaa=np.median(bbbbb.transpose(),axis=0)
    %varexp --plot bbbbb 
    len(arrrxx)
except:
    len(arrrxx)
len(arrrxx)

dNIt=8
Ngroup=3
Lo=1
Nproc=3*Ngroup#*(os.cpu_count())
wrkdir = r""
[hhhao,Arr_AAA]=(hkl.load(wrkdir + "BTC-USD"+".rlf1"))
NIter=100

for hhhai in range(hhhao):  
    hhha=hhhai+1
    for iGr in range(Ngroup):  
        ZDat=Arr_AAA[iGr*NIter*int(Nproc/Ngroup)+max(0,(hhha)-dNIt)*int(Nproc/Ngroup):iGr*NIter*int(Nproc/Ngroup)+(hhha)*int(Nproc/Ngroup)].copy()
        if iGr==0:
            xxxx=ZDat.copy()
        else:
            xxxx=np.concatenate((xxxx, ZDat))
    ZDat=xxxx.copy()                     
    
    anI=len(ZDat)
    for ii in range(3):       
        ar0x=np.median(ZDat,axis=0)
        ar0x_=.4*(np.median(abs((ZDat)-(ar0x)),axis=0))
            
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
        for i in range(anI):  
            ZDat[i]= savgol_filter(ZDat[i], 14, 5)
        P=np.zeros(3,float)
        for i in range(anI):
            dd=ZDat[i][lnn-NNew:].copy()                         
            x=ar0x[lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))].copy()
            ZDat[i][lnn-NNew:]=filterFourierQ(ZDat[i],(ar0x),NNew,1)[lnn-NNew:]
            P[0:2]=np.polyfit(x,ZDat[i][lnn-NNew:lnn-NNew+int(lSrez*(NNew-(lnn-len(ar0x))))],1)
            if not P[0]>0:
                P[0:2]=np.polyfit(dd,ZDat[i][lnn-NNew:],1)
            ZDat[i][lnn-NNew:]=(ZDat[i][lnn-NNew:]-P[1])/P[0]                      
    #ZDat=Arr_AAA[iGr*NIter*int(Nproc/Ngroup)+max(0,(hhh+1)-dNIt)*int(Nproc/Ngroup):iGr*NIter*int(Nproc/Ngroup)+(hhh+1)*int(Nproc/Ngroup)].copy()
    bbbbb=ZDat[:,:].transpose().copy()
    aaaaa=np.median(bbbbb.transpose(),axis=0)
    %varexp --plot bbbbb 


                                    # for jj in range(aMM):    
                                    #     dd1=dd[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)].copy()
                                    #     mdd4_=mdd4[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)].copy()
                                        
                                    #     seqA0=(dd1.reshape(len(dd1)*len(dd1[0])))[1:]*np.ceil(0.5*(1/(mdd4_.reshape(len(dd1)*len(dd1[0]))==1)[0:len(dd1)*len(dd1[0])-1]+1/(mdd4_.reshape(len(dd1)*len(dd1[0]))==1)[1:]))
                                    #     seqA0_=  mdd4_.reshape(len(dd1)*len(dd1[0]))*0
                                    #     seqA0_[0]=1
                                    #     seqA0_[1:]=(abs(seqA0)==np.Inf)+np.isnan(seqA0)
                                    #     seqA0_=seqA0_.reshape((len(dd1),len(dd1[0])))
                                    #     seqA0=np.asarray(list(filter(lambda x: abs(x)!= np.Inf, seqA0)),float) 
                                    #     seqA0=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqA0)),float)    
                                        
                                    #     seqA=(dd1.reshape(len(dd1)*len(dd1[0])))[1:]*np.ceil(0.5*np.fabs(1/(1*(mdd4_.reshape(len(dd1)*len(dd1[0]))==1)[0:len(dd1)*len(dd1[0])-1]-1*(mdd4_.reshape(len(dd1)*len(dd1[0]))==1)[1:])))                                                                                                       
                                    #     seqA_=  mdd4_.reshape(len(dd1)*len(dd1[0]))*0
                                    #     seqA_[0]=1
                                    #     seqA_[1:]=(abs(seqA)==np.Inf)+np.isnan(seqA)
                                    #     seqA_=seqA_.reshape((len(dd1),len(dd1[0])))
                                    #     seqA=np.asarray(list(filter(lambda x: abs(x)!= np.Inf, seqA)),float) 
                                    #     seqA=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqA)),float)    
                                        
                                    #     DD__A=DD_[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)].copy()
                                    #     DD__B=-DD__A[:,::-1].copy()
                                    #     # DD__A=DD__A*(DD__A>0)
                                    #     # DD__B=DD__B*(DD__B>0)
                                        
                                    #     if len(dd1)>1 and len(dd1[0])>=len(dd1):
                                    #         eeB= -(seqA0_)*(DD__A+dd1)-seqA0_*(-XFilter.RALF1FilterX(DD__A+dd1,len(dd1),len(dd1[0]),1,0))#+seqA_*((DD__A))
                                    #         eeA=  (seqA0_)*(DD__B-dd1)+seqA0_*(-XFilter.RALF1FilterX(DD__B-dd1,len(dd1),len(dd1[0]),1,0))#-seqA_*((DD__B))
                                    #         # dd1=(eeA+eeB)/2
                                    #         # eeA=dd1.copy()
                                    #         # eeB=dd1.copy()
                                    #         dd_AA[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]=eeA.copy()#*(eeB>0)*((eeA+eeB)>0)
                                    #         dd_BB[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]=eeB.copy()#*(eeA<0)*((eeA+eeB)<0)
                                      
                                    #     dd2_1=(dd_AA)[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)] 
                                    #     seqB=(dd2_1.reshape(len(dd2_1)*len(dd2_1[0])))[1:]*np.ceil(0.5*np.fabs(1/(1*(mdd4_.reshape(len(dd2_1)*len(dd2_1[0]))==1)[0:len(dd2_1)*len(dd2_1[0])-1]-1*(mdd4_.reshape(len(dd2_1)*len(dd2_1[0]))==1)[1:])))
                                    #     seqB=np.asarray(list(filter(lambda x: abs(x)!= np.Inf, seqB)),float) 
                                    #     seqB=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqB)),float)    
                                    #     dd2_2=(dd_BB)[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]
                                    #     seqC=(dd2_2.reshape(len(dd2_2)*len(dd2_2[0])))[1:]*np.ceil(0.5*np.fabs(1/(1*(mdd4_.reshape(len(dd2_2)*len(dd2_2[0]))==1)[0:len(dd2_2)*len(dd2_2[0])-1]-1*(mdd4_.reshape(len(dd2_2)*len(dd2_2[0]))==1)[1:])))
                                    #     seqC=np.asarray(list(filter(lambda x: abs(x)!= np.Inf, seqC)),float) 
                                    #     seqC=np.asarray(list(filter(lambda x: abs(np.isnan(x))!= 1, seqC)),float)    

                                    #     try:
                                    #         P_1=P.copy()
                                    #         P_2=P.copy()
                                    #         P_1[0:2]=np.polyfit(seqA,seqB,1)
                                    #         P_2[0:2]=np.polyfit(seqA,seqC,1)
                                    #         P_1[0]=np.std(seqB)/np.std(seqA)
                                    #         P_1[1]=np.mean(seqB)-P_1[0]*np.mean(seqA)  
                                    #         P_1[2]=0
                                    #         P_2[0]=np.std(seqC)/np.std(seqA)
                                    #         P_2[1]=np.mean(seqC)-P_2[0]*np.mean(seqA)  
                                    #         P_2[2]=0
                                    #         if 100*scp.pearsonr(seqA,seqB)[0]>0 and 100*scp.pearsonr(seqA,seqC)[0]>0: #and not (abs(P_1[0]-1)>0.5 or abs(P_2[0]-1)>0.5)
                                    #             dd_CC[int(ii*anI/aNN):int((ii+1)*anI/aNN),int(jj*Nf/aMM):int((jj+1)*Nf/aMM)]=0.5*(dd2_1+dd2_2)#0.5*((dd2_1.copy()-P_1[1])/P_1[0]+(dd2_2.copy()-P_2[1])/P_2[0])
                                    #         else:
                                    #             PP=0    
                                    #     except:
                                    #         PP=0    

