import numpy as np
import cv2
from scipy import ndimage
import concurrent.futures
import multiprocessing as mp
from scipy.signal import savgol_filter
from RALf1FiltrVID import RALf1FiltrV,RandomV,filterFourierV
import dill                            
wrkdir = r".\\"
wwrkdir_=r".\W1\\"

filename = wwrkdir_+"globalsavepkl"
anamef="fralf.tmp"
fo = open(anamef, "w")
fo.write(str(0)+'\n')
fo.close()  
hhh=0
try:
    dill.load_session(filename+".ralf")  
except:
    hhh=hhh        
    coef=0.2
    NIt=1
    NumSteps=3
    NCircls=10
    Nproc=int(np.floor(mp.cpu_count()/3))
    Limite=260000
        
    # Create a VideoCapture object and read from input file 
    cap = cv2.VideoCapture(wwrkdir_ +'coronavx.mp4')#or mp4     
    ArrX=[[] for i in range(3)]
    NumFr=int(cap.get(cv2.CAP_PROP_FRAME_COUNT))#CAP_PROP_POS_MSEC CAP_PROP_FRAME_COUNT 
    aDur=int(cap.get(cv2.CAP_PROP_FPS))
    sz1=int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    sz2=int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    gray_sz1=sz1
    while(cap.isOpened()):
        ret, frame = cap.read()    
        if ret==True:
            frame_=frame
            for icl in range(3):
                gray=np.zeros((sz1,sz2),np.uint8)
                for i in range(sz1):
                    for j in range(sz2):
                        gray[i][j]=frame[i][j][icl]
                cv2.imshow('frame', gray)
                gray = ndimage.zoom(gray, coef)        
                ArrX[icl].append(gray)    
                # Press Q on keyboard to  exit 
                if cv2.waitKey(30) & 0xFF == ord('q'): 
                    break   
        else:
            break        
    cap.release() 
    cv2.destroyAllWindows() 
    Arr=np.asarray(ArrX,float)    
    sz1=len(gray)
    sz2=len(gray[0])
    NNew0=int(NumFr*0.4)
    dill.dump_session(filename+".ralf")

while hhh<NumSteps:
    ArrRezMx=np.zeros((3,NumFr,sz1,sz2),float)-np.Inf
    ArrRezMn=np.zeros((3,NumFr,sz1,sz2),float)+np.Inf
    hh=0
    try:
        dill.load_session(filename+("%s_%s"%(hhh,hh-1))+".ralf")        
    except:
        hh=hh
    NNew=int(NNew0*(1-hhh/(NumSteps)))
    Nn0=NumFr-NNew+int(NNew*(1/(NumSteps)))
    while hh<NCircls: 
        SZ=int(sz1*sz2)
        NumFri=RandomV(SZ) 
        NumFrX_=np.zeros(SZ,int)
        NumFrY_=np.zeros(SZ,int)
        for i in range(SZ):
            NumFrY_[i]= np.floor(NumFri[i]/sz1)
            NumFrX_[i]=NumFri[i]-NumFrY_[i]*sz1 
        Ndel=(SZ*NumFr/Limite)
        NChan=int(np.ceil(sz1/Ndel))     
        SZ=NChan*int(np.ceil(sz2*Ndel))
        Arr_=np.zeros((3,SZ,int(NumFr)),float)
        
        for kk in range(int(SZ/(sz1*sz2))+1):
            if kk==0:
                NumFrX= NumFrX_.copy()
                NumFrY=NumFrY_.copy()  
            else:
                NumFrX=np.concatenate((NumFrX, NumFrX_))
                NumFrY=np.concatenate((NumFrY, NumFrY_))    
        for icl in range(3):    
            for i in range(int(NumFr)):
                for j in range(SZ):
                    Arr_[icl][j][i]=Arr[icl][i][NumFrX[j]][NumFrY[j]]            
            
        for l in range(int(np.ceil(sz2*Ndel))):         
            Arr_x=np.zeros((3,NChan,NumFr),float)
            argss=[[] for kk in range(3*Nproc)]
            for kk in range(3*Nproc):
                argss[kk]=[0, "%s"%NChan, "%s"%NNew, "%s"%NIt]
                icl=int(kk/Nproc)
                Arr_x[icl]=Arr_[icl][l*NChan:(l+1)*NChan].copy()    
                for i in range(NChan): 
                    for j in range(NumFr):
                        if j>NumFr-NNew-1:
                            argss[kk].append(0)
                        else:
                            argss[kk].append(Arr_x[icl][i][j]-Arr_x[icl][i][NumFr-NNew-1])                                                       
                    
            if __name__ == '__main__':   
                arezAMx=[]
                # kk=0
                # dd=RALf1FiltrV(argss[kk])
                # for icl in range(3*Nproc):
                #     arezAMx.append(dd)
                executor=concurrent.futures.ThreadPoolExecutor(max_workers=Nproc)
                future_to = {executor.submit(RALf1FiltrV, argss[kk]) for kk in range(3*Nproc)}
                for future in concurrent.futures.as_completed(future_to):                
                    arezAMx.append(future.result())
                del(executor)            
            
            arezAMx=np.asarray(arezAMx,float)
            for icl in range(3):   
                arezAMx_=arezAMx[0+icl*Nproc:Nproc+icl*Nproc].copy()
                arezBMx=np.zeros((NChan,NumFr),float)            
                arezAMxZ=(np.max(arezAMx_,axis=0)+np.min(arezAMx_,axis=0))/2
                
                arezAMxZ=filterFourierV(arezAMxZ,np.reshape(Arr_x[icl],len(Arr_x[icl])*len(Arr_x[icl][0])),NNew,NChan)
                
                for i in range(NChan):
                    arezBMx[i]=arezAMxZ[0+NumFr*i:NumFr+NumFr*i].copy()              
                    # arezBMx[i][NumFr-NNew:NumFr]=arezBMx[i][NumFr-NNew:NumFr]*np.std(Arr_x[icl][i][NumFr-NNew:Nn0])/np.std(arezBMx[i][NumFr-NNew:Nn0])                           
                    # arezBMx[i][NumFr-NNew:NumFr]=arezBMx[i][NumFr-NNew:NumFr]-np.mean(arezBMx[i][NumFr-NNew:Nn0])+np.mean(Arr_x[icl][i][NumFr-NNew:Nn0])
                    arezBMx[i][NumFr-NNew:NumFr]=arezBMx[i][NumFr-NNew:NumFr]+Arr_x[icl][i][NumFr-NNew-1]
                    #arezBMx[i][0:NumFr-NNew]=Arr_0[icl][i][0:NumFr-NNew].copy()
                    #arezBMx[i]= savgol_filter(arezBMx[i], 11, 5) 
                    arezBMx[i][0:NumFr-NNew]=Arr_x[icl][i][0:NumFr-NNew].copy()
                   
                for j in range(l*NChan,(l+1)*NChan):
                    for i in range(int(NumFr)):
                        ArrRezMx[icl][i][NumFrX[j]][NumFrY[j]]=max(ArrRezMx[icl][i][NumFrX[j]][NumFrY[j]],arezBMx[j-l*NChan][i])
                        ArrRezMn[icl][i][NumFrX[j]][NumFrY[j]]=min(ArrRezMn[icl][i][NumFrX[j]][NumFrY[j]],arezBMx[j-l*NChan][i])                 
    
            print ("calculated %s percents"%(int((l+1)/(int(np.ceil(sz2*Ndel)))*1000)/10))
        
        ArrRez=(ArrRezMx+ArrRezMn)/2      
        if hh==0:
            ArrRez0=ArrRez
        else:
            ArrRez0=(ArrRez0*hh+ArrRez)/(hh+1) 
        ZZ=ArrRez0.copy()
        for icl in range(3):
            ffZZ=[]
            for l in range(NumFr-NNew):
                ffZZ.append(np.fft.fft2(ZZ[icl][l]-ZZ[icl][NumFr-NNew-1]))
            ffZZ=np.asarray(np.abs(ffZZ),float)
            affZZ=np.max(ffZZ,0)
            for l in range(NumFr-NNew,NumFr):
                ZZ_=np.fft.fft2(ZZ[icl][l]-ZZ[icl][NumFr-NNew-1])
                aZZ_=np.abs(ZZ_)+1e-32
                mZZ_=np.mean(aZZ_)                
                fZZ=np.fft.ifft2((ZZ_/aZZ_)*(1*(aZZ_>0.62*mZZ_))*affZZ)
                ZZ[icl][l]=fZZ.real+ZZ[icl][NumFr-NNew-1]   
                        
        ArrRez_=ZZ.copy()
        ArrRez_=np.asarray(ArrRez_-(ArrRez_-255)*(ArrRez_>255),np.uint8)
        ArrRez_=np.asarray(ArrRez_-ArrRez_*(ArrRez_<0),np.uint8)
        coefX=gray_sz1/sz1
        out = cv2.VideoWriter(wwrkdir_ +'coronavout.avi',cv2.VideoWriter_fourcc('M','J','P','G'), aDur, (int(sz1*coefX),int(sz2*coefX)))
        for l in range(int(NumFr)):
            frame=frame_
            for icl in range(3):        
                frame[ : , : , icl] = ndimage.zoom(ArrRez_[icl][l], coefX)  
            #frame=cv2.medianBlur(frame,21)
            cv2.imshow('frame', frame)
            # Press Q on keyboard to  exit 
            if cv2.waitKey(300) & 0xFF == ord('q'): 
                break    
            out.write(frame)
        out.release()
        cv2.destroyAllWindows()
        hh=hh+1
        dill.dump_session(filename+("%s_%s"%(hhh,hh-1))+".ralf")              
   
    cap = cv2.VideoCapture(wwrkdir_ +'coronavout.avi')#or mp4     
    ArrX=[[] for i in range(3)]
    NumFr=int(cap.get(cv2.CAP_PROP_FRAME_COUNT))#CAP_PROP_POS_MSEC CAP_PROP_FRAME_COUNT 
    aDur=int(cap.get(cv2.CAP_PROP_FPS))
    sz1=int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    sz2=int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    gray_sz1=sz1
    while(cap.isOpened()):
        ret, frame = cap.read()    
        if ret==True:
            frame_=frame
            for icl in range(3):
                gray=np.zeros((sz1,sz2),np.uint8)
                for i in range(sz1):
                    for j in range(sz2):
                        gray[i][j]=frame[i][j][icl]
                gray = ndimage.zoom(gray, coef)        
                ArrX[icl].append(gray) 
        else:
            break        
    cap.release() 
    Arr=np.asarray(ArrX,float)        
    sz1=len(gray)
    sz2=len(gray[0])
    NNew0=int(NumFr*0.4)    
    hhh=hhh+1
    dill.dump_session(filename+".ralf")    
    
