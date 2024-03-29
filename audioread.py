import concurrent.futures
import pylab as plt 
import numpy as np
import pandas as pd
from PIL import Image
import msvcrt
import cv2 as cv
from RALf1FiltrVID import RALf1FiltrQ,filterFourierQ
import moviepy.editor as mpv
from scipy import stats as scp
  
wrkdir = r".\\"
wwrkdir_=r".\W10\\"
nama='explosion'
nmfile0=nama+'.mp4'
nmfile=nama+'out.mp4'
filename = wwrkdir_+"globalsavepkl"

Lengt=300
Ngroup=6
Nproc=Ngroup*4#(mp.cpu_count()-1)
Lo=0
aTmStop=2
NIt=4
NIter=60
DT=0.1
Nf_K=3
dsiz=5000

if __name__ == '__main__': 
    ImApp=[]       
    audio = mpv.AudioFileClip(wwrkdir_ +nama+".mp4")
    samplerate, samples = audio.fps, np.asarray(audio.to_soundarray(),float)
    siz=len(samples)

    interv="%d"%dsiz
    siz_=int(siz/dsiz)
    samplesx=np.zeros(siz_,float)
    for i in range(siz_):            
        samplesx[i]=(((np.mean(samples[i*dsiz+0:(i+1)*dsiz,0]))))
    samplesx=np.log(abs(samplesx))
    samplesx_=np.mean(np.asarray(list(filter((-np.Inf).__ne__, samplesx)),float))
    samplesx=samplesx-samplesx_
    siz=len(samplesx)
    samplesy=np.concatenate((samplesx,np.zeros(siz,float)))
    samplesx=samplesy[np.asarray(range(siz),int)+siz*(samplesx==-np.Inf)]  
    samplesx=samplesx-np.mean(samplesx)+6*np.std(samplesx)
    arrrxx=np.asarray(samplesx,float)
    arrrxx_=[]
    arrrxx_.append(arrrxx)
    
    esz=len(arrrxx_)
    arr_rezDzRez=[[] for j in range(esz)]
    kkk=0
    out=0
    while kkk <esz:        
        arrr=arrrxx_[kkk].copy()    
        arrr=np.asarray(arrr,float)    
        Lengt=min(Lengt,arrr.size)
        
        Nf=Lengt
        arrr=arrr[0:Nf]
        
        nn=int(Nf*DT)             
        NNew=int(Nf/2)  
        Nf=Nf+nn        
        ar0=np.asarray(arrr[0:])           
        
        arr_z=np.zeros(Nf,float)
        arr_z[0:Nf-NNew]=arrr[0:Nf-NNew].copy()
        arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]
          
        aname=nama
        adat0=Nf-NNew
        
        fig = plt.figure()
        axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
        axes_ = fig.add_axes([0, 0, 0.3, 0.3])     
        axes.plot(ar0, 'r.')
        axes.plot(arr_z, 'go-')  #cut data used for model
        axes.text(4, 4, 'Course = %s, start = %s, step = %s'%(aname,adat0,interv),
                verticalalignment='bottom', horizontalalignment='right',
                transform=axes_.transAxes,color='blue', fontsize=14)        
        fig.savefig(wrkdir +'dynamic.png',dpi=300,transparent=False,bbox_inches = 'tight')
        frame=Image.open(wrkdir +'dynamic.png')
        ImApp.append(frame)
        cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)        
        gray_sz1=len(cimg[0])
        gray_sz2=len(cimg)
        aDur=4
        fourcc = cv.VideoWriter_fourcc(*'MP4V')
        out = cv.VideoWriter(wrkdir + aname+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
        for icl in range(len(ImApp)):
            cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
            out.write(cimgx[0:gray_sz2,0:gray_sz1,:]) 
        out.release()
        plt.show()
       
        arr_rezBz=np.zeros(Nf,float)
        mn=np.mean(arr_z[0:Nf-NNew])
        arr_rezMx=  np.zeros((Ngroup,Nf),float)
        arr_rezMn=  np.zeros((Ngroup,Nf),float)
        Arr_AM=   np.asarray([[[0 for k in range(Nf)] for j in range(NIter)] for i in range(Nproc)],float)
        Arr_BM=   np.asarray([[[0 for k in range(Nf)] for j in range(NIter)] for i in range(Nproc)],float)
        all_rezAz=np.zeros((NIter,Nf),float)
        arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]  
        argss=[[0] for j in range(Nproc)]    
        hh=0
        hhh=0
        hhh_=0
        Koef=-1
        ZZ=0
        key=0
        TKoef=-100
        
        nnn=int(nn/4)
        arr_rezCz=np.zeros(int(nnn*(aTmStop-1)),float)
        while hhh_<aTmStop and not key == 13: 
            Aprocess=[]
            if hhh==NIter:
                if hhh_<aTmStop-1:
                    TKoef=-100
                    hhh=0
                    arr_rezCz[hhh_*nnn:(hhh_+1)*nnn]=arr_rezBz[Nf-NNew:Nf-NNew+nnn].copy()
                    arr_z=np.zeros(Nf+nnn,float)
                    arr_z[0:nnn]=arr_rezBz[0:nnn].copy()
                    arr_z[nnn:Nf]=arr_rezBz[nnn:Nf].copy()
                    arr_z[Nf:]=arr_z[Nf-1]
                    Lengt=Lengt+nnn
                    ar0=arr_z[0:Lengt-nnn].copy()
                    Nf=Nf+nnn
                    arr_rezBz=np.zeros(Nf,float)
                    mn=np.mean(arr_z[0:Nf-NNew])
                    arr_rezMx=  np.zeros((Ngroup,Nf),float)
                    arr_rezMn=  np.zeros((Ngroup,Nf),float)
                    Arr_AM=   np.asarray([[[0 for k in range(Nf)] for j in range(NIter)] for i in range(Nproc)],float)
                    Arr_BM=   np.asarray([[[0 for k in range(Nf)] for j in range(NIter)] for i in range(Nproc)],float)
                    all_rezAz=np.zeros((NIter,Nf),float)
                    arr_z[Nf-NNew:]=arr_z[Nf-NNew-1]  
                    hhh_=hhh_+1
                else:
                    hhh_=hhh_+1
                    ZZ=1
            if ZZ==0:                                                          
                for iProc in range(Nproc):
                    arr_A=arr_z.copy()
                    if Lo:
                        arr_A=np.log(arr_z)
                    Asr=np.mean(arr_A[0:Nf-NNew])
                    arr_A=arr_A-Asr
                    Klg=np.power(10,np.floor(np.log10(np.max(abs(arr_A)))))
                    arr_A=arr_A/Klg
        
                    if hh==0 and iProc==0:
                        file=open('datapy%d.tmp'%(iProc+1),'w')
                        for i in range(Nf):
                            file.write(str(arr_A[i])+'\n')
                        file.close()
        
                        file=open("anum.tmp",'wt')
                        file.write("%d\n"%(NNew))
                        file.write("%d"%(iProc+1))
                        file.close()
        
                    program =wrkdir + "RALF1FiltrX_lg.py"
                    NChan=1
                    argss[iProc]=["python", "%s"%NChan, "%s"%NNew, "%s"%NIt]#"%s"%(iProc+1)]
                    for i in range(Nf):
                        argss[iProc].append(str("%1.3f"%(arr_A[i])))
                        
                arezAMx=[]
                # for iProc in range(Nproc):
                #     arezAMx.append(RALf1FiltrQ(argss[iProc]))

                with concurrent.futures.ThreadPoolExecutor(max_workers=Nproc) as executor:
                    future_to = {executor.submit(RALf1FiltrQ, argss[iProc]) for iProc in range(Nproc)}
                    for future in concurrent.futures.as_completed(future_to):                
                        arezAMx.append(future.result())
        
    
        
                arezAMx= np.asarray(arezAMx,float)*Klg+Asr

                        
                Arr_AAA=np.zeros((Ngroup,int(Nproc*(hhh+1)/Ngroup),Nf),float)  
                for iGr in range(Ngroup):
                    for iProc in range(int(Nproc/Ngroup)):                
                        Arr_AM[iProc+int(iGr*(Nproc/Ngroup))][hhh]= np.asarray(arezAMx[iProc+int(iGr*(Nproc/Ngroup))],float)
                        Arr_AAA[iGr][0*(hhh+1)+iProc*(hhh+1):1*(hhh+1)+iProc*(hhh+1)]=Arr_AM[iProc+int(iGr*(Nproc/Ngroup))][0:hhh+1].copy()
                    
                    Arr_BBB=Arr_AAA[iGr].transpose()
                    for i in range(Nf):
                        arr_rezMx[iGr][i]=max(Arr_BBB[i])
                        arr_rezMn[iGr][i]=min(Arr_BBB[i])
                        
                aMx=arr_rezMx.transpose()
                aMn=arr_rezMn.transpose()                 
                for i in range(Nf-NNew,Nf):
                    arr_rezBz[i]=(max(aMx[i])+min(aMn[i]))/2
                
                if Lo:           
                    arr_rezBz=filterFourierQ(arr_rezBz,np.log(arr_z),NNew,1)
                    arr_rezBz[0:Nf-NNew]=np.log(ar0[0:Nf-NNew]) 
                else:
                    arr_rezBz=filterFourierQ(arr_rezBz,arr_z,NNew,1)  
                    arr_rezBz[0:Nf-NNew]=ar0[0:Nf-NNew].copy()                    
                
                if Lo:
                    #arr_rezBz[Nf-NNew:Nf]=arr_rezBz[Nf-NNew:Nf]*np.std(np.log(ar0[Nf-NNew:len(ar0)-int((len(ar0)-(Nf-NNew))/2)]))/np.std(arr_rezBz[Nf-NNew:len(ar0)-int((len(ar0)-(Nf-NNew))/2)])
                    arr_rezBz[Nf-NNew:Nf]=arr_rezBz[Nf-NNew:Nf]-np.mean(arr_rezBz[Nf-NNew:len(ar0)-int((len(ar0)-(Nf-NNew))/2)])+np.mean(np.log(ar0[Nf-NNew:len(ar0)-int((len(ar0)-(Nf-NNew))/2)]))
     
                else: 
                    #arr_rezBz[Nf-NNew:Nf]=arr_rezBz[Nf-NNew:Nf]*np.std(ar0[Nf-NNew:len(ar0)-int((len(ar0)-(Nf-NNew))/2)])/np.std(arr_rezBz[Nf-NNew:len(ar0)-int((len(ar0)-(Nf-NNew))/2)])                        
                    arr_rezBz[Nf-NNew:Nf]=arr_rezBz[Nf-NNew:Nf]-np.mean(arr_rezBz[Nf-NNew:len(ar0)-int((len(ar0)-(Nf-NNew))/2)])+np.mean((ar0[Nf-NNew:len(ar0)-int((len(ar0)-(Nf-NNew))/2)]))
                
                all_rezAz[hhh]=arr_rezBz.copy()        
                all_rezAz_=all_rezAz[0:hhh+1].transpose()                
                for i in range(Nf):
                    arr_rezBz[i]=np.mean(all_rezAz_[i])
                    
                if Lo:                       
                    arr_rezBz=np.exp(arr_rezBz) 
                    all_rezAz_=np.exp(all_rezAz_)

                mm1=ar0[Nf-NNew:].copy()                            
                mm2=arr_rezBz[Nf-NNew:len(ar0)].copy()   
                if np.std(mm1)>0 and np.std(mm2)>0:
                    Koef=100*scp.spearmanr(mm1,mm2)[0]
                else:
                    Koef=-2
                if (Koef+100)> .62*(TKoef+100):                                        
                    fig = plt.figure()
                    axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
                    axes_ = fig.add_axes([0, 0, 0.3, 0.3])   
                    axes.plot(all_rezAz_,'oy',alpha=0.1)
                    axes.plot(ar0, 'ro-', alpha=0.1)
                    axes.plot(arr_rezBz,'cx-', alpha=0.5)
                    axes.text(4, 4, 'Course = %s, start = %s, step = %s'%(aname,adat0,interv),
                            verticalalignment='bottom', horizontalalignment='right',
                            transform=axes_.transAxes,color='blue', fontsize=14)        
                    axes_.plot(mm1,mm2, 'ok', markersize=3, alpha=0.1)               
                    axes_.text(0.2, 0.6, '%d'%int(np.floor(Koef)),
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes_.transAxes,color='green', fontsize=30)
                    fig.savefig(wrkdir +'dynamic.png',dpi=300,transparent=False,bbox_inches = 'tight')
                    frame=Image.open(wrkdir +'dynamic.png')
                    cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)        
                    gray_sz1=min(gray_sz1,len(cimg[0]))
                    gray_sz2=min(gray_sz2,len(cimg))
                    ImApp.append(frame)
                    out = cv.VideoWriter(wrkdir + aname+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
                    for icl in range(len(ImApp)):
                        cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
                        out.write(cimgx[0:gray_sz2,0:gray_sz1,:]) 
                    out.release()
                    plt.show()
                    hh=hh+1
                    hhh=hhh+1
                    print (hh)
                    #arr_z=arr_rezBz.copy()
                if msvcrt.kbhit():              
                    key = ord(msvcrt.getch())  
    
        arr_rezDz=[]
        file=open(wrkdir + aname+'_rezpy.txt','w')
        for i in range(len(ar0)):
            arr_rezDz.append(ar0[i])
        for i in range(len(arr_rezCz)):
            arr_rezDz.append(arr_rezCz[i])
        k=i
        for i in range(Nf-NNew,Nf):    
            arr_rezDz.append(arr_rezBz[i])            
        arr_rezDz=np.asarray(arr_rezDz,float)
        aln=len(arr_rezDz)-NNew+nnn
        for i in range(aln):     
            file.write(str(arr_rezDz[i])+'    ')
            file.write("\n")
        file.close()
        arr_rezDzRez[kkk]=arr_rezDz.copy()
        df = pd.DataFrame(arr_rezDzRez)
        df.to_excel (wrkdir +r'export_traces.xlsx', index = None, header=False) 
        
        fig = plt.figure()
        axes = fig.add_axes([0.1, 0.1, 1.2, 1.2])
        axes_ = fig.add_axes([0, 0, 0.3, 0.3])
        axes.plot(arrr, 'ro-')
        axes.plot(arr_rezDz, 'cx-', alpha=0.5) #predicted data
        axes.text(4, 4, 'Course = %s, start = %s, step = %s'%(aname,adat0,interv),
                verticalalignment='bottom', horizontalalignment='right',
                transform=axes_.transAxes,color='blue', fontsize=14)                
        axes_.plot(mm1,mm2, 'ok', markersize=3, alpha=0.1)               
        fig.savefig(wrkdir +'dynamic.png',dpi=300,transparent=False,bbox_inches = 'tight')
        frame=Image.open(wrkdir +'dynamic.png')
        cimg = cv.cvtColor(np.array(frame), cv.COLOR_RGB2BGR)        
        gray_sz1=min(gray_sz1,len(cimg[0]))
        gray_sz2=min(gray_sz2,len(cimg))
        for icl in range(10):
            ImApp.append(frame)
        out = cv.VideoWriter(wrkdir + aname+'.mp4',fourcc, aDur, (gray_sz1,gray_sz2))                   
        for icl in range(len(ImApp)):
            cimgx=(cv.cvtColor(np.array(ImApp[icl]), cv.COLOR_RGB2BGR)) 
            out.write(cimgx[0:gray_sz2,0:gray_sz1,:])       
        out.release()
        plt.show()
        kkk=kkk+1    
   
