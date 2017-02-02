################################################################################
# program: vb_utilities.py
# author: Tom Irvine
# Email: tom@vibrationdata.com
# version: 2.2
# date: January 27, 2015
# description:  utility functions
#
################################################################################

from __future__ import print_function

import sys

if sys.version_info[0] == 2:
    print ("Python 2.x")
    from tkFileDialog import askopenfilename 
        
if sys.version_info[0] == 3:
    print ("Python 3.x")    
    from tkinter.filedialog import askopenfilename   


import os
import re
import numpy as np

from scipy import stats

from scipy.signal import lfilter
from math import pi,cos,sin,tan,sqrt

import matplotlib.pyplot as plt

from matplotlib.ticker import ScalarFormatter

###############################################################################

def loglog_plot(f,a,xlab,ylab,f1,f2,title_string,fig_num):

    plt.ion() 

    plt.close(fig_num)        
    plt.figure(fig_num)            
        
    plt.plot(f,a)
    plt.title(title_string)
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.grid(which='both')
    plt.xscale('log')
    plt.yscale('log')
        
    if(abs(f1-10)<0.5 and abs(f2-2000)<4):
            
        ax=plt.gca().xaxis
        ax.set_major_formatter(ScalarFormatter())
        plt.ticklabel_format(style='plain', axis='x', scilimits=(f1,f2))    
              
        extraticks=[10,2000]
        plt.xticks(list(plt.xticks()[0]) + extraticks) 

        
    if(abs(f1-20)<0.5 and abs(f2-2000)<4):
            
        ax=plt.gca().xaxis
        ax.set_major_formatter(ScalarFormatter())
        plt.ticklabel_format(style='plain', axis='x', scilimits=(f1,f2))    
              
        extraticks=[20,2000]
        plt.xticks(list(plt.xticks()[0]) + extraticks)                
        
        
    plt.xlim([f1,f2])        
        
    plt.show()

    fig_num=fig_num+1

    return fig_num

###############################################################################

def interpolate_psd(f,a,s,df):
    
    fmax=max(f)
    
    num=len(f)
    
    i=0
    
    fi=[]
    ai=[]
    
    while(1): 	

        nf=f[0]+df*i
        
        if( nf > fmax ):
            break
        else:
            fi.append(nf)
            i+=1
   
    m=len(fi)
    
    jf=0
    
    for i in range(0,m):
    
        for j in range (jf,num-1):
            
#            print(' i=%d j=%d m=%d num=%d' %(i,j,m,num))
		
            if( ( fi[i] >= f[j] ) and ( fi[i] <= f[j+1] )  ):
                
                q=a[j]*( ( fi[i] / f[j] )**s[j] )
					
                ai.append(q)
                
                jf=j                 
                
                break
    
    return fi,ai
    
###############################################################################    

def spectral_moments_constant_df(fi,aa_psd,df):
    
    m0=0
    m1=0
    m2=0
    m4=0

    for j in range (0,len(fi)):
	
        m0=m0+aa_psd[j]
        m1=m1+aa_psd[j]*fi[j]
        m2=m2+aa_psd[j]*fi[j]**2
        m4=m4+aa_psd[j]*fi[j]**4


    m0=(m0*df)
    m1=(m1*df)
    m2=(m2*df)
    m4=(m4*df)

    vo=sqrt(m2/m0)
    vp=sqrt(m4/m2)

    return m0,m1,m2,m4,vo,vp

###############################################################################

def half_cosine_fade_perc(y,fper):

    fper=fper/100.

    n=len(y)

    na=int(np.ceil(fper*n))
    nb=n-na
    delta=n-1-nb
    
#    print( 'n=%d na=%d nb=%d fper=%g'  %(n,na,nb,fper))

    for i in range(0,na):
        arg=pi*(( (i-1)/float(na-1) )+1) 
        y[i]=y[i]*0.5*(1+(cos(arg)))


    for i in range(nb,n):
        arg=pi*( (i-nb)/float(delta) )
        y[i]=y[i]*(1+cos(arg))*0.5
            
    return y    


###############################################################################

def integrate_th(num,b,dt):
    
    v=np.zeros(num,'f')    
    
    v[0]=b[0]*dt/2

    for i in range(1,int(num-1)):
        v[i]=v[i-1]+b[i]*dt
    
    v[num-1]=v[num-2]+b[num-1]*dt/2
            
    return v        
            

###############################################################################

def read_one_column_from_dialog(label,pt):
    """
    Prompt the user for the input filename.
    The input file must have one column.
    The input file may have an arbitrary number of header and blank lines.
    Return the column as array b.
    Return the total numbers of lines as num.
    """

    while(1):

        input_file_path = askopenfilename(parent=pt,title=label)

        file_path = input_file_path.rstrip('\n')
#
        if not os.path.exists(file_path):
            print ("This file doesn't exist")
#
        if os.path.exists(file_path):
            print ("This file exists")
            print (" ")
            infile = open(file_path,"rb")
            lines = infile.readlines()
            infile.close()

            b = []
            num=0
            for line in lines:
#
                if sys.version_info[0] == 3:            
                    line = line.decode(encoding='UTF-8') 
                    
                if re.search(r"(\d+)", line):  # matches a digit
                    iflag=0
                else:
                    iflag=1 # did not find digit
#
                if re.search(r"#", line):
                    iflag=1
#
                if iflag==0:
                    line=line.lower()
                    if re.search(r"([a-d])([f-z])", line):  # ignore header lines
                        iflag=1
                    else:
                        line = line.replace(","," ")
                        b.append(float(line))
                        num=num+1
            break;

            b=np.array(b)

            print ("\n samples = %d " % num)
            
    return b,num

###############################################################################

def read_two_columns_from_dialog(label,pt):
    """
    Read data from file using a dialog box
    """ 
    while(1):

        input_file_path = askopenfilename(parent=pt,title=label)

        file_path = input_file_path.rstrip('\n')
#
        if not os.path.exists(file_path):
            print ("This file doesn't exist")
#
        if os.path.exists(file_path):
            print ("This file exists")
            print (" ")
            infile = open(file_path,"rb")
            lines = infile.readlines()
            infile.close()

            a = []
            b = []
            num=0
            for line in lines:
#
                if sys.version_info[0] == 3:            
                    line = line.decode(encoding='UTF-8') 
            
                if re.search(r"(\d+)", line):  # matches a digit
                    iflag=0
                else:
                    iflag=1 # did not find digit
#
                if re.search(r"#", line):
                    iflag=1
#
                if iflag==0:
                    line=line.lower()
                    if re.search(r"([a-d])([f-z])", line):  # ignore header lines
                        iflag=1
                    else:
                        line = line.replace(","," ")
                        col1,col2=line.split()
                        a.append(float(col1))
                        b.append(float(col2))
                        num=num+1
            break

            a=np.array(a)
            b=np.array(b)

            print ("\n samples = %d " % num)
            
    return a,b,num
    
###############################################################################

def read_three_columns_from_dialog(label,pt):
    """
    Read data from file using a dialog box
    """ 
    while(1):

        input_file_path = askopenfilename(parent=pt,title=label)

        file_path = input_file_path.rstrip('\n')
#
        if not os.path.exists(file_path):
            print ("This file doesn't exist")
#
        if os.path.exists(file_path):
            print ("This file exists")
            print (" ")
            infile = open(file_path,"rb")
            lines = infile.readlines()
            infile.close()

            a = []
            b = []
            c = []
            
            num=0
            for line in lines:
#
                if sys.version_info[0] == 3:            
                    line = line.decode(encoding='UTF-8') 
            
                if re.search(r"(\d+)", line):  # matches a digit
                    iflag=0
                else:
                    iflag=1 # did not find digit
#
                if re.search(r"#", line):
                    iflag=1
#
                if iflag==0:
                    line=line.lower()
                    if re.search(r"([a-d])([f-z])", line):  # ignore header lines
                        iflag=1
                    else:
                        line = line.replace(","," ")
                        col1,col2,col3=line.split()
                        a.append(float(col1))
                        b.append(float(col2))
                        c.append(float(col3))
                        num=num+1
            break

            a=np.array(a)
            b=np.array(b)
            c=np.array(c)            

            print ("\n samples = %d " % num)
            
    return a,b,c,num
    
################################################################################

def sample_rate_check(a,b,num,sr,dt):
    dtmin=1e+50
    dtmax=0

    srmin=1./dt
    srmax=1./dt

    for i in range(1, num-1):
        diffA=a[i]-a[i-1]
        if diffA<dtmin:
            dtmin=diffA
        if diffA>dtmax:
            dtmax=diffA

    print ("  dtmin = %8.4g sec" % dtmin)
    print ("     dt = %8.4g sec" % dt)
    print ("  dtmax = %8.4g sec \n" % dtmax)

    if(dtmin>1.0e-20):    
        srmax=float(1./dtmin)
 
    if(dtmax>1.0e-20):   
        srmin=float(1./dtmax)

    print ("  srmax = %8.4g samples/sec" % srmax)
    print ("     sr = %8.4g samples/sec" % sr)
    print ("  srmin = %8.4g samples/sec" % srmin)

    if((srmax-srmin) > 0.01*sr):
        print(" ")
        print(" Warning: sample rate difference ")
#        sr = None
#        while not sr:
#            try:
#                print(" Enter new sample rate ")
#                s = stdin.readline()
#                sr=float(s)
#                dt=1/sr
#            except ValueError:
#                print ('Invalid Number')
    return sr,dt

################################################################################

def WriteData2(nn,aa,bb,output_file_path):
    """
    Write two columns of data to an external ASCII text file
    """
    output_file = output_file_path.rstrip('\n')
    outfile = open(output_file,"w")
    

    if(nn!=len(aa)):
        print(' length error' )
        
    if(nn!=len(bb)):
        print(' length error' )       
        
    
    for i in range (0, nn):
        outfile.write(' %11.7e \t %8.4e \n' %  (aa[i],bb[i]))
    outfile.close()

################################################################################

def WriteData3(nn,aa,bb,cc,output_file_path):
    """
    Write three columns of data to an external ASCII text file
    """
    outfile = open(output_file_path,"w")
    for i in range (0, nn):
        outfile.write(' %8.4e \t %8.4e \t %8.4e \n' %  (aa[i],bb[i],cc[i]))
    outfile.close()
    
################################################################################    

def signal_stats(a,b,num):
    """
    a is the time column.
    b is the amplitude column.
    num is the number of coordinates
    Return
          sr - sample rate
          dt - time step
        mean - average
          sd - standard deviation
         rms - root mean square
        skew - skewness
    kurtosis - peakedness
         dur - duration
    """
    sr=0.
    dt=0.
    ave=0.
    sd=0.
    rms=0.
    skewness=0.
    kurtosis=0.
    dur=0.    
    
    if(len(b)==0):
        print('\n ** Error: len(b)=0 ** \n')
        return sr,dt,ave,sd,rms,skewness,kurtosis,dur
    
    bmax=max(b)
    bmin=min(b)

    ave = np.mean(b)

    dur = a[num-1]-a[0];

    dt=dur/float(num-1)
    sr=1/dt


    rms=np.sqrt(np.var(b))
    sd=np.std(b)

    skewness=stats.skew(b)
    kurtosis=stats.kurtosis(b,fisher=False)
    
    mb=max([bmax,abs(bmin)])
    
    crest=mb/sd

    print ("\n max = %8.4g  min=%8.4g \n" % (bmax,bmin))

    print ("         mean = %8.4g " % ave)
    print ("      std dev = %8.4g " % sd)
    print ("          rms = %8.4g " % rms)
    print ("     skewness = %8.4g " % skewness)
    print ("     kurtosis = %8.4g " % kurtosis)
    print (" crest factor = %8.4g " % crest)

    print ("\n  start = %8.4g sec  end = %8.4g sec" % (a[0],a[num-1]))
    print ("    dur = %8.4g sec \n" % dur)
    return sr,dt,ave,sd,rms,skewness,kurtosis,dur

###############################################################################

class BUTTERWORTH:

    def __init__(self,l,f,fh,fl,dt,iband,iphase,y):
    
        self.l=l
        self.f=f
        self.freq=f
        self.fh=fh
        self.fl=fl
        self.dt=dt
        self.iband=iband
        self.iphase=iphase
        self.y=y
        
        self.om=0
        
        self.a=np.zeros((4,4),'f')	
        self.b=np.zeros((4,4),'f')
        
        self.alpha=np.zeros(2*self.l,'f')
        
        self.s=(1+1j)*np.zeros(20,'f')
        
        self.ns=len(y)
        
        self.ik=0
        
        self.yt=np.zeros(self.ns,'f')        
        
            
    def Butterworth_filter_main(self):
        

        if(self.iband !=3):
            BUTTERWORTH.coefficients(self)
    
        if(self.iband == 1 or self.iband ==2):
            BUTTERWORTH.applymethod(self)

        if(self.iband == 3):
            self.f=self.fh
            self.freq=self.f
            
            print("\n Step 1")
            self.iband=2
    
            BUTTERWORTH.coefficients(self)
            BUTTERWORTH.applymethod(self)

            self.f=self.fl
            self.freq=self.f

            print("\n Step 2")
            self.iband=1
    
            BUTTERWORTH.coefficients(self)
            BUTTERWORTH.applymethod(self)   
            
        return self.y        
            
    @classmethod    
    def coefficients(cls,self):
    
        self.a=np.zeros((4,4),'f')	
        self.b=np.zeros((4,4),'f')		
    
#*** normalize the frequency ***

        targ=pi*self.f*self.dt   # radians
    
        print (" targ = %8.4g " %targ)
             
        self.om=tan(targ)   
    
        print("   om = %8.4g " %self.om)

#*** solve for the poles *******

        BUTTERWORTH.poles(self)

#*** solve for alpha values ****

        print("\n alpha ")    
    
        self.alpha=np.zeros(2*self.l,'f')
        self.alpha=2*self.s.real
    
##    for i in range(0,len(alpha)):
##        print ("  %5.3f +j %5.3f " %(alpha[i].real,alpha[i].imag))

#*** solve for filter coefficients **

        if( self.iband == 1 ):
            BUTTERWORTH.lco(self)
        else:
            BUTTERWORTH.hco(self)
    
#*** plot digital transfer function **

#    dtrans();

#*** check stability ****************
    
        BUTTERWORTH.stab(self)
    
    
    
    @classmethod
    def applymethod(cls,self):
        
        if(self.iphase==1):
            self.apply(self)
            self.apply(self) 
        else:	
            self.apply(self)
        
    

    @classmethod
    def stage1(cls,self):
             
        self.yt=np.zeros(self.ns,'f')

        bc=self.b[self.ik][0:3]
        ac=self.a[self.ik][0:3]
        ac[0]=1
    
        self.yt=lfilter(bc, ac, self.y, axis=-1, zi=None)      


 
    @classmethod
    def stage2(cls,self):
    
        self.y=np.zeros(self.ns,'f')
    
        bc=self.b[self.ik][0:3]
        ac=self.a[self.ik][0:3]
        ac[0]=1

        self.y=lfilter(bc, ac, self.yt, axis=-1, zi=None)  
    
        
    @classmethod
    def apply(cls,self):
        
        BUTTERWORTH.coefficients(self)
 
    
        if(self.iphase==1):	

            yr=np.zeros(self.ns,'f')
            for i in range(0,int(self.ns)):
                yr[self.ns-1-i]=self.y[i]

            self.y=yr
            

#  cascade stage 1

        print("\n  stage 1")
        self.ik=1
        BUTTERWORTH.stage1(self)

#  cascade stage 2

        print("  stage 2");
        self.ik=2
        BUTTERWORTH.stage2(self);
       
#  cascade stage 3

        print("  stage 3");
        self.ik=3
        BUTTERWORTH.stage1(self);
    
        self.y=self.yt

    

    @classmethod	
    def stab(cls,self):
    
        a1=0
        d1=0 
        d2=0 
        d3=0
        dlit=0

        at1=0
        at2=0
        als=0.5e-06
        h2=0

        als*=6.
    
        print ("\n stability reference threshold= %14.7e " %als)

        for i in range(1,int((self.l/2)+1)):
        
            at1= -self.a[i][1]
            at2= -self.a[i][2]

#       print("\n\n stability coordinates: (%12.7g, %14.7g) ",at1,at2);
        
            h2=at2
 
            a1=h2-1.
            d3=at1-a1
         
            a1=1.-h2
            d2=a1-at1
            d1=at2+1.
		
#       print("\n d1=%14.5g  d2=%14.5g  d3=%14.5g",d1,d2,d3);

            dlit=d1

            if(dlit > d2):
                dlit=d2
            if(dlit > d3):
                dlit=d3

            print ("\n stage %ld     dlit= %14.5g " %(i, dlit))

            if(dlit > als):
                print (" good stability")  			
				
            if( (dlit < als) and (dlit > 0.)):		  
                print(" marginally unstable ");
            
            if(dlit < 0.):
                print (" unstable ")	  	
                print ("\n")

################################################################################  
    @classmethod	
    def lco(cls,self):
    
        om2=self.om**2

        for k in range(1,int((self.l/2)+1)):
    
            den = om2-self.alpha[k-1]*self.om+1.
		
            self.a[k][0]=0.
            self.a[k][1]=2.*(om2 -1.)/den
            self.a[k][2]=( om2 +self.alpha[k-1]*self.om+ 1.)/den

            self.b[k][0]=om2/den
            self.b[k][1]=2.*self.b[k][0]
            self.b[k][2]=self.b[k][0]

            print ("\n filter coefficients")		
            print (" a[%i][1]=%10.5g  a[%i][2]=%10.5g" %(k,self.a[k][1],k,self.a[k][2]))
            print (" b[%i][0]=%10.5g  b[%i][1]=%10.5g  b[%i][2]=%10.5g" %(k,self.b[k][0],k,self.b[k][1],k,self.b[k][2]))
    
        print ("\n")
        

################################################################################  
    @classmethod	
    def hco(cls,self):
    
        print ("\n filter coefficients")
    
        om2=self.om**2

        for k in range(1,int((self.l/2)+1)):
    
            den = om2-self.alpha[k-1]*self.om+1.    
		
            self.a[k][0]=0.
            self.a[k][1]=2.*(-1.+ om2)/den
            self.a[k][2]=( 1.+self.alpha[k-1]*self.om+ om2)/den

            self.b[k][0]= 1./den;
            self.b[k][1]=-2.*self.b[k][0]
            self.b[k][2]=    self.b[k][0]
        
            print ("\n a[%i][1]=%10.5g  a[%i][2]=%10.5g" %(k,self.a[k][1],k,self.a[k][2]))
            print (" b[%i][0]=%10.5g  b[%i][1]=%10.5g  b[%i][2]=%10.5g" %(k,self.b[k][0],k,self.b[k][1],k,self.b[k][2]))
            print ("\n")
        

################################################################################  
    @classmethod	
    def poles(cls,self):
        arg=0
        a1=0
        a2=complex(0.,0.)
        h=complex(0.,0.)
        theta=complex(0.,0.)
    
        self.s=(1+1j)*np.zeros(20,'f')    
    
#    print("\n  calculate print ");

        print ("\n poles ")
	
        for k in range(0,int(2*self.l)):
            arg=(2.*(k+1) +self.l-1)*pi/(2.*self.l)
            self.s[k]=cos(arg)+sin(arg)*(1j)
            print (" %4.3f  +j %4.3f " %(self.s[k].real,self.s[k].imag))
         
        for i in range(0,201):   
            arg = i/40.        
        
            h=complex( self.s[0].real,( arg - self.s[0].imag  ))

        for j in range(1,int(self.l)):
            
            theta=complex( -self.s[j].real,( arg - self.s[j].imag ))
            
            temp=h*theta
            h=temp
               
            x=1/h
            h=x
               
            a1 = self.freq*arg
	   
            a2=abs(h)            
           
            a3 = a2**2
            
#            fprint(pFile[3]," %lf %lf %lf \n", a1, a2, a3);     