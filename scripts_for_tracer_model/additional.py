'''
Additional functions for supporting ray-tracing
'''
import os
import numpy as N

def float_to_binary(data,filename):
    output_file=open(filename,'wb')
    data.tofile(output_file)
    output_file.close()

def binary_to_float(filename):
    '''
    To convert the binary data of rays bundle to float
    filename - directory of the data file
    '''
    float_array=N.fromfile(filename,dtype=float)
    return float_array

def zip_files(filesname,zipdir):
    '''
    zip the saved ray bundles 
    '''        
    os.system('zip -j %s.zip %s'%(zipdir, filesname))
    os.system('rm -R %s'%filesname)

def confidence_interval(x,i, Q, A):
    '''
    calculate confidence interval CI
    ref :https://en.wikipedia.org/wiki/Standard_deviation

    Arguements:
    i - the i-th iteration
    x - the current instant value from the i-th iteration
    Q/A - see in the ref, the initial value should be 0        


    Return:
    ci - the confidence interval
    '''

    Q=Q+float(i)/float(i+1)*(x-A)**2
    
    if i ==0:
        ci=8888
    else:
        sd=(Q/float(i))**0.5
        ci=3.*sd/(float(i+1)**0.5)

    return ci, Q


def preprocess_CSR(csr, source='CA'):   
    '''
    pre proceed CSR to true value
    source - 'CA' from Charles; or 'Tonatiuh' from Tonatiuh
    '''
    if source=='CA':
        if csr<=0.1:
            CSR = -2.245e+03*csr**4.+5.207e+02*csr**3.-3.939e+01*csr**2.+1.891e+00*csr+8e-03
        else:
            CSR = 1.973*csr**4.-2.481*csr**3.+0.607*csr**2.+1.151*csr-0.020   

    elif source=='Tonatiuh':

        if csr>0.145:
            CSR = -0.04419909985804843 + csr * (1.401323894233574 + csr * (-0.3639746714505299 + csr * (-0.9579768560161194 + 1.1550475450828657 * csr)))
		
        elif csr>0.035 and csr <=0.145:
            CSR=0.022652077593662934 + csr * (0.5252380349996234 + (2.5484334534423887 - 0.8763755326550412 * csr)) * csr
		
        else:
            CSR = 0.004733749294807862 + csr* (4.716738065192151 + csr * (-463.506669149804 + csr * ( 24745.88727411664+csr * (-606122.7511711778 + 5521693.445014727 * csr) ) ) )

    
    return CSR

  
def bin_element(m, n, W, H):
    '''
    Generate bins for fluxmap

    m- number of elements on the y direction
    n- number of elements on the x direction
    W- width -x direction
    H- height -y direction
    '''
    dx=W/float(m)
    dy=H/float(n)

    x_bins=N.linspace(-W/2.,W/2.,m+1)
    y_bins=N.linspace(-H/2.,H/2.,n+1)

    bin_area=dx*dy
    bins=N.array([x_bins,y_bins])
    fluxmap=N.zeros((m,n))

    return bins, bin_area, fluxmap


def bin_radius(w, h, p):
    '''
    Generate bins in radius direction

    w - receiver width
    h - receiver height
    p - number of elements in the R direction
    '''
    R=N.sqrt(w**2+h**2)

    r_bins=N.linspace(0.,R, p+1)
    bin_area=N.pi*(r_bins[1:]**2-r_bins[:-1]**2)
    fluxmap=N.zeros(p)

    return r_bins, bin_area, fluxmap








