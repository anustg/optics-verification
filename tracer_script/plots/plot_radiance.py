'''
processing the theoretical solution of case A1 and A2
'''

import numpy as N
import matplotlib.pyplot as plt

def pillbox(sigma, theta=None, plot=True):

    if plot:
        theta=N.linspace(0., 0.01,100000)
    L=N.array([])
    for i in xrange(len(theta)):    
        if theta[i]<=sigma:
            L=N.append(L,1./sigma)

        else:
            L=N.append(L,0.)

    if plot:
        return L, theta
    else:
        return L

def guassian(sigma, theta):
    L=1./(N.sqrt(2.*N.pi)*sigma)*N.exp(-theta**2/2./sigma**2)

    dt=theta[1]-theta[0]
    scale=N.sum(L)*dt

    return L, scale, theta

def buie(theta, CSR):
    L=N.array([])
    for i in xrange(len(theta)):    
        if theta[i] <4.65e-3:
            f=N.cos(0.326*theta[i] *1000.)/N.cos(0.308*theta[i] *1000.)
            L=N.append(L, f)
        else:
            kappa=0.9*N.log(13.5*CSR)*N.power(CSR,-0.3)
            gamma=2.2*N.log(0.52*CSR)*N.power(CSR,0.43)-0.1

            f=N.exp(kappa)*N.power(theta[i] *1000.,gamma)
            L=N.append(L, f)

    dt=theta[1]-theta[0]
    scale=N.sum(L)*dt

    return L, scale, theta

def buie_tonatiuh(theta, CSR):
    if CSR>0.145:
        CSR = -0.04419909985804843 + CSR * (1.401323894233574 + CSR * (-0.3639746714505299 + CSR * (-0.9579768560161194 + 1.1550475450828657 * CSR)))

    elif CSR>0.035 and CSR <=0.145:
        CSR=0.022652077593662934 + CSR * (0.5252380349996234 + (2.5484334534423887 - 0.8763755326550412 * CSR)) * CSR

    else:
        CSR = 0.004733749294807862 + CSR * (4.716738065192151 + CSR * (-463.506669149804 + CSR * ( 24745.88727411664+CSR * (-606122.7511711778 + 5521693.445014727 * CSR ) ) ) )


    L=N.array([])
    for i in xrange(len(theta)):    
        if theta[i] <4.65e-3:
            f=N.cos(0.326*theta[i] *1000.)/N.cos(0.308*theta[i] *1000.)
            L=N.append(L, f)
        else:
            kappa=0.9*N.log(13.5*CSR)*N.power(CSR,-0.3)
            gamma=2.2*N.log(0.52*CSR)*N.power(CSR,0.43)-0.1

            f=N.exp(kappa)*N.power(theta[i] *1000.,gamma)
            L=N.append(L, f)

    dt=theta[1]-theta[0]
    scale=N.sum(L)*dt

    return L, scale, theta



def buie_CA(theta, CSR):
    if CSR<=0.1:
        CSR = -2.245e+03*CSR**4.+5.207e+02*CSR**3.-3.939e+01*CSR**2.+1.891e+00*CSR+8e-03
    else:
        CSR = 1.973*CSR**4.-2.481*CSR**3.+0.607*CSR**2.+1.151*CSR-0.020

    L=N.array([])
    for i in xrange(len(theta)):    
        if theta[i] <4.65e-3:
            f=N.cos(0.326*theta[i] *1000.)/N.cos(0.308*theta[i] *1000.)
            L=N.append(L, f)
        else:
            kappa=0.9*N.log(13.5*CSR)*N.power(CSR,-0.3)
            gamma=2.2*N.log(0.52*CSR)*N.power(CSR,0.43)-0.1

            f=N.exp(kappa)*N.power(theta[i] *1000.,gamma)
            L=N.append(L, f)

    dt=theta[1]-theta[0]
    scale=N.sum(L)*dt

    return L, scale


    
    

def get_radiance(casename, tool):
    '''
    casename - str, eg. '1.2.1'
    tool - str,'tracer', 'tonatiuh','soltrace','solarpilot', 'solstice' or'heliosim'
    
    '''
    data=N.loadtxt('./results/A/%s/A_%s.csv'%(tool, casename), delimiter=',', skiprows=1)

    x=data[:,0].astype(float)
    y=data[:,1].astype(float)
    q=data[:,2].astype(float) #kW/m2


    hst_w=10.
    hst_h=10.
    DNI=1. #kW/m2     
    Q=hst_w*hst_h*DNI #kW

    H=500.
    m=100
    n=100

    # solid angle
    r=N.sqrt(x**2+y**2)
    theta=N.arctan(r/H)

    binarea=hst_w*hst_h/float(m*n)

    R_sq=x**2+y**2+H**2

    solid_angle= binarea*N.cos(theta)/R_sq

    # radiance
    Li=q/solid_angle

    # normalise the radiance
    p=50
    bin_theta=N.linspace(0.,0.012,p+1)
    d_theta=bin_theta[1]-bin_theta[0]

    L=N.zeros(p)
    t=N.zeros(p)

    for j in xrange(len(bin_theta)-1):
        for i in xrange(len(Li)):
            if (theta[i]>=bin_theta[j] and theta[i]<bin_theta[j+1]):
                L[j]=(L[j]*t[j]+Li[i])/(t[j]+1.)
                t[j]+=1.

    th0=N.linspace(d_theta/2., 0.012-d_theta/2., p)

    total=N.sum(L*d_theta)
    L0=L/total

    return L0, th0


def plot_1_1():
    'pillbox slope error'
    tools=N.array([ 'soltrace','tracer', 'solstice', 'heliosim'])
    Tools=N.array([ 'SolTrace', 'Tracer', 'SOLSTICE', 'Heliosim'])
    colors=N.array(['orange','g', 'y','r'])
    symbl=N.array(['o', '^', 's', 'd'])

    idx={}
    idx[0]=N.arange(0,47,4).astype(int)
    idx[1]=N.arange(1,48,4).astype(int)
    idx[2]=N.arange(2,49,4).astype(int)
    idx[3]=N.arange(3,50,4).astype(int)

    plt.figure(1, dpi=1000)  

    for i in xrange(1,4):
        casename='1.1.%s'%i
     
        if i==1:
            for j in xrange(len(tools)):
                tool=tools[j]
                L, theta=get_radiance(casename, tool)
                plt.plot(theta[idx[j]], L[idx[j]], symbl[j], color='black', label=Tools[j])

            sigma=float(i)*1e-3*2.
            theoret, th=pillbox(sigma) 
            plt.plot(th, theoret, '-',color='black', label='Theoretical')

        else:
            for j in xrange(len(tools)):
                tool=tools[j]
                L, theta=get_radiance(casename, tool)
                plt.plot(theta[idx[j]], L[idx[j]], symbl[j], color='black')

            sigma=float(i)*1e-3*2.
            theoret, th=pillbox(sigma) 
            plt.plot(th, theoret, '-',color='black')
    plt.xticks(size=18)   
    plt.yticks(size=18)
    plt.axis(fontsize=18)
    plt.xlabel('$\\theta$', fontsize=18)
    plt.ylabel('Normalised Radiance',  fontsize=18)
    plt.title(' Pillbox slope error', fontsize=18)
    plt.legend(fontsize=18)
    plt.savefig(open('./plots/radiance/A_1.1.png', 'w'), bbox_inches='tight')
    plt.close()

    # plot differences
    for i in xrange(1,4):
        plt.figure(i+1, dpi=1000)
        casename='1.1.%s'%i 
        sigma=float(i)*1e-3*2.
        theoret=pillbox(sigma, theta, False) 

        for j in xrange(len(tools)):
            tool=tools[j]
            L, theta=get_radiance(casename, tool)
            diff=L-theoret
            plt.plot(theta, diff, color=colors[j], label=Tools[j])

        plt.xlabel('$\\theta$')
        plt.ylabel('Differences of the normalised Radiance from simulations \n  to the theoretical value')
        plt.title(' Pillbox slope error \n %s mrad'%i)

        plt.legend()
        plt.savefig(open('./plots/radiance/A_%s_diff.png'%casename, 'w'), bbox_inches='tight')
        plt.close()
   
    
def plot_1_2():
    'normal slope error'
    tools=N.array([ 'tonatiuh','soltrace','tracer', 'solstice', 'heliosim'])
    Tools=N.array([ 'Tonatiuh','SolTrace','Tracer', 'SOLSTICE', 'Heliosim'])
    colors=N.array(['b','orange','g', 'y','r'])
    symbl=N.array(['v','o', '^', 's', 'd'])

    idx={}
    idx[0]=N.arange(0,46,5).astype(int)
    idx[1]=N.arange(1,47,5).astype(int)
    idx[2]=N.arange(2,48,5).astype(int)
    idx[3]=N.arange(3,49,5).astype(int)
    idx[4]=N.arange(4,50,5).astype(int)

    plt.figure(1, dpi=1000)  

    for i in xrange(1,4):
        casename='1.2.%s'%i
     
        if i==1:

            for j in xrange(len(tools)):
                tool=tools[j]
                L, theta=get_radiance(casename, tool)
                sigma=float(i)*1e-3*2.
                theoret, scale, th=guassian(sigma, theta) 
                L*=scale
                plt.plot(theta[idx[j]], L[idx[j]], symbl[j], color='black', label=Tools[j])

            plt.plot(th, theoret, '-',color='black', label='Theoretical')

        else:
 
            for j in xrange(len(tools)):
                tool=tools[j]
                L, theta=get_radiance(casename, tool)
                sigma=float(i)*1e-3*2.
                theoret, scale, th=guassian(sigma, theta)
                L*=scale
                plt.plot(theta[idx[j]], L[idx[j]], symbl[j], color='black')

            plt.plot(th, theoret, '-',color='black')

    plt.xticks(size=18)   
    plt.yticks(size=18)
    plt.axis(fontsize=18)
    plt.xlabel('$\\theta$',fontsize=18)
    plt.ylabel('Normalised Radiance',fontsize=18)
    plt.title(' Gaussian slope error',fontsize=18)

    plt.legend(fontsize=18)
    plt.savefig(open('./plots/radiance/A_1.2.png', 'w'), bbox_inches='tight')
    plt.close()

    # plot differences
    for i in xrange(1,4):
        plt.figure(i+1, dpi=1000)
        casename='1.2.%s'%i 
        sigma=float(i)*1e-3*2.
        theoret,scale, th=guassian(sigma, theta) 

        for j in xrange(len(tools)):
            tool=tools[j]
            L, theta=get_radiance(casename, tool)
            L*=scale
            diff=L-theoret
            plt.plot(theta, diff, color=colors[j], label=Tools[j])

        plt.xlabel('$\\theta$')
        plt.ylabel('Differences of the normalised Radiance from simulations \n  to the theoretical value')
        plt.title(' Gaussian slope error \n %s mrad'%i)

        plt.legend(loc='best')
        plt.savefig(open('./plots/radiance/A_%s_diff.png'%casename, 'w'), bbox_inches='tight')
        plt.close()

def plot_2_1():
    'pillbox sunshape'
    tools=N.array([ 'tonatiuh','soltrace','tracer', 'solstice', 'heliosim'])
    Tools=N.array([ 'Tonatiuh','SolTrace','Tracer', 'SOLSTICE', 'Heliosim'])
    colors=N.array(['b','orange','g', 'y','r'])
    symbl=N.array(['v','o', '^', 's', 'd'])

    idx={}
    idx[0]=N.arange(0,46,5).astype(int)
    idx[1]=N.arange(1,47,5).astype(int)
    idx[2]=N.arange(2,48,5).astype(int)
    idx[3]=N.arange(3,49,5).astype(int)
    idx[4]=N.arange(4,50,5).astype(int)

    plt.figure(1, dpi=1000)  

    casename='2.1'
 
    for j in xrange(len(tools)):
        tool=tools[j]
        L, theta=get_radiance(casename, tool)
        plt.plot(theta[idx[j]], L[idx[j]],  symbl[j], color='black', label=Tools[j])

    sigma=4.e-3
    theoret, th=pillbox(sigma) 
    plt.plot(th, theoret, '-',color='black', label='Theoretical')

    plt.xticks(size=18)   
    plt.yticks(size=18)
    plt.axis(fontsize=18)

    plt.xlabel('$\\theta$',fontsize=18)
    plt.ylabel('Normalised Radiance',fontsize=18)
    plt.title(' Pillbox sunshape',fontsize=18)

    plt.legend(fontsize=18)
    plt.savefig(open('./plots/radiance/A_2.1.png', 'w'), bbox_inches='tight')
    plt.close()

    # plot differences
    plt.figure(2, dpi=1000)
    casename='2.1' 
    sigma=4.e-3
    theoret=pillbox(sigma, theta, False) 

    for j in xrange(len(tools)):
        tool=tools[j]
        L, theta=get_radiance(casename, tool)
        diff=L-theoret
        plt.plot(theta, diff, color=colors[j], label=Tools[j])

    plt.xlabel('$\\theta$')
    plt.ylabel('Differences of the normalised Radiance from simulations \n  to the theoretical value')
    plt.title(' Pillbox Sunshape')

    plt.legend(loc='best')
    plt.savefig(open('./plots/radiance/A_%s_diff.png'%casename, 'w'), bbox_inches='tight')
    plt.close()

def plot_2_2():
    'Gaussian sunshape'
    tools=N.array([ 'soltrace','tracer', 'solstice', 'heliosim'])
    Tools=N.array([ 'SolTrace','Tracer', 'SOLSTICE', 'Heliosim'])
    colors=N.array(['orange','g', 'y','r'])
    symbl=N.array(['o', '^', 's', 'd'])

    idx={}
    idx[0]=N.arange(0,47,4).astype(int)
    idx[1]=N.arange(1,48,4).astype(int)
    idx[2]=N.arange(2,49,4).astype(int)
    idx[3]=N.arange(3,50,4).astype(int)

    plt.figure(1, dpi=1000)  

    casename='2.2'
 
    for j in xrange(len(tools)):
        tool=tools[j]
        L, theta=get_radiance(casename, tool)
        L/=2.
        plt.plot(theta[idx[j]], L[idx[j]],  symbl[j], color='black', label=Tools[j])

    sigma=4.e-3
    theoret,scale, th=guassian(sigma, theta) 
    plt.plot(th, theoret, '-',color='black', label='Theoretical')

    plt.xticks(size=18)   
    plt.yticks(size=18)
    plt.axis(fontsize=18)

    plt.xlabel('$\\theta$', fontsize=18)
    plt.ylabel('Normalised Radiance', fontsize=18)
    plt.title(' Gaussian sunshape', fontsize=18)

    plt.legend(fontsize=18)
    plt.savefig(open('./plots/radiance/A_2.2.png', 'w'), bbox_inches='tight')
    plt.close()

    # plot differences
    plt.figure(2, dpi=1000)
    casename='2.2' 
    sigma=4.e-3
    theoret, scale, th=guassian(sigma, theta) 

    for j in xrange(len(tools)):
        tool=tools[j]
        L, theta=get_radiance(casename, tool)
        L/=2.
        diff=L-theoret
        plt.plot(theta, diff, color=colors[j], label=Tools[j])

    plt.xlabel('$\\theta$')
    plt.ylabel('Differences of the normalised Radiance from simulations \n  to the theoretical value')
    plt.title(' Gaussian Sunshape')

    plt.legend(loc='best')
    plt.savefig(open('./plots/radiance/A_%s_diff.png'%casename, 'w'), bbox_inches='tight')
    plt.close()

    
def plot_2_3():
    'buie sunshape'
    # without CSR -- Tracer
    # tonatiuh CSR -- tonatiuh and solstice
    # CA CSR -- heliosim ?

    tools=N.array([ 'tonatiuh','soltrace','tracer', 'solstice', 'heliosim'])
    Tools=N.array([ 'Tonatiuh','SolTrace','Tracer', 'SOLSTICE', 'Heliosim'])
    colors=N.array(['b','orange','g', 'y','r'])
    symbl=N.array(['v','o', '^', 's', 'd'])

    idx={}
    idx[0]=N.arange(0,46,5).astype(int)
    idx[1]=N.arange(1,47,5).astype(int)
    idx[2]=N.arange(2,48,5).astype(int)
    idx[3]=N.arange(3,49,5).astype(int)
    idx[4]=N.arange(4,50,5).astype(int)

    plt.figure(1, dpi=1000)  

    for i in xrange(1,4):
        casename='2.3.%s'%i
     
        if i==1:
            CSR=float(i)*1e-2

            for j in xrange(len(tools)):
                tool=tools[j]
                L, theta=get_radiance(casename, tool)
                theoret, scale, th=buie(theta, CSR) 
                theoret2, scale2, th2=buie_tonatiuh(theta, CSR) 
                if tool==('tonatiuh' or 'solstice'):
                    L*=scale2
                else:
                    L*=scale
                plt.plot(theta[idx[j]], L[idx[j]], symbl[j], color='black', label=Tools[j])

            plt.plot(th, theoret, '-',color='black', label='w/o $\chi$ callibration')
            plt.plot(th2, theoret2, '-',color='red', label='with $\chi$ callibration')


        else:
            CSR=float(i)*1e-2

            for j in xrange(len(tools)):
                tool=tools[j]
                L, theta=get_radiance(casename, tool)
                theoret, scale, th=buie(theta, CSR) 
                theoret2, scale2, th2=buie_tonatiuh(theta, CSR) 
                if tool==('tracer' or 'heliosim'):
                    L*=scale
                else:
                    L*=scale2
                plt.plot(theta[idx[j]], L[idx[j]],   symbl[j], color='black')

            plt.plot(th, theoret, '-',color='black')
            plt.plot(th2, theoret2, '-',color='red')

    plt.xticks(size=18)   
    plt.yticks(size=18)
    plt.axis(fontsize=18)

    plt.xlabel('$\\theta$',fontsize=18)
    plt.ylabel('Normalised Radiance',fontsize=18)
    plt.title(' Buie sunshape',fontsize=18)
    ax=plt.gca()
    ax.set_yscale('log')
    plt.legend()
    plt.savefig(open('./plots/radiance/A_2.3.png', 'w'), bbox_inches='tight')
    plt.close()

    # plot differences
    for i in xrange(1,4):
        plt.figure(i+1, dpi=1000)
        casename='2.3.%s'%i 
        CSR=float(i)*1e-2
        theoret, scale, th=buie(theta, CSR) 

        for j in xrange(len(tools)):
            tool=tools[j]
            L, theta=get_radiance(casename, tool)
            L*=scale
            diff=L-theoret
            plt.plot(theta, diff, color=colors[j], label=Tools[j])

        plt.xlabel('$\\theta$')
        plt.ylabel('Differences of the normalised Radiance from simulations \n  to the theoretical value')
        plt.title(' Buie sunshape \n CSR %s'%(i*1.e-2))

        plt.legend(loc='best')
        plt.savefig(open('./plots/radiance/A_%s_diff.png'%casename, 'w'), bbox_inches='tight')
        plt.close()

def plot_buie():
    'buie sunshape'
    # without CSR -- Tracer
    # tonatiuh CSR -- tonatiuh and solstice
    # CA CSR -- heliosim ?

    tools=N.array([ 'tonatiuh','tracer', 'solstice', 'heliosim'])
    Tools=N.array([ 'Tonatiuh','Tracer', 'SOLSTICE', 'Heliosim'])
    #tools=N.array([ 'tracer'])
    #Tools=N.array([ 'Tracer'])

    colors=N.array(['b','g', 'y','r'])

    plt.figure(1, dpi=1000)  

    for i in xrange(1,4):
        casename='2.3.%s'%i
     
        if i==1:
            CSR=float(i)*1e-2

            for j in xrange(len(tools)):
                tool=tools[j]

                if tool==('heliosim'):
                    L, theta=get_radiance(casename, tool)
                    theoret, scale=buie(theta, CSR) 
                    L*=scale
                elif tool=='tracer':
                    L, theta=get_radiance(casename, tool)
                    theoret3, scale3=buie_CA(theta, CSR)  
                    L*=scale3
                else:
                    L, theta=get_radiance(casename, tool)
                    theoret2, scale2=buie_tonatiuh(theta, CSR)
                    L*=scale2
                plt.plot(theta, L, '.', color=colors[j], label=Tools[j])

            plt.plot(theta, theoret, '-',color='black', label='w/o $\chi$ callibration')
            plt.plot(theta, theoret2, '-',color='red', label='with $\chi$ callibration')
            plt.plot(theta, theoret3, '-',color='blue', label='with CA $\chi$ callibration')

        else:
            CSR=float(i)*1e-2

            for j in xrange(len(tools)):
                tool=tools[j]

                if tool==('heliosim'):
                    L, theta=get_radiance(casename, tool)
                    theoret, scale=buie(theta, CSR) 
                    L*=scale
                elif tool=='tracer':
                    L, theta=get_radiance(casename, tool)
                    theoret3, scale3=buie_CA(theta, CSR)  
                    L*=scale3
                else:
                    L, theta=get_radiance(casename, tool)
                    theoret2, scale2=buie_tonatiuh(theta, CSR)
                    L*=scale2
                plt.plot(theta, L, '.', color=colors[j])

            plt.plot(theta, theoret, '-',color='black')
            plt.plot(theta, theoret2, '-',color='red')
            plt.plot(theta, theoret3, '-',color='blue',)

    plt.xlabel('$\\theta$')
    plt.ylabel('Normalised Radiance')
    plt.title(' Buie sunshape')
    ax=plt.gca()
    ax.set_yscale('log')
    plt.legend()
    plt.savefig(open('./plots/radiance/A_2.3_buie.png', 'w'), bbox_inches='tight')
    plt.close()



if __name__=='__main__':
    #plot_1_1()
    #plot_1_2()
    #plot_2_1()    
    #plot_2_2() 
    plot_2_3() 
    #plot_buie()





        
    
