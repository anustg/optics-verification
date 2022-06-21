import numpy as N
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
#from repository.interpolation import interpolation
import matplotlib.cm as cm


def plot_fluxmaps(W,H,casename,reftool):
    '''   
    W - width of the receiver
    H - height of the receiver
    casename -case name
    '''  

    rd=casename[0]  
    
    if 'A_1.1' in casename:
        tools=N.array(['soltrace', 'tracer_test','solstice', 'heliosim'])

    elif casename=='A_2.2':
        tools=N.array(['soltrace', 'tracer_test','solstice', 'heliosim','solarpilot'])

    else:
        tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    spec_case=N.array(['_CSR_CA', '_CSR_tonatiuh'])


    fluxfile='../results/data/%s/%s/%s.csv' #(rd, tool, casename) 
  

    # the referece flux
    if casename=='A_2.2' or casename=='A_1.1.1' or casename=='A_1.1.2' or casename=='A_1.1.3':
        flux_res=N.loadtxt(fluxfile%(rd, reftool, casename), skiprows=1, delimiter=',')
    else:    
        flux_res=N.loadtxt(fluxfile%(rd, 'tonatiuh', casename), skiprows=1, delimiter=',')              
    flux_ref=flux_res[:,-1].astype(float)
    flux_ref=flux_ref.reshape(100,100)

    # get the range, max and min
    maxf=0.
    minf=10000.
    w_ratios=[]
    for t in xrange(len(tools)):
        if tools[t]!='solarpilot':
            flux_res=N.loadtxt(fluxfile%(rd, tools[t], casename), skiprows=1,dtype=float, delimiter=',')   

            m=100
            n=100
            fluxmap=flux_res[:,2].reshape(m,n)

            w_ratios.append(1)

            diff=(fluxmap-flux_ref)/(N.max(flux_ref))*100.

            if maxf<=N.max(diff):
                maxf=N.max(diff)
            if minf>=N.min(diff):
                minf=N.min(diff)
          

        else:
            flux_res=N.loadtxt(fluxfile%(rd, tools[t], casename), skiprows=1,dtype=float, delimiter=',')   

            m=100
            n=100
            fluxmap=flux_res[:,2].reshape(m,n)

            w_ratios.append(1)

            diff=(fluxmap-flux_ref)/(N.max(flux_ref))*100.
            rang_1=N.max(diff)
            rang_2=N.min(diff)


            

    if abs(minf)>abs(maxf):
        rang=round(abs(minf), 2)
    else:
        rang=round(abs(maxf),2)
  
    # plot the figure
    plt.figure(1,dpi=1000)
    fig,axes=plt.subplots(nrows=1,ncols=len(tools),figsize=(2.7*len(tools),1.7))
    gs=gridspec.GridSpec(1,len(tools),width_ratios=w_ratios,height_ratios=[1])

    for i in xrange(len(tools)):

        flux_res=N.loadtxt(fluxfile%(rd, tools[i], casename), skiprows=1, delimiter=',')

        m=100
        n=100
        X=flux_res[:,0].reshape(m,n)
        Y=flux_res[:,1].reshape(m,n)
        fluxmap=flux_res[:,2].reshape(m,n)

        ax=plt.subplot(gs[i])  
        if casename=='A_2.2' or casename=='A_1.1.1' or casename=='A_1.1.2' or casename=='A_1.1.3':
            if tools[i]==reftool:
                im1=ax.pcolormesh(X,Y,fluxmap, cmap=cm.hot)
                ref=N.max(fluxmap)
            else:
                diff=(fluxmap-flux_ref)/(N.max(flux_ref))*100.
                im=ax.pcolormesh(X,Y, diff,vmax=rang, vmin=-rang, cmap=cm.seismic)
        else:
            if tools[i]=='tonatiuh':   
                im1=ax.pcolormesh(X,Y,fluxmap, cmap=cm.hot)
                ref=N.max(fluxmap)
            else:
                diff=(fluxmap-flux_ref)/(N.max(flux_ref))*100.
                im=ax.pcolormesh(X,Y, diff,vmax=rang, vmin=-rang, cmap=cm.seismic)
                print 'tool:', tools[t]
                print 'max diff', N.max(abs(diff))
                print ''


        plt.xticks((-W/2., 0, W/2.), size=14)   
        plt.yticks((-H/2., 0, H/2.), size=14)
        plt.axis('equal',fontsize=14)
        
        tt=tools[i]
        if tt=='tracer_test':
            xlb='Tracer'
        elif tt=='heliosim':
            xlb='Heliosim'
        elif tt=='tonatiuh':
            xlb='Tonatiuh'
        elif tt=='solstice':
            xlb='SOLSTICE'
        elif tt=='solarpilot':
            xlb='SolarPILOT'
        else:
            xlb='SolTrace'
        plt.xlabel(xlb,fontsize=14)
        plt.subplots_adjust(wspace=0.5)

    fig.subplots_adjust(right=0.8)
    cbar_ax=fig.add_axes([0.82,0.2,0.02,0.6])
    clb=fig.colorbar(im,cax=cbar_ax,ticks=N.linspace(-rang,rang,3))
    clb.ax.set_title('  Difference %',x=1,y=1.2)
    clb.ax.tick_params(labelsize=14)


    fig.subplots_adjust(left=0.15)
    cbar_ax=fig.add_axes([0.1,0.2,0.02,0.6])
    clb=fig.colorbar(im1,cax=cbar_ax,ticks=N.linspace(0.,int(ref),3))
    clb.ax.set_title('  $kW/m^2$',x=0.07, y=1.1)
    clb.ax.tick_params(labelsize=14) 
    cbar_ax.yaxis.set_ticks_position('left')   


    savefile='./plots/fluxmap_diff'
    if not os.path.exists(savefile):
        os.makedirs(savefile)
    plt.savefig(open(savefile+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()



def collect_fluxmap():

    reftool='soltrace'


    # A   
    #A1
    for i in xrange(2):
        for j in xrange(3):
            W=8.
            H=8.
            casename='A_1.%s.%s'%(i+1,j+1)
            print ''
            print casename
            plot_fluxmaps(W,H,casename,reftool)



    #A2.3
    for j in xrange(3):
        W=8.
        H=8.
        casename='A_2.3.%s'%(j+1)
        print ''
        print casename
        plot_fluxmaps(W,H,casename,reftool)


    #A2 and A3
    for j in xrange(2,4):
        for i in xrange(2):
            W=8.
            H=8.
            casename='A_%s.%s'%(j,i+1)
            print ''
            print casename
            plot_fluxmaps(W,H,casename,reftool)

    # B
    for h in xrange(2):
        for i in xrange(2):
            for j in xrange(4):
                W=8.
                H=6.
                casename='B_%s.%s.%s'%(h+1,i+1,j+1)
                print ''
                print casename
                plot_fluxmaps(W,H,casename,reftool)


    # C
    for i in xrange(2):
        for j in xrange(2):
            W=8.
            H=6.
            casename='C_%s.%s'%(i+1,j+1)
            print ''
            print casename
            plot_fluxmaps(W,H,casename,reftool)
 

def plot_buie(W,H,casename):
    '''   
    W - width of the receiver
    H - height of the receiver
    casename -case name
    '''  

    rd='A'
    
    tools=N.array(['tonatiuh',  'tracer_test', 'solstice', 'heliosim'])
    spec_case=N.array(['_CSR_CA', '_CSR_tonatiuh'])


    fluxfile='../results/data/%s/%s/%s.csv' #(rd, tool, casename) 
  
    w_ratios=[1,1]
    maxf=0.
    minf=10000.

    # the referece flux
    flux_res=N.loadtxt(fluxfile%(rd, 'tonatiuh', casename), skiprows=1, delimiter=',')              
    flux_ref=flux_res[:,-1].astype(float)
    flux_ref=flux_ref.reshape(100,100)

    # get the range, max and min
    maxf=0.
    minf=10000.
    for t in xrange(len(tools)):
        flux_res=N.loadtxt(fluxfile%(rd, tools[t], casename), skiprows=1,dtype=float, delimiter=',')   

        m=100
        n=100
        fluxmap=flux_res[:,2].reshape(m,n)

        w_ratios.append(1)

        diff=fluxmap-flux_ref
        if maxf<=N.max(diff):
            maxf=N.max(diff)
        if minf>=N.min(diff):
            minf=N.min(diff) 

    if abs(minf)>abs(maxf):
        rang=round(abs(minf), 2)
    else:
        rang=round(abs(maxf),2)


    plt.figure(1,dpi=1000)
    fig,axes=plt.subplots(nrows=1,ncols=len(tools)+2,figsize=(2.7*(len(tools)+2),1.7))
    gs=gridspec.GridSpec(1,len(tools)+2,width_ratios=w_ratios,height_ratios=[1])

    for i in xrange(len(tools)):

        flux_res=N.loadtxt(fluxfile%(rd, tools[i], casename), skiprows=1, delimiter=',')
 
        m=100
        n=100

        X=flux_res[:,0].reshape(m,n)
        Y=flux_res[:,1].reshape(m,n)
        fluxmap=flux_res[:,2].reshape(m,n) 

        ax=plt.subplot(gs[i])   
        if tools[i]=='tonatiuh':   
            im1=ax.pcolormesh(X,Y,fluxmap)
            ref=N.max(fluxmap)
        else:
            im=ax.pcolormesh(X,Y, fluxmap-flux_ref,vmax=rang, vmin=-rang)

        plt.xticks((-W/2., 0, W/2.), size=14)   
        plt.yticks((-H/2., 0, H/2.), size=14)
        plt.axis('equal',fontsize=14)
        
        tt=tools[i]
        if tt=='tracer_test':
            xlb='Tracer'
        elif tt=='heliosim':
            xlb='Heliosim'
        elif tt=='tonatiuh':
            xlb='Tonatiuh'
        elif tt=='solstice':
            xlb='SOLSTICE'
        elif tt=='solarpilot':
            xlb='SolarPILOT'
        else:
            xlb='SolTrace'
        plt.xlabel(xlb,fontsize=14)
        plt.subplots_adjust(wspace=0.5)

    for i in xrange(len(spec_case)):
        flux_res=N.loadtxt(fluxfile%(rd, 'tracer_test', casename+spec_case[i]), skiprows=1, delimiter=',')
 
        m=100
        n=100

        X=flux_res[:,0].reshape(m,n)
        Y=flux_res[:,1].reshape(m,n)
        fluxmap=flux_res[:,2].reshape(m,n) 

        ax=plt.subplot(gs[4+i])     
        im=ax.pcolormesh(X,Y,fluxmap-flux_ref,vmin=-rang,vmax=rang)
        plt.xticks((-W/2., 0, W/2.), size=14)   
        plt.yticks((-H/2., 0, H/2.), size=14)
        plt.axis('equal',fontsize=14)
        
        xlb='Tracer'+spec_case[i]
        plt.xlabel(xlb,fontsize=14)
        plt.subplots_adjust(wspace=0.5)
            

    fig.subplots_adjust(right=0.8)
    cbar_ax=fig.add_axes([0.82,0.2,0.02,0.6])
    clb=fig.colorbar(im,cax=cbar_ax,ticks=N.linspace(-rang,rang,3))
    clb.ax.set_title('  $kW/m^2$',x=0.85,y=1.06)
    clb.ax.tick_params(labelsize=14)

    fig.subplots_adjust(left=0.15)
    cbar_ax=fig.add_axes([0.1,0.2,0.02,0.6])
    clb=fig.colorbar(im1,cax=cbar_ax,ticks=N.linspace(0.,int(ref),3))
    clb.ax.set_title('  $kW/m^2$',x=0.07, y=1.06)
    clb.ax.tick_params(labelsize=14) 
    cbar_ax.yaxis.set_ticks_position('left')   


    savefile='./plots/fluxmap_diff'
    if not os.path.exists(savefile):
        os.makedirs(savefile)
    plt.savefig(open(savefile+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()

if __name__=='__main__':
    collect_fluxmap()
    #for j in xrange(3):
    #    W=8.
    #    H=8.
    #    casename='A_2.3.%s'%(j+1)
    #    print ''
    #    print casename
    #    plot_buie(W,H,casename)



