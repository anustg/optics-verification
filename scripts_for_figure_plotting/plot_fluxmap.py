import numpy as N
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
#from repository.interpolation import interpolation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm

def plot_A11(W, H):

    for k in xrange(1,4):
        casename='A_1.1.%s'%k

        tools=N.array(['tonatiuh','soltrace','tracer_test', 'solstice', 'heliosim', 'solarpilot'])
        fluxfile='../results/data/A/%s/%s.csv' #(rd, tool, casename) 
      
        w_ratios=[1,1,1,1,1,1]
        maxf=0.
        minf=10000.
        for t in xrange(len(tools)-1):
            if tools[t]!='tonatiuh':
                flux_res=N.loadtxt(fluxfile%(tools[t], casename), skiprows=1, delimiter=',')
                fluxmap=flux_res[:,-1].astype(float)
                if maxf<N.max(fluxmap):
                    maxf=N.max(fluxmap)
                if minf>N.min(fluxmap):
                    minf=N.min(fluxmap)

        plt.figure(1,dpi=1000)
        fig,axes=plt.subplots(nrows=1,ncols=6,figsize=(2.7*6,1.7))
        gs=gridspec.GridSpec(1,len(tools),width_ratios=w_ratios,height_ratios=[1])

        for i in xrange(len(tools)):

            if tools[i]=='solarpilot':
                ax=plt.subplot(gs[i])    
                plt.annotate('N/A', xy=(-0.009,-0.004)) 
                plt.xticks([], []) 
                plt.yticks([], [])
                plt.axis('equal',fontsize=14)

            else:

                flux_res=N.loadtxt(fluxfile%( tools[i], casename), skiprows=1, delimiter=',')
                m=100
                n=100
                X=flux_res[:,0].reshape(m,n)
                Y=flux_res[:,1].reshape(m,n)
                fluxmap=flux_res[:,2].reshape(m,n) 

                ax=plt.subplot(gs[i])     
                im=ax.pcolormesh(X,Y,fluxmap,cmap=cm.hot, vmin=minf,vmax=maxf)
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
                xlb='\n SolarPILOT'
            else:
                xlb='SolTrace'
            plt.xlabel(xlb,fontsize=14)
            plt.subplots_adjust(wspace=0.5)

        fig.subplots_adjust(right=0.8)
        cbar_ax=fig.add_axes([0.82,0.2,0.02,0.6])
        clb=fig.colorbar(im,cax=cbar_ax,ticks=N.linspace(int(minf),int(maxf),3))
        clb.ax.set_title('  $kW/m^2$',y=1.05)
        clb.ax.tick_params(labelsize=14)
        savefile='./plots/fluxmap'
        if not os.path.exists(savefile):
            os.makedirs(savefile)
        plt.savefig(open(savefile+'/%s.png'%casename,'w'), bbox_inches='tight')
        plt.close()

def plot_A22(W, H):

    casename='A_2.2'

    tools=N.array(['tonatiuh','soltrace','tracer_test', 'solstice', 'heliosim', 'solarpilot'])
    fluxfile='../results/data/A/%s/%s.csv' #(rd, tool, casename) 
  
    w_ratios=[1,1,1,1,1,1]
    maxf=0.
    minf=10000.
    for t in xrange(len(tools)-1):
        if tools[t]!='tonatiuh'and tools[t]!='solarpilot':
            flux_res=N.loadtxt(fluxfile%(tools[t], casename), skiprows=1, delimiter=',')
            fluxmap=flux_res[:,-1].astype(float)
            if maxf<N.max(fluxmap):
                maxf=N.max(fluxmap)
            if minf>N.min(fluxmap):
                minf=N.min(fluxmap)

    plt.figure(1,dpi=1000)
    fig,axes=plt.subplots(nrows=1,ncols=6,figsize=(2.7*6,1.7))
    gs=gridspec.GridSpec(1,len(tools),width_ratios=w_ratios,height_ratios=[1])

    for i in xrange(len(tools)):

        if tools[i]=='tonatiuh':
            ax=plt.subplot(gs[i])    
            plt.annotate('N/A', xy=(-0.009,-0.004)) 
            plt.xticks([], []) 
            plt.yticks([], [])
            plt.axis('equal',fontsize=14)

        else:

            flux_res=N.loadtxt(fluxfile%( tools[i], casename), skiprows=1, delimiter=',')
            m=100
            n=100
            X=flux_res[:,0].reshape(m,n)
            Y=flux_res[:,1].reshape(m,n)
            fluxmap=flux_res[:,2].reshape(m,n) 

            ax=plt.subplot(gs[i])     
            im=ax.pcolormesh(X,Y,fluxmap,cmap=cm.hot,vmin=minf,vmax=maxf)
            plt.xticks((-W/2., 0, W/2.), size=14)   
            plt.yticks((-H/2., 0, H/2.), size=14)
            plt.axis('equal',fontsize=14)
        
        tt=tools[i]
        if tt=='tracer_test':
            xlb='Tracer'
        elif tt=='heliosim':
            xlb='Heliosim'
        elif tt=='tonatiuh':
            xlb='\n Tonatiuh'
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
    clb=fig.colorbar(im,cax=cbar_ax,ticks=N.linspace(int(minf),int(maxf),3))
    clb.ax.set_title('  $kW/m^2$',y=1.05)
    clb.ax.tick_params(labelsize=14)
    savefile='./plots/fluxmap'
    if not os.path.exists(savefile):
        os.makedirs(savefile)
    plt.savefig(open(savefile+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()



def plot_fluxmaps(W,H,casename):
    '''   
    W - width of the receiver
    H - height of the receiver
    casename -case name
    '''  

    rd=casename[0]  
    tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])


    fluxfile='../results/data/%s/%s/%s.csv' #(rd, tool, casename) 
  
    w_ratios=[]
    maxf=0.
    minf=10000.
    for t in xrange(len(tools)):

        if tools[t]!='tonatiuh' and tools[t]!='solarpilot':
            flux_res=N.loadtxt(fluxfile%(rd, tools[t], casename), skiprows=1, delimiter=',')
            
            fluxmap=flux_res[:,-1].astype(float)

            if maxf<N.max(fluxmap):
                maxf=N.max(fluxmap)
            if minf>N.min(fluxmap):
                minf=N.min(fluxmap)
        w_ratios.append(1)

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
        im=ax.pcolormesh(X,Y,fluxmap,cmap=cm.hot,vmin=minf,vmax=maxf)
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
    clb=fig.colorbar(im,cax=cbar_ax,ticks=N.linspace(int(minf),int(maxf),3))
    clb.ax.set_title('  $kW/m^2$',y=1.05)
    clb.ax.tick_params(labelsize=14)
    savefile='./plots/fluxmap'
    if not os.path.exists(savefile):
        os.makedirs(savefile)
    plt.savefig(open(savefile+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()



def collect_fluxmap(case):

    if case=='A_1.2':
        for j in xrange(3):
            W=8.
            H=8.
            casename='A_1.2.%s'%(j+1)
            print ''
            print casename
            plot_fluxmaps(W,H,casename)

    elif case=='A_2.1':
        W=8.
        H=8.
        casename='A_2.1'
        print ''
        print casename
        plot_fluxmaps(W,H,casename)
      

    elif case=='A_2.3':
        for j in xrange(3):
            W=8.
            H=8.
            casename='A_2.3.%s'%(j+1)
            print ''
            print casename
            plot_fluxmaps(W,H,casename)

    elif case=='A_3':
        for i in xrange(2):
            W=8.
            H=8.
            casename='A_3.%s'%(i+1)
            print ''
            print casename
            plot_fluxmaps(W,H,casename)

    elif case=='B':
        # B
        for h in xrange(2):
            for i in xrange(2):
                for j in xrange(4):
                    W=8.
                    H=6.
                    casename='B_%s.%s.%s'%(h+1,i+1,j+1)
                    print ''
                    print casename
                    plot_fluxmaps(W,H,casename)

    elif case=='C':
        # C
        for i in xrange(2):
            for j in xrange(2):
                W=8.
                H=6.
                casename='C_%s.%s'%(i+1,j+1)
                print ''
                print casename
                plot_fluxmaps(W,H,casename)

def plot_tonatiuh_A1_1():
    '''   
    W - width of the receiver
    H - height of the receiver
    casename -case name
    '''  
    W=8.
    H=8.

    fluxfile='../results/data/A/tonatiuh/A_1.1.%s.csv' 
  
    plt.figure(1, dpi=1000)

    fig,axes=plt.subplots(nrows=1,ncols=3,figsize=(2.7*3.,1.7))
    gs=gridspec.GridSpec(1,3,width_ratios=[1,1,1],height_ratios=[1])

    for i in xrange(1,4):

        flux_res=N.loadtxt(fluxfile%(i), skiprows=1, delimiter=',')       
        fluxmap=flux_res[:,-1].astype(float)
 
        m=100
        n=100

        X=flux_res[:,0].reshape(m,n)
        Y=flux_res[:,1].reshape(m,n)
        fluxmap=flux_res[:,2].reshape(m,n) 

        ax=plt.subplot(gs[i-1])     
        im=ax.pcolormesh(X,Y,fluxmap)
        plt.xticks((-W/2., 0., W/2.), size=14)   
        plt.yticks((-H/2., 0., H/2.), size=14)
        #ax.set_ylim([-H/2., H/2.])
        #ax.set_xlim([-W/2., W/2.])    
        #plt.axis('equal',fontsize=14)
        plt.xlabel('%s mrad'%(int(i)),fontsize=14)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        clb=fig.colorbar(im, cax=cax, orientation='vertical')
        #clb.ax.set_title(' $kW/m^2$',y=1.01)
        #clb.ax.tick_params(labelsize=10)

        plt.subplots_adjust(wspace=0.5)

    plt.setp(axes, xticks=(-W/2., 0., W/2.), yticks=[-H/2., 0., H/2.])


    #fig.subplots_adjust(right=0.8)
    #cbar_ax=fig.add_axes([0.82,0.2,0.02,0.6])
    #clb=fig.colorbar(im,cax=cbar_ax,ticks=N.linspace(int(minf),int(maxf),3))

    st=plt.suptitle('Pillbox Slope Error \n Tonatiuh', fontsize=14)
    #figure.tight_layout()
    st.set_y(1.3)
    fig.subplots_adjust(top=1.)
    savefile='./plots/fluxmap'
    if not os.path.exists(savefile):
        os.makedirs(savefile)
    plt.savefig(open(savefile+'/tonatiuh_A1.1.png','w'), bbox_inches='tight')
    plt.close()

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
    for t in xrange(len(tools)):

        if tools[t]!='tonatiuh':
            flux_res=N.loadtxt(fluxfile%(rd, tools[t], casename), skiprows=1, delimiter=',')
            
            fluxmap=flux_res[:,-1].astype(float)

            if maxf<N.max(fluxmap):
                maxf=N.max(fluxmap)
            if minf>N.min(fluxmap):
                minf=N.min(fluxmap)
        w_ratios.append(1)

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
        im=ax.pcolormesh(X,Y,fluxmap,vmin=minf,vmax=maxf)
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
        im=ax.pcolormesh(X,Y,fluxmap,vmin=minf,vmax=maxf)
        plt.xticks((-W/2., 0, W/2.), size=14)   
        plt.yticks((-H/2., 0, H/2.), size=14)
        plt.axis('equal',fontsize=14)
        
        xlb='Tracer'+spec_case[i]
        plt.xlabel(xlb,fontsize=14)
        plt.subplots_adjust(wspace=0.5)
            

    fig.subplots_adjust(right=0.8)
    cbar_ax=fig.add_axes([0.82,0.2,0.02,0.6])
    clb=fig.colorbar(im,cax=cbar_ax,ticks=N.linspace(int(minf),int(maxf),3))
    clb.ax.set_title('  $kW/m^2$',y=1.05)
    clb.ax.tick_params(labelsize=14)
    savefile='./plots/fluxmap'
    if not os.path.exists(savefile):
        os.makedirs(savefile)
    plt.savefig(open(savefile+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()

if __name__=='__main__':
    plot_A11(W=8, H=8)
    plot_A22(W=8, H=8)
    collect_fluxmap('A_1.2')
    collect_fluxmap('A_2.1')
    collect_fluxmap('A_2.3')
    collect_fluxmap('A_3')
    collect_fluxmap('B')
    collect_fluxmap('C')




