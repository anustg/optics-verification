import numpy as N
import os
import matplotlib.pyplot as plt
from matplotlib import ticker


def plot_energy(casename, title):

    '''
    casename:
    A_1 - slope error
    A_2 - sunshape
    A_2.3 - Buie sunshape
    A_3 -sunshape + slope error
    B_1 - pillbox sunshape
    B_2 - buie sunshape
    '''

    res_dir='../results/data/%s/%s/summary.csv' #tool
    tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    spec_case=N.array(['_CSR_CA', '_CSR_tonatiuh'])


    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(casename[0],tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):

            if casename=='A_2':
                if ('A_2.1' in res[j,0]) or ('A_2.2' in res[j,0]):
                    RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

            else:
                if casename in res[j,0]:
                    RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)

 
    ymax=N.max(RESULTS[tools[i]][:,1].astype(float))*1.05
    ymin=N.min(RESULTS[tools[i]][:,2].astype(float))*0.8

    #============
    if casename=='A_1':
        NAME=N.array(['A 1.1.1 \n1 mrad \nPillbox ', 'A 1.1.2 \n2 mrad \nPillbox ' , 'A 1.1.3 \n3 mrad \nPillbox ', 'A 1.2.1 \n1 mrad \nNormal ', 'A 1.2.2 \n2 mrad \nNormal ', 'A 1.2.3 \n3 mrad \nNormal ' ])
    elif casename=='A_2':
        NAME=N.array(['A 2.1 \n4 mrad \nPillbox ', 'A 2.2 \n4 mrad \nGaussian '])
    elif casename=='A_3':
        NAME=N.array(['A 3.1 \nNormal slope error \nPillbox sunshape ', 'A 3.2 \n Normal slope error \nBuie sunshape ' ])
    else:
        #NAME=RESULTS[tools[0]][:,0]
        NAME=N.array([])
        for i in xrange(2):
            for j in xrange(4):
                NAME=N.append(NAME, '%s.%s.%s'%(title, i+1, j+1))


  
    n=len(RESULTS[tools[0]])
    width=0.1
    ind=N.arange(0.1,(n)*0.8,0.8)

    abs_color=N.array(['b','orange','g','y','r','m'])
    spil_color=N.array(['0.82','0.82','0.82','0.82','0.82','0.82'])
    leg=[]


    plt.figure(1,dpi=1000)
    for i in xrange(len(tools)):
        
        #print tools[i]
        #print RESULTS[tools[i]][:,2].astype(float)


        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,color=abs_color[i],label=tools[i], align='edge',edgecolor=['black']*n)
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,2].astype(float),color=spil_color[i], align='edge',edgecolor=['black']*n)
        plt.ylim(ymin, ymax)  
        plt.xlim([0,n*0.8])
        leg.append(p1)
    leg.append(p2)


    plt.ylabel('Energy (kW)')
    plt.title(title)
    plt.xticks(ind+width*2.5,NAME)

    plt.legend(leg,['$Q_{abs}$ (Tonatiuh)','$Q_{abs}$ (SolTrace)','$Q_{abs}$ (Tracer)','$Q_{abs}$ (SOLSTICE)','$Q_{abs}$ (Heliosim)','$Q_{abs}$ (SolarPILOT)','$Q_{spil}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()


def plot_energy_buie(casename, title):

    '''
    casename:
    A_2.3 - Buie sunshape
    '''

    res_dir='../results/data/A/%s/summary.csv' #tool
    tools=N.array(['tonatiuh','solstice','tracer_test', 'heliosim',  ])


    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):
            if casename in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)

    total=RESULTS[tools[i]][:,2].astype(float)
    ymax=N.max(total)*1.05
    ymin=N.min(total)*0.8

    #============
    NAME=RESULTS[tools[0]][:,0]
    
  
    n=len(RESULTS[tools[0]])
    ind=N.arange(n)
    width=0.12
    abs_color=N.array(['b','y','g','r','m','darkorange'])
    spil_color=N.array(['0.8','0.82','0.82','0.82','0.84','0.84'])
    leg=[]

    plt.figure(1,dpi=1000)
    for i in xrange(len(tools)):
        if tools[i]=='tracer_test':
            CA=N.array([])
            Ton=N.array([])
            non=N.array([])
            for res in RESULTS[tools[i]]:
                if 'CA' in res[0]:
                    CA=N.append(CA, res)
                elif 'tonatiuh' in res[0]:
                    Ton=N.append(Ton, res)
                else:
                    non=N.append(non, res)

            CA=CA.reshape(len(CA)/11, 11)
            Ton=Ton.reshape(len(Ton)/11, 11)
            non=non.reshape(len(non)/11, 11)

            p1=plt.bar(ind+width*(i),Ton[:,2].astype(float),width,color=abs_color[i],label=tools[i]+'_CSR_Tonatiuh', align='center',edgecolor=['black']*n)
            p2=plt.bar(ind+width*(i),Ton[:,3].astype(float),width,bottom=Ton[:,2].astype(float),color=spil_color[i], align='center',edgecolor=['black']*n)
            leg.append(p1)

            p1=plt.bar(ind+width*(i+1),CA[:,2].astype(float),width,color=abs_color[i],label=tools[i]+'_CSR_CA', align='center',edgecolor=['black']*n)
            p2=plt.bar(ind+width*(i+1),CA[:,3].astype(float),width,bottom=CA[:,2].astype(float),color=spil_color[i+1], align='center',edgecolor=['black']*n)
            leg.append(p1)


            p1=plt.bar(ind+width*(i+2),non[:,2].astype(float),width,color=abs_color[i],label=tools[i], align='center',edgecolor=['black']*n)
            p2=plt.bar(ind+width*(i+2),non[:,3].astype(float),width,bottom=non[:,2].astype(float),color=spil_color[i+2], align='center',edgecolor=['black']*n)



        else:
            if i==3:
                j=i+2
            else:
                j=i

            p1=plt.bar(ind+width*j,RESULTS[tools[i]][:,2].astype(float),width,color=abs_color[i],label=tools[i], align='center',edgecolor=['black']*n)
            p2=plt.bar(ind+width*j,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,2].astype(float),color=spil_color[i], align='center',edgecolor=['black']*n)

        leg.append(p1)
    leg.append(p2)
    plt.ylim(ymin, ymax)  

    plt.ylabel('Energy (kW)')
    plt.title(title)
    plt.xticks(ind+width*2.5,NAME)
    plt.legend(leg,['$Q_{abs}$ (Tonatiuh)','$Q_{abs}$ (SOLSTICE)','$Q_{abs}$ (Tracer_CSR)','$Q_{abs}$ (Tracer_CSR2)','$Q_{abs}$ (Tracer)','$Q_{abs}$ (Heliosim)','$Q_{spil}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)



    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()


def plot_energy_C(title):

    res_dir='../results/data/C/%s/summary.csv' #tool

    tools=N.array(['tonatiuh','soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])

    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):           
            if 'C_1.1' in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])
        for j in xrange(len(res)):      
            if 'C_2.1' in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])
        for j in xrange(len(res)):      
            if 'C_1.2' in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])
        for j in xrange(len(res)):      
            if 'C_2.2' in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)
        #print ''
        #print RESULTS[tools[i]]
          
    n=len(RESULTS[tools[0]])
    ind=N.arange(n)
    width=0.11
    abs_color=N.array(['b','orange','g','y','r','m'])

    NAME=N.array(['C 1.1 \nSolar Noon \n(Pillbox)', 'C 2.1 \nSolar Noon \n(Buie)','C 1.2 \n Morning \n(Pillbox)', 'C 2.2 \n Morning \n(Buie)'   ])    
    leg=[]

    plt.figure(1,dpi=1000)
    for i in xrange(len(tools)):
        #abs
        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,7].astype(float),width,color=abs_color[i],label=tools[i], align='center',edgecolor=['black']*n)
        # refl
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,6].astype(float),width,bottom=RESULTS[tools[i]][:,7].astype(float),color='cyan', align='center',edgecolor=['black']*n)
        # spil
        p3=plt.bar(ind+width*i,RESULTS[tools[i]][:,5].astype(float),width,bottom=RESULTS[tools[i]][:,7].astype(float)+RESULTS[tools[i]][:,6].astype(float),color='grey', align='center',edgecolor=['black']*n)
        # field absorb
        p4=plt.bar(ind+width*i,RESULTS[tools[i]][:,4].astype(float),width,bottom=RESULTS[tools[i]][:,7].astype(float)+RESULTS[tools[i]][:,6].astype(float)+RESULTS[tools[i]][:,5].astype(float),color='lightpink', align='center',edgecolor=['black']*n)        
        # block
        p5=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,7].astype(float)+RESULTS[tools[i]][:,6].astype(float)+RESULTS[tools[i]][:,5].astype(float)+RESULTS[tools[i]][:,4].astype(float),color='limegreen', align='center',edgecolor=['black']*n)
        # shade+cos
        p6=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,bottom=RESULTS[tools[i]][:,7].astype(float)+RESULTS[tools[i]][:,6].astype(float)+RESULTS[tools[i]][:,5].astype(float)+RESULTS[tools[i]][:,4].astype(float)+RESULTS[tools[i]][:,3].astype(float),color='yellow', align='center',edgecolor=['black']*n)


        #plt.ylim(60.,101.)  
        leg.append(p1)
    leg.append(p2)
    leg.append(p3)
    leg.append(p4)
    leg.append(p5)
    leg.append(p6)

    plt.ylabel('Energy (kW)')
    plt.title(title)
    plt.xticks(ind+width*2., NAME)

    plt.legend(leg,['$Q_{abs}$ (Tonatiuh)','$Q_{abs}$ (SolTrace)','$Q_{abs}$ (Tracer)','$Q_{abs}$ (SOLSTICE)','$Q_{abs}$ (Heliosim)','$Q_{abs}$ (SolarPILOT)','$Q_{refl}$','$Q_{spil}$', '$Q_{field,abs}$', '$Q_{block}$', '$Q_{shade&cos}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/C.png','w'), bbox_inches='tight')
    plt.close()


def plot_A11(title):

    ft=16.

    casename='A_1.1'
    res_dir='../results/data/A/%s/summary.csv' #tool
    tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim'])
    lab=['Tonatiuh', 'SolTrace','Tracer', 'SOLSTICE','Heliosim']

    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):

            if casename in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)


    qmax=N.max(RESULTS[tools[i]][:,1].astype(float))
    qmin=N.min(RESULTS[tools[i]][:,2].astype(float))
    ymax=qmax*1.05
    ymin=qmin*0.8

    #============

    NAME=N.array([' 1 mrad \n A 1.1.1', ' 2 mrad \n A 1.1.2  ' , ' 3 mrad \n A 1.1.3 ' ])

    n=len(RESULTS[tools[0]])
    width=0.12
    ind=N.arange(0.1,(n)*0.8,0.8)

    leg=[]

    plt.figure(1,dpi=1000)
    fig, ax=plt.subplots()
    for i in xrange(len(tools)):
        
        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,color='0.9',label=tools[i], align='edge',edgecolor=['black']*n)
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,2].astype(float),color='0.2', align='edge',edgecolor=['black']*n)
        plt.ylim(ymin, ymax)  
        plt.xlim([0,n*0.8])

        for j in xrange(n):   
            ax.text(ind[j]+width*i+width/2.*1.22, ymin+(qmax-ymin)/(2.*n)*(i+1), lab[i],horizontalalignment='center',
        verticalalignment='center',rotation=90,fontsize=12,color='red')

    
    leg.append(p1)
    leg.append(p2)

    plt.ylabel('Energy (kW)', fontsize=ft)
    plt.title(title, fontsize=ft)
    plt.xticks(ind+width*2.5,NAME,fontsize=ft)
    plt.yticks(fontsize=ft)

    plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=ft)


    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()

def plot_A12(title):

    ft=16.

    casename='A_1.2'
    res_dir='../results/data/A/%s/summary.csv' #tool
    tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    lab=['Tonatiuh', 'SolTrace','Tracer', 'SOLSTICE','Heliosim','SolarPILOT']

    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):

            if casename in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)


    qmax=N.max(RESULTS[tools[i]][:,1].astype(float))
    qmin=N.min(RESULTS[tools[i]][:,2].astype(float))
    ymax=qmax*1.05
    ymin=qmin*0.8

    #============

    NAME=N.array(['  1 mrad \n A 1.2.1', ' 2 mrad \n A 1.2.2  ' , ' 3 mrad \n A 1.2.3 ' ])

    n=len(RESULTS[tools[0]])
    width=0.11
    ind=N.arange(0.08,(n)*0.8,0.8)

    leg=[]

    plt.figure(1,dpi=1000)
    fig, ax=plt.subplots()
    for i in xrange(len(tools)):
        
        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,color='0.9',label=tools[i], align='edge',edgecolor=['black']*n)
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,2].astype(float),color='0.2', align='edge',edgecolor=['black']*n)
        plt.ylim(ymin, ymax)  
        plt.xlim([0,n*0.8])

        for j in xrange(n):   
            ax.text(ind[j]+width*i+width/2.*1.22, ymin+(qmax-ymin)/(3.*n)*(i+1), lab[i],horizontalalignment='center',
        verticalalignment='center',rotation=90,fontsize=12,color='red')

    
    leg.append(p1)
    leg.append(p2)

    plt.ylabel('Energy (kW)', fontsize=ft)
    plt.title(title, fontsize=ft)
    plt.xticks(ind+width*3.,NAME,fontsize=ft)
    plt.yticks(fontsize=ft)

    plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=ft)


    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()

def plot_A2(title):

    ft=16.

    casename='A_2'
    res_dir='../results/data/A/%s/summary.csv' #tool
    tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    lab=['Tonatiuh', 'SolTrace','Tracer', 'SOLSTICE','Heliosim','SolarPILOT']

    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):
            if 'A_2.3' not in res[j,0]:
                if casename in res[j,0]:
                    RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)


    qmax=N.max(RESULTS[tools[i]][:,1].astype(float))
    qmin=N.min(RESULTS[tools[i]][:,2].astype(float))
    ymax=qmax*1.05
    ymin=qmin*0.8

    #============

    NAME=N.array(['Pillbox \n A 2.1', ' Normal \n A 2.2  ' ])

    n=len(RESULTS[tools[0]])
    width=0.12
    ind=N.arange(0.08,(n)*0.7,0.7)

    leg=[]

    plt.figure(1,dpi=1000)
    fig, ax=plt.subplots()
    for i in xrange(len(tools)):
        
        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,color='0.9',label=tools[i], align='edge', edgecolor=['black']*n)
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,2].astype(float),color='0.2', align='edge',edgecolor=['black']*n)
        plt.ylim(ymin, ymax)  
        plt.xlim([0,n*0.8])

        if tools[i]=='tonatiuh':
            for j in xrange(1):   
                ax.text(ind[j]+width*i+width/2.*1.22, ymin+(qmax-ymin)/(3.*n)*(i+1), lab[i],horizontalalignment='center',
            verticalalignment='center',rotation=90,fontsize=12,color='red')

        else:
            for j in xrange(n):   
                ax.text(ind[j]+width*i+width/2.*1.2, ymin+(qmax-ymin)/(4.*n)*(i+1), lab[i],horizontalalignment='center',
            verticalalignment='center',rotation=90,fontsize=12,color='red')

    
    leg.append(p1)
    leg.append(p2)

    plt.ylabel('Energy (kW)', fontsize=ft)
    plt.title(title, fontsize=ft)
    plt.xticks([ind[0]+width*3., ind[1]+width*3.5], NAME,fontsize=ft)
    plt.yticks(fontsize=ft)

    plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=ft)


    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()

def plot_A23(title):

    ft=16.

    casename='A_2.3'
    res_dir='../results/data/A/%s/summary.csv' #tool
    tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    lab=['Tonatiuh', 'SolTrace','Tracer', 'SOLSTICE','Heliosim','SolarPILOT']

    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):

            if casename in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)


    qmax=N.max(RESULTS[tools[i]][:,1].astype(float))
    qmin=N.min(RESULTS[tools[i]][:,2].astype(float))
    ymax=qmax*1.05
    ymin=qmin*0.8

    #============

    NAME=N.array(['  CSR 0.01 \n A 2.3.1', ' CSR 0.02 \n A 2.3.2  ' , ' CSR 0.03 \n A 2.3.3 ' ])

    n=len(RESULTS[tools[0]])
    width=0.11
    ind=N.arange(0.08,(n)*0.8,0.8)

    leg=[]

    plt.figure(1,dpi=1000)
    fig, ax=plt.subplots()
    for i in xrange(len(tools)):
        
        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,color='0.9',label=tools[i], align='edge',edgecolor=['black']*n)
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,2].astype(float),color='0.2', align='edge',edgecolor=['black']*n)
        plt.ylim(ymin, ymax)  
        plt.xlim([0,n*0.8])

        for j in xrange(n):   
            ax.text(ind[j]+width*i+width/2.*1.22, ymin+(qmax-ymin)/(2.5*n)*(i+1), lab[i],horizontalalignment='center',
        verticalalignment='center',rotation=90,fontsize=12,color='red')

    
    leg.append(p1)
    leg.append(p2)

    plt.ylabel('Energy (kW)', fontsize=ft)
    plt.title(title, fontsize=ft)
    plt.xticks(ind+width*3.,NAME,fontsize=ft)
    plt.yticks(fontsize=ft)

    plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=ft)


    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()

def plot_A3(title):

    ft=16.

    casename='A_3'
    res_dir='../results/data/A/%s/summary.csv' #tool
    tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    lab=['Tonatiuh', 'SolTrace','Tracer', 'SOLSTICE','Heliosim','SolarPILOT']

    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):

            if casename in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)


    qmax=N.max(RESULTS[tools[i]][:,1].astype(float))
    qmin=N.min(RESULTS[tools[i]][:,2].astype(float))
    ymax=qmax*1.05
    ymin=qmin*0.8

    #============

    NAME=N.array(['A 3.1 \nNormal slope error \nPillbox sunshape ', 'A 3.2 \n Normal slope error \nBuie sunshape ' ])

    n=len(RESULTS[tools[0]])
    width=0.11
    ind=N.arange(0.08,(n)*0.8,0.8)

    leg=[]

    plt.figure(1,dpi=1000)
    fig, ax=plt.subplots()
    for i in xrange(len(tools)):
        
        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,color='0.9',label=tools[i], align='edge',edgecolor=['black']*n)
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,2].astype(float),color='0.2', align='edge',edgecolor=['black']*n)
        plt.ylim(ymin, ymax)  
        plt.xlim([0,n*0.8])

        for j in xrange(n):   
            ax.text(ind[j]+width*i+width/2.*1.22, ymin+(qmax-ymin)/(3.6*n)*(i+1), lab[i],horizontalalignment='center',
        verticalalignment='center',rotation=90,fontsize=12,color='red')

    
    leg.append(p1)
    leg.append(p2)

    plt.ylabel('Energy (kW)', fontsize=ft)
    plt.title(title, fontsize=ft)
    plt.xticks(ind+width*3.,NAME,fontsize=ft)
    plt.yticks(fontsize=ft)

    plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=ft)


    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()

def plot_B11(title):

    ft=16.

    casename='B_1.1'
    res_dir='../results/data/B/%s/summary.csv' #tool
    tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    lab=['Tonatiuh', 'SolTrace','Tracer', 'SOLSTICE','Heliosim','SolarPILOT']

    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):

            if casename in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)


    qmax=N.max(RESULTS[tools[i]][:,1].astype(float))
    qmin=N.min(RESULTS[tools[i]][:,2].astype(float))
    ymax=qmax*1.05
    ymin=qmin*0.8

    #============

    NAME=N.array([])

    for j in xrange(4):
        NAME=N.append(NAME, 'B 1.1.%s'%( j+1))


    n=len(RESULTS[tools[0]])
    width=0.12
    ind=N.arange(0.08,(n)*0.8,0.8)

    leg=[]

    plt.figure(1,dpi=1000)
    fig, ax=plt.subplots()
    for i in xrange(len(tools)):
        
        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,color='0.9',label=tools[i], align='edge',edgecolor=['black']*n)
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,2].astype(float),color='0.2', align='edge',edgecolor=['black']*n)
        plt.ylim(ymin, ymax)  
        plt.xlim([0,n*0.8])

        for j in xrange(n):   
            ax.text(ind[j]+width*i+width/2.*1.22, ymin+(qmax-ymin)/(3.2*n)*(i+1)+1.5, lab[i],horizontalalignment='center',
        verticalalignment='center',rotation=90,fontsize=10,color='red')

    
    leg.append(p1)
    leg.append(p2)

    plt.ylabel('Energy (kW)', fontsize=ft)
    plt.title(title, fontsize=ft)
    plt.xticks(ind+width*3.,NAME,fontsize=ft)
    plt.yticks(fontsize=ft)

    #plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=ft)
    plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],loc='best',fontsize=ft)


    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()

def plot_B12(title):

    ft=16.

    casename='B_1.2'
    res_dir='../results/data/B/%s/summary.csv' #tool
    tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    lab=['Tonatiuh', 'SolTrace','Tracer', 'SOLSTICE','Heliosim','SolarPILOT']

    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):

            if casename in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)


    qmax=N.max(RESULTS[tools[i]][:,1].astype(float))
    qmin=N.min(RESULTS[tools[i]][:,2].astype(float))
    ymax=qmax*1.05
    ymin=0.

    #============

    NAME=N.array([])

    for j in xrange(4):
        NAME=N.append(NAME, 'B 1.2.%s'%( j+1))


    n=len(RESULTS[tools[0]])
    width=0.12
    ind=N.arange(0.08,(n)*0.8,0.8)

    leg=[]

    plt.figure(1,dpi=1000)
    fig, ax=plt.subplots()
    for i in xrange(len(tools)):
        
        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,color='0.9',label=tools[i], align='edge',edgecolor=['black']*n)
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,2].astype(float),color='0.2', align='edge',edgecolor=['black']*n)
        plt.ylim(ymin, ymax)  
        plt.xlim([0,n*0.8])

        for j in xrange(n):   
            ax.text(ind[j]+width*i+width/2.*1.22, ymin+(qmax-ymin)/(7.*n)*(i+1)+5., lab[i],horizontalalignment='center',
        verticalalignment='center',rotation=90,fontsize=10,color='red')

    
    leg.append(p1)
    leg.append(p2)

    plt.ylabel('Energy (kW)', fontsize=ft)
    plt.title(title, fontsize=ft)
    plt.xticks(ind+width*3.,NAME,fontsize=ft)
    plt.yticks(fontsize=ft)

    #plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=ft)
    plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],loc='best',fontsize=ft)


    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()

def plot_B21(title):

    ft=16.

    casename='B_2.1'
    res_dir='../results/data/B/%s/summary.csv' #tool
    tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    lab=['Tonatiuh', 'SolTrace','Tracer', 'SOLSTICE','Heliosim','SolarPILOT']

    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):

            if casename in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)


    qmax=N.max(RESULTS[tools[i]][:,1].astype(float))
    qmin=N.min(RESULTS[tools[i]][:,2].astype(float))
    ymax=qmax*1.05
    ymin=qmin*0.8

    #============

    NAME=N.array([])

    for j in xrange(4):
        NAME=N.append(NAME, 'B 2.1.%s'%( j+1))


    n=len(RESULTS[tools[0]])
    width=0.12
    ind=N.arange(0.08,(n)*0.8,0.8)

    leg=[]

    plt.figure(1,dpi=1000)
    fig, ax=plt.subplots()
    for i in xrange(len(tools)):
        
        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,color='0.9',label=tools[i], align='edge',edgecolor=['black']*n)
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,2].astype(float),color='0.2', align='edge',edgecolor=['black']*n)
        plt.ylim(ymin, ymax)  
        plt.xlim([0,n*0.8])

        for j in xrange(n):   
            ax.text(ind[j]+width*i+width/2.*1.22, ymin+(qmax-ymin)/(3.2*n)*(i+1)+1.5, lab[i],horizontalalignment='center',
        verticalalignment='center',rotation=90,fontsize=10,color='red')

    
    leg.append(p1)
    leg.append(p2)

    plt.ylabel('Energy (kW)', fontsize=ft)
    plt.title(title, fontsize=ft)
    plt.xticks(ind+width*3.,NAME,fontsize=ft)
    plt.yticks(fontsize=ft)

    #plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=ft)
    plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],loc='best',fontsize=ft)

    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()

def plot_B22(title):

    ft=16.

    casename='B_2.2'
    res_dir='../results/data/B/%s/summary.csv' #tool
    tools=N.array(['tonatiuh', 'soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    lab=['Tonatiuh', 'SolTrace','Tracer', 'SOLSTICE','Heliosim','SolarPILOT']

    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):

            if casename in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)


    qmax=N.max(RESULTS[tools[i]][:,1].astype(float))
    qmin=N.min(RESULTS[tools[i]][:,2].astype(float))
    ymax=qmax*1.05
    ymin=0.

    #============

    NAME=N.array([])

    for j in xrange(4):
        NAME=N.append(NAME, 'B 2.2.%s'%( j+1))


    n=len(RESULTS[tools[0]])
    width=0.12
    ind=N.arange(0.08,(n)*0.8,0.8)

    leg=[]

    plt.figure(1,dpi=1000)
    fig, ax=plt.subplots()
    for i in xrange(len(tools)):
        
        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,color='0.9',label=tools[i], align='edge',edgecolor=['black']*n)
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,2].astype(float),color='0.2', align='edge',edgecolor=['black']*n)
        plt.ylim(ymin, ymax)  
        plt.xlim([0,n*0.8])

        for j in xrange(n):   
            ax.text(ind[j]+width*i+width/2.*1.22, ymin+(qmax-ymin)/(7.*n)*(i+1)+5., lab[i],horizontalalignment='center',
        verticalalignment='center',rotation=90,fontsize=10,color='red')

    
    leg.append(p1)
    leg.append(p2)

    plt.ylabel('Energy (kW)', fontsize=ft)
    plt.title(title, fontsize=ft)
    plt.xticks(ind+width*3.,NAME,fontsize=ft)
    plt.yticks(fontsize=ft)

    #plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=ft)
    plt.legend(leg,['$Q_{abs}$','$Q_{spil}$'],loc='best',fontsize=ft)

    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/%s.png'%casename,'w'), bbox_inches='tight')
    plt.close()

def plot_C(title):

    res_dir='../results/data/C/%s/summary.csv' #tool

    tools=N.array(['tonatiuh','soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    lab=['Tonatiuh', 'SolTrace','Tracer', 'SOLSTICE','Heliosim','SolarPILOT']

    RESULTS={}
    for i in xrange(len(tools)):
        RESULTS[tools[i]]=N.array([])
        res=N.loadtxt(res_dir%(tools[i]), dtype=str, delimiter=',')
        for j in xrange(len(res)):           
            if 'C_1.1' in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])
        for j in xrange(len(res)):      
            if 'C_2.1' in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])
        for j in xrange(len(res)):      
            if 'C_1.2' in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])
        for j in xrange(len(res)):      
            if 'C_2.2' in res[j,0]:
                RESULTS[tools[i]]=N.append(RESULTS[tools[i]], res[j])

        RESULTS[tools[i]]=RESULTS[tools[i]].reshape(len(RESULTS[tools[i]])/11, 11)
        #print ''
        #print RESULTS[tools[i]]
          
    n=len(RESULTS[tools[0]])
    ind=N.arange(0.2,(n)*0.8,0.8)
    width=0.11

    NAME=N.array(['C 1.1 \nSolar Noon \n(Pillbox)', 'C 2.1 \nSolar Noon \n(Buie)','C 1.2 \n Morning \n(Pillbox)', 'C 2.2 \n Morning \n(Buie)'   ])   
    leg=[]

    plt.figure(1,dpi=1000)
    fig, ax=plt.subplots()
    for i in xrange(len(tools)):
        #abs
        p1=plt.bar(ind+width*i,RESULTS[tools[i]][:,7].astype(float),width,color='0.9',label=tools[i], align='edge',edgecolor=['black']*n)
        # refl
        p2=plt.bar(ind+width*i,RESULTS[tools[i]][:,6].astype(float),width,bottom=RESULTS[tools[i]][:,7].astype(float),color='0.75', align='edge',edgecolor=['black']*n)
        # spil
        p3=plt.bar(ind+width*i,RESULTS[tools[i]][:,5].astype(float),width,bottom=RESULTS[tools[i]][:,7].astype(float)+RESULTS[tools[i]][:,6].astype(float),color='0.55', align='edge',edgecolor=['black']*n)
        # field absorb
        p4=plt.bar(ind+width*i,RESULTS[tools[i]][:,4].astype(float),width,bottom=RESULTS[tools[i]][:,7].astype(float)+RESULTS[tools[i]][:,6].astype(float)+RESULTS[tools[i]][:,5].astype(float),color='0.35', align='edge',edgecolor=['black']*n)        
        # block
        p5=plt.bar(ind+width*i,RESULTS[tools[i]][:,3].astype(float),width,bottom=RESULTS[tools[i]][:,7].astype(float)+RESULTS[tools[i]][:,6].astype(float)+RESULTS[tools[i]][:,5].astype(float)+RESULTS[tools[i]][:,4].astype(float),color='1.', align='edge',edgecolor=['black']*n)
        # shade+cos
        p6=plt.bar(ind+width*i,RESULTS[tools[i]][:,2].astype(float),width,bottom=RESULTS[tools[i]][:,7].astype(float)+RESULTS[tools[i]][:,6].astype(float)+RESULTS[tools[i]][:,5].astype(float)+RESULTS[tools[i]][:,4].astype(float)+RESULTS[tools[i]][:,3].astype(float),color='0.15', align='edge',edgecolor=['black']*n)

        for j in xrange(n):   
            ax.text(ind[j]+width*i+width/2.*1.22, 52000./(5.*n)*(i+1)+2000., lab[i],horizontalalignment='center',
        verticalalignment='center',rotation=90,fontsize=10,color='red')


        #plt.ylim(60.,101.)  
    leg.append(p1)
    leg.append(p2)
    leg.append(p3)
    leg.append(p4)
    leg.append(p5)
    leg.append(p6)

    plt.ylabel('Energy (kW)')
    plt.title(title)
    plt.xticks(ind+width*3., NAME)

    plt.legend(leg,['$Q_{abs}$','$Q_{refl}$','$Q_{spil}$', '$Q_{field,abs}$', '$Q_{block}$', '$Q_{shade&cos}$'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    save_en='./plots/energy-bar'
    if not os.path.exists(save_en):
        os.makedirs(save_en)
    plt.savefig(open(save_en+'/C.png','w'), bbox_inches='tight')
    plt.close()


if __name__=='__main__':
    plot_A11('Pillbox Slope Error')
    plot_A12('Normal Slope Error')
    plot_A2('Sunshape')    
    plot_A23('Buie Sunshape') 
    plot_A3('Combination of Sunshape and Slope Error') 
    plot_B11('Pillbox Sunshape \n Solar Noon') 
    plot_B12('Pillbox Sunshape \n Morning') 
    plot_B21('Buie Sunshape \n Solar Noon') 
    plot_B22('Buie Sunshape \n Morning') 
    plot_C('Round C') 
    #plot_energy('A_2', 'Sunshape')
    #plot_energy('A_2.3', 'Buie Sunshape')   
    #plot_energy_buie('A_2.3', 'Buie Sunshape')   
    #plot_energy('A_3', 'Combination of Sunshape and Slope Error')   
    #plot_energy('B_1', 'B1')
    #plot_energy('B_2','B2')        
    #plot_energy_C('Round C')
     
