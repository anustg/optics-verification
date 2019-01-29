import numpy as N
import os
import matplotlib.pyplot as plt
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from matplotlib import gridspec


def plot_energy_C():
    res_dir='../results/data/C/%s/summary.csv' #tool

    tools=N.array(['tonatiuh','soltrace', 'tracer_test', 'solstice', 'heliosim','solarpilot'])
    tools_name=N.array(['Tonatiuh','SolTrace', 'Tracer', 'Solstice', 'Heliosim','SolarPILOT'])

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

    NAME=N.array(['C 1.1 (Pillbox)', 'C 2.1 (Buie)','C 1.2 (Pillbox)', 'C 2.2 (Buie)'])  

    plt.figure(1, dpi=1000)
    fig, axs = plt.subplots(2, 2, figsize=(12, 12), sharey='row')
    print N.shape(axs)    
    c=0
    fts=26
    mks=12


    for i in xrange(len(axs)):
        for j in xrange(len(axs[0])):
            summary=N.array([])
            x=N.arange(1,6)
            zero=N.zeros(5)
            ax=axs[i,j]
           
            # abs
            D=N.array([])
            tonatiuh=RESULTS['tonatiuh'][c][7].astype(float)
            for n in tools[1:]:
                d=(RESULTS[n][c][7].astype(float)-tonatiuh)/tonatiuh*100.
                D=N.append(D, d)
            ax.plot(x, D, 'o',markersize=mks, label='$Q_{abs}$')
            summary=N.append(summary, D)

            # refl
            D=N.array([])
            tonatiuh=RESULTS['tonatiuh'][c][6].astype(float)
            for n in tools[1:]:
                d=(RESULTS[n][c][6].astype(float)-tonatiuh)/tonatiuh*100.
                D=N.append(D, d)
            ax.plot(x, D, 'v',markersize=mks, label='$Q_{refl}$')
            summary=N.append(summary, D)

            # spil
            D=N.array([])
            tonatiuh=RESULTS['tonatiuh'][c][5].astype(float)
            for n in tools[1:]:
                d=(RESULTS[n][c][5].astype(float)-tonatiuh)/tonatiuh*100.
                D=N.append(D, d)
            ax.plot(x, D, '^',markersize=mks,label='$Q_{spil}$')
            summary=N.append(summary, D)

            # field abs
            D=N.array([])
            tonatiuh=RESULTS['tonatiuh'][c][4].astype(float)
            for n in tools[1:]:
                d=(RESULTS[n][c][4].astype(float)-tonatiuh)/tonatiuh*100.
                D=N.append(D, d)
            ax.plot(x, D, 'd',markersize=mks,label='$Q_{hstat,abs}$')
            summary=N.append(summary, D)

            # block
            D=N.array([])
            tonatiuh=RESULTS['tonatiuh'][c][3].astype(float)
            for n in tools[1:]:
                d=(RESULTS[n][c][3].astype(float)-tonatiuh)/tonatiuh*100.
                D=N.append(D, d)
            ax.plot(x, D, 's',markersize=mks,label='$Q_{block}$')
            summary=N.append(summary, D)

            # shade+cos
            D=N.array([])
            tonatiuh=RESULTS['tonatiuh'][c][2].astype(float)
            for n in tools[1:]:
                d=(RESULTS[n][c][2].astype(float)-tonatiuh)/tonatiuh*100.
                D=N.append(D, d)
            ax.plot(x, D, '*',markersize=mks,label='$Q_{shad}$+$Q_{cos}$')
            summary=N.append(summary, D)

            ax.set_xlim([0.8,5.2])
            ax.grid()
            ax.plot(x, zero, 'r--')
            ax.set_xticklabels(tools_name[:], rotation=45., fontsize=20)
            ax.set_title(NAME[c], fontsize=fts)

            ax.tick_params(axis='y', labelsize=fts)

            if c==1:
                ax.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0, fontsize=28)  
            c+=1 
            summary=summary.reshape(len(summary)/5, 5)
            summary=N.vstack((tools[1:], summary))

            led=N.array(['%','Qabs', 'Qrefl', 'Qspil', 'Qfieldabs', 'Qblok', 'Qshad&cos'])
            led=led.reshape(7, 1)
      
            summary=N.hstack((led, summary))
            N.savetxt('./plots/summary_C%s.%s_en_different.csv'%(i+1, j+1), summary, fmt='%s', delimiter=',') 
      
    fig.subplots_adjust(hspace=0.53)
    plt.subplots_adjust(right=0.85)
    plt.subplots_adjust(left=0.15)
    fig.text(0.5, 0.94, 'Solar Noon', ha='center',fontsize=26)
    fig.text(0.5, 0.45, 'Morning', ha='center',fontsize=26)
    fig.text(0.04, 0.5, 'Difference compared to Tonatiuh (%)', va='center', rotation='vertical', fontsize=fts)
    plt.savefig(open('./plots/C_diff.png','w'),bbox_inches='tight')

if __name__=='__main__':
    plot_energy_C()
        
        
            
