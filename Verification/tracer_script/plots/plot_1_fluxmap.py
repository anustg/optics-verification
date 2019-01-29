import numpy as N
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def plot_fluxmap(casename, savefolder, W, H):
    '''
    casename - str, eg.'A1.1.1'
    savefolder - str, the directory of the folder for saving and storing the results
    W - float, target width (m)
    H - float, target height (m)
    '''
    fluxfile=N.loadtxt('%s/%s/fluxmap.csv'%(savefolder, casename), delimiter=',', skiprows=1)
    X=fluxfile[:, 0].astype(float).reshape(100, 100)
    Y=fluxfile[:, 1].astype(float).reshape(100, 100)
    flux=fluxfile[:, -1].astype(float).reshape(100, 100)

    plt.figure(1)
    plt.pcolormesh(X,Y,flux,cmap=cm.hot)
    plt.xticks((-W/2., 0, W/2.), size=14)   
    plt.yticks((-H/2., 0, H/2.), size=14)
    plt.axis('equal',fontsize=14)
    plt.colorbar()
    plt.savefig(open(savefolder+'/'+casename+'/fluxmap.png','w'), dpi=1000, bbox_inches='tight')
    plt.close()


if __name__=="__main__":
    plot_fluxmap('B1.2.1', './tracer_results', 8., 6.)

    
