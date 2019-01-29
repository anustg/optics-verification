import numpy as N
import os

def gen_case_name():
    name=N.array([])

    # A1
    for i in xrange(2):
        for j in xrange(3):
            case='A1.%s.%s'%(i+1,j+1)
            name=N.append(name, case)

    # A2.1 and A2.2
    for i in xrange(2):
        case='A2.%s'%(i+1)
        name=N.append(name, case)

    # A2.3
    for i in xrange(3):
        case='A2.3.%s'%(i+1)
        name=N.append(name, case)

    # A3.1 and A3.2
    for i in xrange(2):
        case='A3.%s'%(i+1)
        name=N.append(name, case)

    # Round B
    for h in xrange(2):
        for i in xrange(2):
            for j in xrange(4):
                case='B%s.%s.%s'%(h+1,i+1,j+1)
                name=N.append(name, case)

    # Round C
    for i in xrange(2):
        for j in xrange(2):
            case='C%s.%s'%(i+1,j+1)
            name=N.append(name, case)	

    return name

def get_results():
    CASE=gen_case_name()
    for i in xrange(len(CASE)):
        if i ==0:
            resfile='./tracer_results/%s/results.csv'%CASE[i]
            res=N.loadtxt(resfile, delimiter=',', dtype=str)
            title=N.array(['', CASE[i]]).reshape(1, 2)
            results=N.vstack((title, res))
        else:
            resfile='./tracer_results/%s/results.csv'%CASE[i]
            res=N.loadtxt(resfile, delimiter=',', dtype=str)
            title=N.array([CASE[i]])
            value=N.append(title, res[:,1])
            value=value.reshape(len(value), 1)
            results=N.hstack((results, value ))

    N.savetxt('./tracer_results/summary.csv', results, fmt='%s', delimiter=',')
            
def org_files(initial=False):

    if initial:
        os.makedirs('../results/data/A/tracer_test')
        os.makedirs('../results/data/B/tracer_test')    
        os.makedirs('../results/data/C/tracer_test')
    
    CASE=gen_case_name()

    title1=N.array(['Case name','Qin (kW)','Qabs (kW)','Qspil (kW)','Qrefl (kW)','Effective rays','Time (min)','num of process','CI_flux (kW/m2)','CI_en (kW)','From source'])
    title2=N.array(['Case name','Qall (kW)','Qshade+cos (kW)','Qblock (kW)','Qfield_abs (kW)','Qspil (kW)','Qrefl (kW)','Qabs (kW)','Effective rays','Time','From source'])
    A=N.array([])
    B=N.array([])
    C=N.array([])
    results_A=title1
    results_B=title1
    results_C=title2
    source='tracer_test'
    for i in xrange(len(CASE)):
        case=CASE[i]
        fluxfile='./tracer_results/%s/fluxmap.csv'%case
        newfile='../results/data/%s/tracer_test/%s_%s.csv'%(case[0], case[0], case[1:])
        os.system('cp %s %s'%(fluxfile, newfile))

        # RD A
        if case[0]=='A':
            casename='%s_%s'%(case[0], case[1:])
            result=N.loadtxt('./tracer_results/%s/results.csv'%case, delimiter=',', dtype=str)
            Qabs=result[7,1].astype(float)
            Qspil=result[5,1].astype(float)
            Qrefl=result[6,1].astype(float)
            Qin=Qabs+Qspil+Qrefl
            eff_rays=result[22,1].astype(float)#million
            time=result[27,1] #min
            num_process='-'
            CI_flux=result[16,1]
            CI_en=result[17,1]
            results_A=N.append(results_A,(casename,Qin,Qabs,Qspil,Qrefl,eff_rays,time,num_process,CI_flux,CI_en,source))
                       
        # RD B
        elif case[0]=='B':
            casename='%s_%s'%(case[0], case[1:])
            result=N.loadtxt('./tracer_results/%s/results.csv'%case, delimiter=',', dtype=str)
            Qabs=result[7,1].astype(float)
            Qspil=result[5,1].astype(float)
            Qrefl=result[6,1].astype(float)
            Qin=Qabs+Qspil+Qrefl
            eff_rays=result[22,1].astype(float)#million
            time=result[27,1] #min
            num_process='-'
            CI_flux=result[16,1]
            CI_en=result[17,1]
            results_B=N.append(results_B,(casename,Qin,Qabs,Qspil,Qrefl,eff_rays,time,num_process,CI_flux,CI_en,source))

        # RD C
        else:
            casename='%s_%s'%(case[0], case[1:])
            result=N.loadtxt('./tracer_results/%s/results.csv'%case, delimiter=',', dtype=str)
            Qall=result[0,1].astype(float) 
            Qshade=result[1,1].astype(float) 
            Qblock=result[2,1].astype(float) 
            Qfield=result[3,1].astype(float) 
            Qspil=result[5,1].astype(float) 
            Qrefl=result[6,1].astype(float)       
            Qabs=result[7,1].astype(float)
            eff_rays=result[22,1].astype(float)#million
            time=result[27,1] #min
            results_C=N.append(results_C,(casename,Qall, Qshade, Qblock, Qfield, Qspil, Qrefl, Qabs, eff_rays,time,source))
        
    results_A=results_A.reshape(len(results_A)/11,11)                
    results_B=results_B.reshape(len(results_B)/11,11) 
    results_C=results_C.reshape(len(results_C)/11,11) 
    
    N.savetxt('../results/data/A/tracer_test/summary.csv', results_A, fmt='%s', delimiter=',')
    N.savetxt('../results/data/B/tracer_test/summary.csv', results_B, fmt='%s', delimiter=',')
    N.savetxt('../results/data/C/tracer_test/summary.csv', results_C, fmt='%s', delimiter=',')

if __name__=='__main__':
    get_results()
    org_files()

    
