import numpy as N


def get_energy(e, pst_rec, hst_w, hst_h, hst_num, DNI, reflectance):
    '''
    Use bundle tree structure to get 
    incident energy, spillage, reflection and absorption

    Arguements:

    e - TracerEngine
    pst_rec - float, receiver position
    hst_w -float, the width of heliostat
    hst_h -float, the height of heliostat
    hst_num - int, number of heliostat
    DNI - float, direct normal intensity (W/m2)
    reflectance - float, the reflectivity of heliostat

    Returns:
    Q_total_inc - float, the total theoretical incident power on the field (kW)
    Q_to_rec00 - float, the power from the source that arrive directly onto the receiver (kW) 
    Q_shade - float, the losses due to heliostats shading (kW)
    Q_hst_abs - float, the losses due to non-deal reflection of heliostats (kW)
    Q_block - float, the losses due to heliostats blockage (kW)
    Q_in - float, the incident power to the receiver (kW), which is equal to Q_spil+Q_refl+Q_abs
    Q_spil - float, the power irradiate from the heliostats towards to the receiver but missed the receiver, (kW)
    Q_refl - float, the power reflected by the receiver (kW)
    Q_abs - float, the power from the heliostat and absorbed by the receiver (kW)
    eff_rays - int, the number of effective rays that hit an object from the primary shoot
    success - bool, it is True if rays reached the receiver from the reflection of helisotats 
    '''

    energy_bund={}
    hits=len(e.tree._bunds)
    for h in xrange(hits):
        energy_bund[h]=e.tree._bunds[h].get_energy()/1000. #kW

    if hits>2:

        # energy incident to the heliostat field
        A_heliostat=hst_w*hst_h*float(hst_num)
        Q_total_inc=A_heliostat*DNI/1000. #kW

        z1=e.tree._bunds[1].get_vertices()[2]
        z2=e.tree._bunds[2].get_vertices()[2]
    
        # energy incident on the heliostat field (after shading and cosine effect)
        Q_hst_t=0. 
        idx_0=e.tree._bunds[1].get_parents()
        N_hst_t=0
        for l in xrange(len(z1)):
            if z1[l]<pst_rec:
                # total energy emitted from heliostats in bundle[1]
                Q_hst_t+=energy_bund[1][l]
                N_hst_t+=1

        #energy incident on the heliostat field (after shading)--
        Q_inc_field=Q_hst_t/(reflectance)

        # method #2 for Q_inc_field
        Q_inc_field2=0.
        for i in xrange(len(idx_0)):
            if z1[i]<pst_rec:
                Q_inc_field2+=energy_bund[0][idx_0[i]]
        #print '--------------------------'
        #print 'Q_inc_field 1', Q_inc_field
        #print 'Q_inc_field 2', Q_inc_field2
        #print '----------------------------'      


        # absorption by heliostats
        Q_hst_abs=Q_inc_field-Q_hst_t 

        # blocking method #1
        Q_hst_hst=0.
        index_no_spil=e.tree._bunds[2].get_parents()
        for l in xrange(len(z2)):
            if z2[l]<pst_rec:
                index_to_hst=index_no_spil[l]
                if z1[index_to_hst]<pst_rec:
                    Q_hst_hst+=energy_bund[1][index_to_hst]
        Q_block=Q_hst_hst
                
        #----spillage------
        Q_hst_nl=0.	
        for index in index_no_spil:	
            if z1[index]<pst_rec:
                # energy emitted from heliostats in bundle[1] that didn't miss an object
                Q_hst_nl+=energy_bund[1][index]
        Q_spil=Q_hst_t-Q_hst_nl

        #----incident------
        # Q_to_rec_00
        # from bundle[0](from the sun) directly to the receiver
        over_rec1=z1>pst_rec
        N_to_rec00=N.sum(over_rec1)
        N_total00=e.tree._bunds[0].get_num_rays()
        Q_to_rec00=float(N_to_rec00)/float(N_total00)*N.sum(energy_bund[0])

        # shading and cosine
        Q_shade=Q_total_inc-Q_inc_field-Q_to_rec00 

        # Q_to_rec_11
        # from bundle[1](mainly from heliostats) to the receiver
        Q_rec_rec=0.
        Q_hst_rec=0.
        for l in xrange(len(z2)):
            if z2[l]>pst_rec:
                index_to_rec=index_no_spil[l]
                if z1[index_to_rec]>pst_rec:
                    # from receiver to receiver in bundle [1]
                    Q_rec_rec+=energy_bund[1][index_to_rec]
                else:
                    # from heliostat to receiver in bundle [1]
                    Q_hst_rec+=energy_bund[1][index_to_rec]
        Q_to_rec11=Q_hst_rec

        Q_in=Q_to_rec00+Q_to_rec11+Q_spil

        #print 'spillage'
        #print ' 1', Q_spil
        #print ' 2', Q_hst_t-Q_block-Q_hst_rec
        #print '--------------'

       # blocking method #2
        Q_block2=Q_hst_t-Q_in+Q_to_rec00                   
        #print ''
        #print 'block'
        #print ' 1', Q_block
        #print ' 2', Q_block2
        #print '--------------'
        #print Q_to_rec00
        #print '----------------'

     
        #----reflection------
        Q_rec_t=0.
        Q_rec_rec_h=0.
        Q_received=0.
        for h in xrange(1,hits):
            z_h=e.tree._bunds[h].get_vertices()[2]
            index_h=e.tree._bunds[h].get_parents()
            for l in xrange(len(z_h)):
                if z_h[l]>pst_rec:
                    #energy emitted from receiver in bund[h]
                    Q_rec_t+=energy_bund[h][l]
                    Q_received+=energy_bund[h-1][index_h[l]]
            if h<hits-1:
                # z vertices of rays in bund[h+1],
                # includes those from both rec and helios
                z_h1=e.tree._bunds[h+1].get_vertices()[2]
                # index of rays in bund[h] irradiate to both rec and helios
                index_nl_h=e.tree._bunds[h+1].get_parents()

                for l in xrange(len(z_h1)):
                    if z_h1[l]>pst_rec:
                        index_to_rec_h=index_nl_h[l]
                        if z_h[index_to_rec_h]>pst_rec:	
                            Q_rec_rec_h+=energy_bund[h][index_to_rec_h]
                        else:
                            if h>1:
                                # from heliostat to receiver in bundle [h]
                                Q_rec_rec_h+=energy_bund[h][index_to_rec_h]  

        Q_refl=Q_rec_t-Q_rec_rec_h

        Q_abs=Q_received-Q_rec_t

        # the effective rays (in bundle [1])
        eff_rays=len(e.tree._bunds[1].get_vertices()[2])      
        
        #print '************************************************'
        #print 'Q_sun', Q_total_inc
        #print 'Q_shade',Q_shade
        #print 'Q_abs_field',Q_hst_abs
        #print 'Q_Block', Q_block
        #print 'Qspil', Q_spil
        #print 'Qin', Q_in
        #print 'Qrefl', Q_refl
        #print 'Qabs',Q_abs
        #print 'close1?', Q_total_inc+Q_to_rec00-Q_shade-Q_block-Q_hst_abs-Q_in
        #print 'close2?',Q_in-Q_spil-Q_refl-Q_abs
        #print '' 
        #print '************************************************' 
        success=True    
    else:
        Q_total_inc, Q_to_rec00, Q_shade, Q_hst_abs, Q_block, Q_in, Q_spil, Q_refl, Q_abs, eff_rays =N.zeros(10)
        success=False              


    return Q_total_inc, Q_to_rec00, Q_shade, Q_hst_abs, Q_block, Q_in, Q_spil, Q_refl, Q_abs, eff_rays, success


 


