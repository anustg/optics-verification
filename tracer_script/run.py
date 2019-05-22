import math
import numpy as N
import time
import os


from tracer.models.heliostat_field import HeliostatGenerator, KnownField, FlatOneSidedReceiver, MountReceiver, solar_vector, TowerScene
from tracer.sources import rect_ray_bundle, buie_integration
from tracer.tracer_engine import TracerEngine
from tracer.CoIn_rendering.rendering import *

from parameters import Parameters
from additional import confidence_interval, bin_element, bin_radius
from EnergyBalance import get_energy


class RayTracing(Parameters):
    '''
    use the TowerScene to do ray-tracing simulations
    '''

    def __init__(self, num_rays, rendering, case,hemisphere='North'):
        '''
        num_rays - int, number of rays
        rendering - boolean, visualise the scene
        '''

        Parameters.__init__(self)
        self.num_rays=num_rays
        self.rendering=rendering
        self.case=case       
        self.hemisphere=hemisphere

        # set the convergent standard
        self.convergent='flux' # 'energy' or 'flux' or 'iters' or 'rays'
        self.cvg_value=5.e-2
        self.savefolder='./tracer_results/%s'%case
        if not os.path.exists(self.savefolder):
            os.makedirs(self.savefolder)

        self.gen_plant()

    def gen_plant(self):
        # define the case
        if self.case[0]=='A':
            if self.case[1]=='1':
                if self.case[:4]=='A1.1':
                    distribution='pillbox'
                    err=float(self.case[-1])*1.e-3
                    self.A1(err, distribution)
                elif self.case[:4]=='A1.2':
                    distribution='normal'
                    err=float(self.case[-1])*1.e-3
                    self.A1(err, distribution)
            elif self.case[1]=='2':
                if self.case[3]=='1':
                    sunshape='pillbox'
                    self.A2(sunshape)
                elif self.case[3]=='2':
                    sunshape='gaussian'
                    self.A2(sunshape)
                else:
                    sunshape='buie'
                    CSR=float(self.case[-1])*1.e-2
                    self.A2(sunshape, CSR)
            else:
                if self.case[-1]=='1':
                    sunshape='pillbox'
                else:
                    sunshape='buie'
                self.A3(sunshape)

        elif self.case[0]=='B':
            if self.case[1]=='1':
                sunshape='pillbox'
            else:
                sunshape='buie'

            if self.case[3]=='1':
                # solar noon
                sun_az=0.
                sun_zenith=12.
            else:
                # morning
                sun_az=-104.
                sun_zenith=68.

            if self.case[-1]=='1':
                hst_pos=N.r_[0., 46.5, 0.]
            elif self.case[-1]=='2':
                hst_pos=N.r_[0., 536.9, 0.]
            elif self.case[-1]=='3':
                hst_pos=N.r_[-324.3, 427.9, 0.]
            else:
                hst_pos=N.r_[252.5, 118.1, 0.]

            self.B(hst_pos, sunshape, sun_az, sun_zenith)

        elif self.case[0]=='C':
            if self.case[1]=='1':
                sunshape='pillbox'
            else:
                sunshape='buie'

            if self.case[3]=='1':
                # solar noon
                sun_az=0.
                sun_zenith=12.
            else:
                # morning
                sun_az=-104.
                sun_zenith=68.

            self.C(sunshape, sun_az, sun_zenith)

        self.output_parameters(self.savefolder)

        #--------------
        #    field
        #--------------
        heliostat=HeliostatGenerator(self.hst_w, self.hst_h, absorptivity=1.-self.reflectivity, sigma_xy=self.slope_error, slope=self.slope_type,curved=self.curved, one_mirror=self.oneMirror, index=self.index, pos=self.pos, foc=self.foc)
        # layout and field
        layout=KnownField(self.hst_file,self.pos,self.foc)

        #---------------
        #   receiver
        #---------------
        self.loc_z=self.rec_loc[-1]
        receiver=FlatOneSidedReceiver(self.rec_w,self.rec_h,self.absorptivity)
        mount_rec=MountReceiver(self.rec_loc,self.rec_rot)

        #--------------
        #     solar
        #--------------
        self.sun_vec=solar_vector(self.sun_az, self.sun_zenith)
        if self.sunshape=='buie':
            self.BUIE=buie_integration(CSR=self.sun_width,  preproc_CSR='CA')
        else:
            self.BUIE=None

        # Creating the tower scene
        #------------------------------
        tower_scene=TowerScene(self.sun_vec, self.oneMirror)
        tower_scene(heliostat, layout, self.aiming_mode, self.tracking_mode,receiver,mount_rec)

        self.system=tower_scene.system
        self.pos=heliostat.pos

        if self.oneMirror:
            self.num_hst=1
        else:
            self.num_hst=len(self.pos)


    def gen_rays(self):
        #===================================================
        if self.oneMirror:
            field_centre=self.pos
            x=self.hst_w*1.2
            y=self.hst_h*1.2
            centre = N.c_[10.*self.sun_vec + field_centre]

        else:
            t_pos = self.pos.T        
            xc_min = t_pos[0][N.argmin(t_pos[0])]
            xc_max = t_pos[0][N.argmax(t_pos[0])]
            yc_min = t_pos[1][N.argmin(t_pos[1])]
            yc_max = t_pos[1][N.argmax(t_pos[1])]

            if yc_min==yc_max:
	            y_dist=self.hst_height
            else:
	            y_dist = yc_max - yc_min
            if xc_min==xc_max:
	            x_dist=self.hst_width
            else:
	            x_dist = xc_max - xc_min

            xc_cent = (xc_min + xc_max) / 2.
            yc_cent = (yc_min + yc_max) / 2.
            field_centre = N.r_[xc_cent, yc_cent, 0.]

            # disc source covering entire field area:
            radius = 1.10 * math.sqrt((x_dist/2.)**2 + (y_dist/2.)**2)
            x=x_dist*1.2
            y=y_dist*1.2
            centre = N.c_[500.*self.sun_vec + field_centre]

        direction = N.array(-self.sun_vec)	
        rays=rect_ray_bundle(self.num_rays, centre, direction, x, y, self.sunshape, ang_rang=self.sun_width, flux=self.DNI, BUIE=self.BUIE)
        return rays


    def trace(self):
        '''
        Raytrace method.

        Raytraces successive bundles and stores the results of the fluxmap on the receiver.
        '''
        #==============================
        # Initialise the simulation
        #==============================

        # record processing time
        timer_rays=0.
        timer_mcrt=0.
        timer_postprocess=0.
        timer_total=time.time()

        # define bins for fluxmap
        bins, bin_area, fluxmap = bin_element(self.m, self.n, self.rec_w, self.rec_h)
        rbins, rb_area, r_flux = bin_radius(self.rec_w, self.rec_h, self.p)
        dx=self.rec_w/float(self.m)
        dy=self.rec_h/float(self.n)
        xf=N.linspace(-self.rec_w/2.+dx/2.,self.rec_w/2.-dx/2.,self.m)
        yf=N.linspace(-self.rec_h/2.+dy/2.,self.rec_h/2.-dy/2.,self.n)

        # initialise energy and Confident Interval (CI)
        Q_rec_00=0.
        Q_shade=0.
        Q_hst_abs=0.
        Q_block=0
        Q_in=0.
        Q_refl=0.
        Q_spil=0.
        Q_abs=0.

        abs_front=0.
        abs_back=0.

        ci_flux=10000
        ci_in=10000
        ci_shade=10000
        ci_field=10000
        ci_block=10000

        X_flux=0
        X_in=0
        X_shade=0
        X_block=0
        X_field=0

        Q_pk=N.array(())
        Q_in_all=N.array(())
        Q_shade_all=N.array(())
        Q_field_all=N.array(())
        Q_block_all=N.array(())

        CI_flux=N.array(([]))
        CI_in=N.array(([]))
        CI_shade=N.array(([]))
        CI_field=N.array(([]))
        CI_block=N.array(([]))

        i=0
        eff_rays=0.
        peak=300.
        Q_in=5000.

        if self.convergent=='energy':
            a=ci_shade
            b=self.cvg_value*Q_shade                             
        elif self.convergent=='flux':
            a=ci_flux
            b=self.cvg_value*peak
        elif self.convergent=='iters':
            a=self.cvg_value
            b=i
        elif self.convergent=='rays':
            a=self.cvg_value
            b=eff_rays

        #====================
        #   Ray Tracing
        #====================

        # define the tracer engine
        sys=time.time()
        e = TracerEngine(self.system)
        e.minerer=1e-10
        timer_sys=time.time()-sys

        # iterations on ray-tracing

        while (a>b or i<=10):
            print ''
            print ''
            print 'ITERATION ', i+1 

            # generate rays
            ray_time=time.time()
            rays=self.gen_rays()
            timer_rays+=time.time()-ray_time

            # MC ray-tracing
            mcrt=time.time()
            e.ray_tracer(rays, reps=100, min_energy=1e-10)
            timer_mcrt+=time.time()-mcrt

            if self.rendering:
                trace_scene=Renderer(e)
                trace_scene.show_rays(resolution=10)

            # Get results
            # flux for each surface from Absorption Account
            postprocess=time.time()
            front =self.system._objects[0].get_surfaces()[0]
            back = self.system._objects[0].get_surfaces()[1]

            # the absorbed energy from absorption account
            en, pts = front.get_optics_manager().get_all_hits()
            en2, pts2 = back.get_optics_manager().get_all_hits()	
            abs_front=(abs_front*float(i)+N.sum(en)/1000.)/float(i+1)
            abs_back =(abs_back*float(i)+N.sum(en2)/1000.)/float(i+1)

            # flux distribution on the front surface
            x, y = front.global_to_local(pts)[:2]   
            H, xbins, ybins = N.histogram2d(x, y, bins, weights=en)             
            fluxmap_i=H/(1000.*bin_area)

            fluxmap = (fluxmap*float(i)+fluxmap_i)/float(i+1)
            peak=N.max(fluxmap)
            average=N.sum(fluxmap)/float(self.m)/float(self.n)

            r=N.sqrt(x**2+y**2)
            rh,redge=N.histogram(r,rbins,weights=en)
            r_flux=(r_flux*float(i)+rh/(1000.*rb_area))/float(i+1)
            radius=redge[1:]-(redge[1]-redge[0])/2.  

            # get the field-wide energy
            Q_total_inc,Q_rec00_i, Q_shade_i,Q_hst_abs_i, Q_block_i, Q_in_i,Q_spil_i,Q_refl_i,Q_abs_i,eff_rays_i,success=get_energy(e, self.loc_z-self.rec_h*3., self.hst_w, self.hst_h, self.num_hst, self.DNI, self.reflectivity)
            if success:
                Q_rec_00=(Q_rec_00*float(i)+Q_rec00_i)/float(i+1)
                Q_shade=(Q_shade*float(i)+Q_shade_i)/float(i+1) 
                Q_hst_abs=(Q_hst_abs*float(i)+Q_hst_abs_i)/float(i+1)  
                Q_block=(Q_block*float(i)+Q_block_i)/float(i+1)            
                Q_in=(Q_in*float(i)+Q_in_i)/float(i+1)
                Q_spil=(Q_spil*float(i)+Q_spil_i)/float(i+1)
                Q_refl=(Q_refl*float(i)+Q_refl_i)/float(i+1)
                Q_abs=(Q_abs*float(i)+Q_abs_i)/float(i+1)     
                eff_rays+=eff_rays_i  

                Q_abs_acc=abs_front+abs_back

                ci_flux, X_flux=confidence_interval(fluxmap_i,i, X_flux,fluxmap)
                Ci_flux=N.max(ci_flux)

                ci_in, X_in=confidence_interval(Q_in_i, i, X_in, Q_in)
                ci_shade, X_shade=confidence_interval(Q_shade_i, i, X_shade, Q_shade)
                ci_block, X_block=confidence_interval(Q_block_i, i, X_block, Q_block)
                ci_field, X_field=confidence_interval(Q_hst_abs_i, i, X_field,Q_hst_abs)
         
                if (i%10==0 and i!=0):             
                   
                    CI_in=N.append(CI_in, ci_in)
                    CI_shade=N.append(CI_shade, ci_shade)
                    CI_field=N.append(CI_field, ci_field)
                    CI_block=N.append(CI_block, ci_block)

                    Q_in_all=N.append(Q_in_all, Q_in)
                    Q_shade_all=N.append(Q_shade_all, Q_shade)
                    Q_field_all=N.append(Q_field_all, Q_hst_abs)
                    Q_block_all=N.append(Q_block_all, Q_block)

                    Q_pk=N.append(Q_pk, peak)
                    CI_flux=N.append(CI_flux, Ci_flux)

                # save results
                if (i%200==0 and i!=0):
    
                    results=N.array((
                     ['Q_total_field_area (kW)',Q_total_inc],
                     ['Q_shade_and_cosine (kW)',Q_shade],
                     ['Q_block (kW)',Q_block],
                     ['Q_field_abs (kW)',Q_hst_abs],
                     ['', ''],
		             ['Q_spil (kW)',Q_spil],
		             ['Q_refl (kW)',Q_refl],
		             ['Q_abs (kW)',Q_abs],
                     ['', ''],
		             ['receiver efficiency ',Q_abs/Q_in],
		             ['F_spil 100%',Q_spil/Q_in],
		             ['F_refl 100%',Q_refl/Q_in],
                     ['', ''],
		             ['Peak Flux (kW/m2)',peak],
                     ['Average Flux (kW/m2)',average],
                     ['', ''],
			         ['CI_flux',Ci_flux],
                     ['CI_en',ci_in],
                     ['CI_shade',ci_shade],
                     ['CI_field',ci_field],
                     ['CI_block',ci_block],
                     ['', ''],
			         ['effecitve rays (million)',eff_rays/1.e6],
			         ['Time gen system (min)',timer_sys/60.],
			         ['Time gen rays (min)',timer_rays/60.],
			         ['Time MCRT (min)', timer_mcrt/60.],
		             ['Time postprocessing (min)',timer_postprocess/60.],
			         ['Total time (min)',total_time/60.]))
                    N.savetxt(self.savefolder+'/results.csv', results, fmt='%s', delimiter=',')
                    
                    title=N.array(['x(m)', 'y(m)', 'flux (kW/m2)'])

                    flux=fluxmap.T
                    # view the flux distribution at the position faces to the tower
                    if self.hemisphere=='North':
                        flux=N.fliplr(flux)
                        flux=N.flipud(flux)
                    flux=flux.flatten()

                    X,Y=N.meshgrid(xf, yf)
                    X=X.flatten()
                    Y=Y.flatten()
                    
                    fluxfile=N.hstack((X,Y))
                    fluxfile=N.hstack((fluxfile, flux))
                    fluxfile=fluxfile.reshape(3,len(fluxfile)/3)
                    fluxfile=fluxfile.T
                    title=N.array(['x (m)','y (m)', 'flux (kW/m2)'])
                    fluxfile=N.vstack((title, fluxfile))             

                    N.savetxt(self.savefolder+'/fluxmap.csv', fluxfile, fmt='%s', delimiter=',')
                                        
                timer_postprocess+=time.time()-postprocess
                total_time=time.time()-timer_total
                print ' Case', self.case
                print ' time_sys ', timer_sys/60., 'min'
                print ' time rays', timer_rays/60., 'min'
                print ' time_mcrt', timer_mcrt/60., 'min'
                print ' time processing ', timer_postprocess/60., 'min'   
                print ' total time', total_time/60., 'min'               
                print ' '
                print ' Q_sun:', Q_total_inc
                print ' Q_shade:', Q_shade, '+/-', ci_shade
                print ' Q_abs_field:',Q_hst_abs, '+/-', ci_field
                print ' Q_block:', Q_block, '+/-', ci_block
                print ' Q_in:',Q_in, '+/-', ci_in
                print ' Q_abs:', Q_abs
                print ' Q_spil:', Q_spil
                print ' Q_ref:',Q_refl
                print ' close1??:', Q_total_inc-Q_shade-Q_block-Q_hst_abs-Q_in
                print ' close2??:', Q_in-Q_abs-Q_spil-Q_refl
                print ' '
                print ' Peak flux (kW/m2):',peak, '+/-', Ci_flux
                print ' AVG flux(kW/m2): ', N.max(average)  
                print ' '
                print ' efficiency=Q_abs/Q_in:', Q_abs/Q_in
                print ' spil%= Q_spil/Q_in:', Q_spil/Q_in
                print ' ref%=  Q_ref/Q_in:',  Q_refl/Q_in
                print ''
                print ' Q_acc',Q_abs_acc
                print ' Q_tree',Q_abs
                print ''           
                print ' effective rays',eff_rays

            #==============================================================
            e.tree._bunds = []
            for clear in xrange(len(e._asm.get_surfaces())):
                e._asm.get_surfaces()[clear].get_optics_manager().reset()
            #============================================================== 
           
            i+=1
            if self.convergent=='energy':
                a=ci_shade
                b=self.cvg_value*Q_shade                             
            elif self.convergent=='flux':
                a=Ci_flux
                b=self.cvg_value*peak
            elif self.convergent=='iters':
                b=i
            elif self.convergent=='rays':
                b=eff_rays

		
if __name__=='__main__':
	
	rt=RayTracing(num_rays=1000, rendering=False, case='C1.2')
	rt.trace()	
