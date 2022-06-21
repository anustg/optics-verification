import numpy as N
from scipy.constants import degree


class Parameters:

    def __init__(self):
        pass


    def A1(self, err, dist):
        '''
        err - float, the deviation of slope error, rad
        dist - str, distribution of the slope error, 'normal' or 'pillbox' or 'perfect'
        '''
        # Sun
        #---------------
        self.DNI=1000.
        self.sunshape='collimated' #'pillbox' or 'gaussian' or 'buie' or 'collimated'
        self.sun_width=0.
        self.sun_az=180.*degree
        self.sun_zenith=0.*degree

        # Field
        #-----------------
        self.hst_w=10.
        self.hst_h=10.
        self.foc=500.
        self.reflectivity=1.
        self.slope_type=dist #'normal' or 'pillbox' or 'perfect'
        self.slope_error=err #rad (0 is a perfect mirror)
        self.curved=True # or False: flat mirror
        self.oneMirror=True # or True: for simulating just one mirror
        self.hst_file=None # director of the csv or None        
        self.index=-1 #
        self.pos=N.r_[0., 0., 0.] 
        self.tracking_mode='AzEl' # 'AzEl'or 'TiltRoll'
        self.aiming_mode='SinglePoint' # 'MultiFixed' or 'SinglePoint'

        # Receiver
        #-------------------
        self.rec_w=8.
        self.rec_h=8.
        self.absorptivity=1.
        self.rec_loc=N.r_[0.,0.,500.] # receiver location
        self.rec_rot=N.r_[180.,0.,0.]*degree # receiver rot through x, y, z (rad)
        
        # binning of fluxmap
        #------------------------
        self.m=100 #number of element in x direction
        self.n=100 # in y direction
        self.p=50  # in radius direction

    def A2(self, sunshape, CSR='-'):
        '''
        sunshape - str, 'pillbox' or 'gaussian' or 'buie' or 'collimated'
        CSR- if it's Buie sunshape
        '''
        # Sun
        #---------------
        self.DNI=1000.
        self.sunshape=sunshape #'pillbox' or 'gaussian' or 'buie' or 'collimated'
        if sunshape=='buie':
            self.sun_width=CSR
        else:
            self.sun_width=4.e-3
        self.sun_az=180.*degree
        self.sun_zenith=0.*degree

        # Field
        #-----------------
        self.hst_w=10.
        self.hst_h=10.
        self.foc=500.
        self.reflectivity=1.
        self.slope_type='perfect' #'normal' or 'pillbox' or 'perfect'
        self.slope_error=0 #rad (0 is a perfect mirror)
        self.curved=True # or False: flat mirror
        self.oneMirror=True # or True: for simulating just one mirror
        self.hst_file=None # director of the csv or None        
        self.index=-1 #
        self.pos=N.r_[0., 0., 0.] 
        self.tracking_mode='AzEl' # 'AzEl'or 'TiltRoll'
        self.aiming_mode='SinglePoint' # 'MultiFixed' or 'SinglePoint'

        # Receiver
        #-------------------
        self.rec_w=8.
        self.rec_h=8.
        self.absorptivity=1.
        self.rec_loc=N.r_[0.,0.,500.] # receiver location
        self.rec_rot=N.r_[180.,0.,0.]*degree # receiver rot through x, y, z (rad)
        
        # binning of fluxmap
        #------------------------
        self.m=100 #number of element in x direction
        self.n=100 # in y direction
        self.p=50  # in radius direction

    def A3(self, sunshape):
        '''
        sunshape - str, 'pillbox' or 'gaussian' or 'buie' or 'collimated'
        CSR- if it's Buie sunshape
        '''
        # Sun
        #---------------
        self.DNI=1000.
        self.sunshape=sunshape #'pillbox' or 'gaussian' or 'buie' or 'collimated'
        if sunshape=='pillbox':
            self.sun_width=4.65e-3
        else:
            self.sun_width=0.02
        self.sun_az=180.*degree
        self.sun_zenith=0.*degree

        # Field
        #-----------------
        self.hst_w=10.
        self.hst_h=10.
        self.foc=500.
        self.reflectivity=1.
        self.slope_type='normal' #'normal' or 'pillbox' or 'perfect'
        self.slope_error=2.e-3 #rad (0 is a perfect mirror)
        self.curved=True # or False: flat mirror
        self.oneMirror=True # or True: for simulating just one mirror
        self.hst_file=None # director of the csv or None        
        self.index=-1 #
        self.pos=N.r_[0., 0., 0.] 
        self.tracking_mode='AzEl' # 'AzEl'or 'TiltRoll'
        self.aiming_mode='SinglePoint' # 'MultiFixed' or 'SinglePoint'

        # Receiver
        #-------------------
        self.rec_w=8.
        self.rec_h=8.
        self.absorptivity=1.
        self.rec_loc=N.r_[0.,0.,500.] # receiver location
        self.rec_rot=N.r_[180.,0.,0.]*degree # receiver rot through x, y, z (rad)
        
        # binning of fluxmap
        #------------------------
        self.m=100 #number of element in x direction
        self.n=100 # in y direction
        self.p=50  # in radius direction

    def B(self, hst_position, sunshape, sun_az, sun_zenith):
        '''
        sunshape - str, 'pillbox' or 'gaussian' or 'buie' or 'collimated'
        CSR- if it's Buie sunshape
        '''
        # Sun
        #---------------
        self.DNI=1000.
        self.sunshape=sunshape #'pillbox' or 'gaussian' or 'buie' or 'collimated'
        if sunshape=='pillbox':
            self.sun_width=4.65e-3
        else:
            self.sun_width=0.02
        self.sun_az=sun_az*degree
        self.sun_zenith=sun_zenith*degree

        # Field
        #-----------------
        self.hst_w=10.
        self.hst_h=10.
        self.reflectivity=1.
        self.slope_type='normal' #'normal' or 'pillbox' or 'perfect'
        self.slope_error=2.e-3 #rad (0 is a perfect mirror)
        self.curved=True # or False: flat mirror
        self.oneMirror=True # or True: for simulating just one mirror
        self.hst_file=None # director of the csv or None        
        self.index=-1 #
        self.pos=hst_position 
        self.tracking_mode='AzEl' # 'AzEl'or 'TiltRoll'
        self.aiming_mode='SinglePoint' # 'MultiFixed' or 'SinglePoint'
        xp=self.pos[0]
        yp=self.pos[1]
        zp=self.pos[2]

        self.foc=N.sqrt(xp**2+yp**2+(zp-62.)**2)

        # Receiver
        #-------------------
        self.rec_w=8.
        self.rec_h=6.
        self.absorptivity=1.
        self.rec_loc=N.r_[0.,0.,62.] # receiver location
        self.rec_rot=N.r_[-90.,0.,0.]*degree # receiver rot through x, y, z (rad)
        
        # binning of fluxmap
        #------------------------
        self.m=100 #number of element in x direction
        self.n=100 # in y direction
        self.p=50  # in radius direction

    def C(self, sunshape, sun_az, sun_zenith):
        '''
        sunshape - str, 'pillbox' or 'gaussian' or 'buie' or 'collimated'
        CSR- if it's Buie sunshape
        '''
        # Sun
        #---------------
        self.DNI=1000.
        self.sunshape=sunshape #'pillbox' or 'gaussian' or 'buie' or 'collimated'
        if sunshape=='pillbox':
            self.sun_width=4.65e-3
        else:
            self.sun_width=0.02
        self.sun_az=sun_az*degree
        self.sun_zenith=sun_zenith*degree

        # Field
        #-----------------
        self.hst_w=10.
        self.hst_h=10.
        self.reflectivity=0.95
        self.slope_type='normal' #'normal' or 'pillbox' or 'perfect'
        self.slope_error=2.e-3 #rad (0 is a perfect mirror)
        self.curved=True # or False: flat mirror
        self.oneMirror=False # or True: for simulating just one mirror
        self.hst_file='../case_details/Round3_layout.csv' # director of the csv or None        
        self.index=-1 #
        self.pos='-'
        self.foc='-' 
        self.tracking_mode='AzEl' # 'AzEl'or 'TiltRoll'
        self.aiming_mode='SinglePoint' # 'MultiFixed' or 'SinglePoint'

        # Receiver
        #-------------------
        self.rec_w=8.
        self.rec_h=6.
        self.absorptivity=0.9
        self.rec_loc=N.r_[0.,0.,62.] # receiver location
        self.rec_rot=N.r_[-90.,0.,0.]*degree # receiver rot through x, y, z (rad)
        
        # binning of fluxmap
        #------------------------
        self.m=100 #number of element in x direction
        self.n=100 # in y direction
        self.p=50  # in radius direction


    def output_parameters(self, savefolder):
        '''
        Backup the simulated case parameters
        Arguement: 
            savefolder: str, the directory for saving the file
        '''

        constants=N.array((['Source (sun)',''],
            ['azimuth angle (deg):', self.sun_az/degree],
            ['zenith angle (deg):',self.sun_zenith/degree],
            ['sunshape:', self.sunshape],
            ['sun width parameter:',self.sun_width],
            ['DNI (W/m2):', self.DNI],
            [' ',' '],
            ['Heliostat',''],
            ['mirror width (m):', self.hst_w],
            ['mirror height (m):',self.hst_h],
            ['mirror reflectivity:',self.reflectivity],
            ['slope distribution:',self.slope_type],
            ['slope error:',self.slope_error],
            ['tracing system:',self.tracking_mode],
            ['aiming mode:',self.aiming_mode],
            ['',''],
            ['Receiver:',''],	
            ['W (m)', self.rec_w],
            ['H (m)',self.rec_h],
            ['rec absorptivity:',self.absorptivity],
            ['location:',str(self.rec_loc)],
            ['rotation:',str(self.rec_rot/degree)],		
            ['flux map x:',self.m],
            ['flum map y:',self.n],
            ['flux radius r:',self.p]))

        N.savetxt(savefolder+'/parameters.csv', constants, fmt='%s', delimiter=',')




