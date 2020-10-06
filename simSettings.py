import oceanState, oceanStateSplit

class SimSettings:

    def __init__(self):


        # Run duration:
        self.tEnd = 12*3600#3600*24*10

        # Save settings:
        self.saveIntS = 1200#1200 #180
        self.saveFile = 'C:/temp/real_gin_par.nc'
        self.saveAvgIntS = 3600*24*5
        self.saveAvgFile = 'C:/temp/aver1.nc'
        self.saveSubsetFile = 'C:/temp/real_gin_subset.nc'


        # Activate/deactivate features:
        self.modeSplittingOn = True
        self.useRk = True
        self.freshwaterOn = False
        self.coriolisOn = True
        self.atmoOn = True
        self.advectiveTermsOn = True # 1 to include advective terms in momentum equation, 0 to disable
        self.biharmonic = True # 1 to use biharmonic friction of horizontal velocities. 0 to use Smagorinsky
        self.trcVertMix = True
        self.trcHorizMix = True
        self.passiveTracer = True
        self.recordAverages = True
        self.implicitVerticalDiff = True # True to use implicit scheme for vertical diff og momentum in split

        # Mixing parameters:
        self.calcMixingInterval = 1
        self.KBi = 1e11 # Biharmonic coefficient. Should be overridden by scenario since it is resolution dependent
        self.A_z = 10*1e-2 # Vertical eddy viscosity used for vertical diffusion of momentum
        self.Ri0 = 0.65 # Threshold Richardson number value. Taken from SINMOD
        self.G_vmix = 30 # Shape parameter for Richardson number vertical mixing scheme
        self.KVm = 3e-2 # [m2/s] Maximum vertical diffusivity
        # Bottom friction:
        self.C_b = 2.5e-3 # taken from SINMOD's splitt_v7
        self.CM = 1.0 # For Smagorinsky
        self.CH = 0.4 # For Smagorinsky
        self.CM2D = 1.0  # For Smagorinsky
        # Coriolis parameters:
        self.omega = 0.7292e-4 # Earth rotational speed(rad / s)
        self.latitude = 73
        self.phi0 = self.latitude*3.1415/180.0
        if not self.coriolisOn:
            self.omega = 0

        # Atmosphere / waves:
        self.atmoUpdateInterval = 300 # Update interval for atmo data(s)
        self.p_atm0 = 101000 # Atmospheric pressure(Pa)
        self.H_wave = 0.5 # Wave height(m)
        self.T_wave = 3 # Wave period(s)
        self.windStressCoeff = 0.00015 #0.0002 # Multiply wind speed (m/s) by this coeff to get wind stress

        self.rho_0 = 1023.6

    # Instantiate ocean object. The type of ocean object is decided based on whether
    # mode splitting is activated or not.
    def getOcean(self, imax, jmax, kmax):
        if self.modeSplittingOn:
            return oceanStateSplit.OceanStateSplit(imax, jmax, kmax)
        else:
            return oceanState.OceanState(imax, jmax, kmax)