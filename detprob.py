import sys,os,contextlib,io
import numpy as np
import cPickle as pickle
import astropy.cosmology
import scipy.stats
import scipy.interpolate

@contextlib.contextmanager
def nostdout():
    ''' Locally suppress stoud print. See https://stackoverflow.com/a/2829036'''
    save_stdout = sys.stdout
    sys.stdout = io.BytesIO()
    yield
    sys.stdout = save_stdout

with nostdout(): # pycbc prints a very annyoing new line when imported
    import pycbc.waveform
    import pycbc.filter
    import pycbc.psd


def plotting():
    from matplotlib import use #Useful when working on SSH
    use('Agg')
    from matplotlib import rc #Set plot defaults
    font = {'family':'serif','serif':['cmr10'],'weight' : 'medium','size' : 16}
    rc('font', **font)
    rc('text',usetex=True)
    rc('figure',max_open_warning=1000)
    #matplotlib.rcParams['xtick.top'] = True
    #matplotlib.rcParams['ytick.right'] = True
    global plt
    import matplotlib.pyplot as plt


def singleton(class_):
    ''' Implement singleton design pattern in python.
    See https://stackoverflow.com/a/6810621 '''

    class class_w(class_):
        _instance = None
        def __new__(class2, *args, **kwargs):
            if class_w._instance is None:
                class_w._instance = super(class_w, class2).__new__(class2, *args, **kwargs)
                class_w._instance._sealed = False
            return class_w._instance
        def __init__(self, *args, **kwargs):
            if self._sealed:
                return
            super(class_w, self).__init__(*args, **kwargs)
            self._sealed = True
    class_w.__name__ = class_.__name__
    return class_w





@singleton
class Pomega(object):

    def __init__(self,  mcn=int(1e8),
                        mcbins=int(1e5),
                        binfile='Pwint.pkl'):

        self._interpolate = None
        self.mcn = mcn
        assert isinstance(self.mcn,(int,long))
        self.mcbins = mcbins
        assert isinstance(self.mcbins,(int,long))
        self.binfile=binfile

    def montecarlo_samples(self,mcn):

        # Sample the angles according to their pdf
        theta=np.arccos(np.random.uniform(-1,1,mcn))   # Polar
        phi = np.pi*np.random.uniform(-1,1,mcn)        # Azimuthal
        psi = np.pi*np.random.uniform(-1,1,mcn)        # Azimuthal
        iota = np.arccos(np.random.uniform(-1,1,mcn))  # Polar

        # Antenna patterns. Eq. (57) in arxiv:0903.0338
        Fplus  = 0.5*(1.+np.cos(theta)**2)*np.cos(2.*phi)*np.cos(2.*psi) - np.cos(theta)*np.sin(2.*phi)*np.sin(2.*psi)
        Fcross = 0.5*(1.+np.cos(theta)**2)*np.cos(2.*phi)*np.sin(2.*psi) + np.cos(theta)*np.sin(2.*phi)*np.cos(2.*psi)

        # Projection parameter. Eq. (3.31) in arxiv:gr-qc/9301003 but define omega=Theta/4
        littleomega = ( Fplus**2*(1.+np.cos(iota)**2)**2/4. +Fcross**2*np.cos(iota)**2 )**0.5

        return littleomega if len(littleomega)>1 else littleomega[0]


    def interpolate(self):

        if self._interpolate is None:

            # Takes some time. Store a pickle...
            if not os.path.isfile(self.binfile):

                print('['+self.__class__.__name__+'] Interpolating...')
                hist = np.histogram(self.montecarlo_samples(self.mcn),bins=self.mcbins)
                hist_dist = scipy.stats.rv_histogram(hist)
                with open(self.binfile, 'wb') as f: pickle.dump(hist_dist, f) #

            with open(self.binfile, 'rb') as f: hist_dist = pickle.load(f)

            self._interpolate = hist_dist.sf # sf give the cdf P(>w) instead of P(<w)

        return self._interpolate

    @classmethod
    def eval(cls,w):
        istance=cls()
        interpolant = istance.interpolate()
        interpolated_values = interpolant(w)

        return interpolated_values# if len(interpolated_values)>1 else interpolated_values[0]


def compare_Pw():

    plotting()
    # Download the file from Emanuele Berti's website if it does not exist
    if not os.path.isfile("Pw_single.dat"):
        import urllib
        urllib.urlretrieve('http://www.phy.olemiss.edu/~berti/research/Pw_single.dat', "Pw_single.dat")

    wEm,PwEm=np.loadtxt("Pw_single.dat",unpack=True)
    wmy = np.linspace(-0.1,1.1,1000)
    Pwmy=Pomega.eval(wmy)
    f, ax = plt.subplots(2, sharex=True)

    ax[0].plot(wEm,PwEm,label="Davide")
    ax[0].plot(wmy,Pwmy,label="Emanuele")
    Pwmy = Pomega.eval(wEm)
    ax[1].plot(wEm,np.abs(PwEm-Pwmy),c='C2')#*2./(PwEm+Pwmy))
    ax[1].set_xlabel('$\omega$')
    ax[0].set_ylabel('$P(\omega)$')
    ax[1].set_ylabel('Residuals')
    ax[0].legend()
    plt.savefig(sys._getframe().f_code.co_name+".pdf",bbox_inches='tight')


@singleton
class detprob(object):

    def __init__(self,  approximant='IMRPhenomD',
                        psd='aLIGOZeroDetHighPower',
                        f_lower=10.,
                        delta_f=1./40.,
                        snr_threshold=8.,
                        mc1d=int(1000),
                        binfile='Pdetint.pkl',
                        screen=False):

        self.approximant=approximant
        self.f_lower=f_lower
        self.delta_f=delta_f
        self.psd=psd
        self.snr_threshold=snr_threshold
        self.mc1d=mc1d
        assert isinstance(self.mc1d,(int,long))
        self.binfile=binfile
        self.screen=screen

        self._interpolate = None

        #self.m1=m1 # Solar mass
        #self.m2=m2 # Solar mass
        #self.z=z   # Redshift
        #self.lumdist = astropy.cosmology.Planck15.luminosity_distance(z).value  # Mpc


    def snr(self,m1_vals,m2_vals,z_vals):

        if not hasattr(m1_vals, "__len__"): m1_vals=[m1_vals]
        if not hasattr(m2_vals, "__len__"): m2_vals=[m2_vals]
        if not hasattr(z_vals, "__len__"): z_vals=[z_vals]

        snr=[]
        for m1,m2,z in zip(m1_vals,m2_vals,z_vals):
            lum_dist = astropy.cosmology.Planck15.luminosity_distance(z).value # luminosity distance in Mpc
            hp, hc = pycbc.waveform.get_fd_waveform(approximant=self.approximant,
                                                    mass1=m1*(1.+z),
                                                    mass2=m2*(1.+z),
                                                    delta_f=self.delta_f,
                                                    f_lower=self.f_lower,
                                                    distance=lum_dist)
            evaluatedpsd = pycbc.psd.analytical.from_string(self.psd,len(hp), self.delta_f, self.f_lower)
            snr_one=pycbc.filter.sigma(hp, psd=evaluatedpsd, low_frequency_cutoff=self.f_lower)
            if self.screen==True:
                print("  m1="+str(m1)+" m1="+str(m2)+" z="+str(z)+" SNR="+str(snr_one))
            snr.append(snr_one ) # use hp only because I want optimally oriented sources

        return np.array(snr) if len(snr)>1 else snr[0]


    def compute(self,m1,m2,z):
        snr = self.snr(m1,m2,z)
        return Pomega.eval(self.snr_threshold/snr)



    def interpolate(self):
        if self._interpolate is None:

            # Takes some time. Store a pickle...
            if not os.path.isfile(self.binfile):

                print('['+self.__class__.__name__+'] Interpolating...')

                # See https://stackoverflow.com/a/30059599
                m1_grid = np.linspace(1,100,self.mc1d)
                m2_grid = np.linspace(1,100,self.mc1d)
                z_grid  = np.linspace(1e-4,3,self.mc1d)

                values = np.zeros([len(m1_grid),len(m2_grid),len(z_grid)])

                for i,m1 in enumerate(m1_grid):
                    for j,m2 in enumerate(m2_grid):
                        for k,z in enumerate(z_grid):
                            values[i,j,k] = self.compute(m1,m2,z)

                interpolant = scipy.interpolate.RegularGridInterpolator(points=points,values=values,fill_value=None)

                with open(self.binfile, 'wb') as f: pickle.dump(interpolant, f)

            with open(self.binfile, 'rb') as f: interpolant = pickle.load(f)

            self._interpolate = interpolant

        return self._interpolate

    @classmethod
    def eval(cls,m1,m2,z):

        if not hasattr(m1, "__len__"): m1=[m1]
        if not hasattr(m2, "__len__"): m2=[m2]
        if not hasattr(z, "__len__"): z=[z]

        istance=cls()
        interpolant = istance.interpolate()

        interpolated_values = interpolant(np.transpose([m1,m2,z]))

        return interpolated_values if len(interpolated_values)>1 else interpolated_values[0]

dp=detprob(screen=True)
print dp.compute(10,10,0.1)
print dp.eval(m1=10,m2=10,z=0.1)
print dp.eval(m1=[10,10],m2=[10,10],z=[0.1,0.1])






#
#


#
#
#
# def compute_my_detprob():
#
#
#
#
# def Pwfloor(w,interp):
#     '''Return 0 if SNR is below treshold.'''
#     if w>=1. :
#         return 0.
#     elif w<=0. :
#         return 1.
#     else:
#         return float(interp(w))
#
#
# def getSNR(m1,m2,z):
#     ''' Compute the SNR using pycbc'''
#
#     lum_dist = cosmo.luminosity_distance(z).value # luminosity distance in Mpc
#     print(z,lum_dist)
#     # get waveform
#     flow=10
#     hp, hc = pycbc.waveform.get_fd_waveform(approximant='IMRPhenomD',mass1=m1*(1.+z),mass2=m2*(1.+z),delta_f=1./40.,f_lower=flow,distance=lum_dist)
#     psd = pycbc.psd.aLIGOZeroDetHighPower(len(hp), hp.delta_f, flow)
#     snr = pycbc.filter.sigma(hp, psd=psd, low_frequency_cutoff=flow) # use hp only because I want optimally oriented sources
#
#     return snr
#
# def compute_detected_probability(snr,SNRthr=8):
#
#     pfint= Pwint() # Interpolate the peanut factor function  for all
#     #detprob= np.array([Pwfloor(SNRthr/x,pfint) for x in snr]) # Compute dection probabilities
#     detprob= Pwfloor(SNRthr/snr,pfint)# Compute dection probabilities
#
#     return detprob
