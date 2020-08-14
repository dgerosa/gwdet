'''
gwdet: detectability of gravitational-wave signals from compact binary coalescences
https://github.com/dgerosa/gwdet
'''

from __future__ import print_function,division
import sys
import os
import contextlib
import io
import urllib
import time
import warnings
import six
from six.moves import cPickle as pickle
import multiprocessing

import requests
import numpy as np
import astropy.cosmology
import scipy.stats
import scipy.interpolate
import pathos.multiprocessing


__author__ = "Davide Gerosa"
__license__ = "MIT"
__version__ = "0.1.1"
__email__ = "dgerosa@caltech.edu"
this_module='gwdet'

@contextlib.contextmanager
def nostdout():
    ''' Locally suppress stoud print. See https://stackoverflow.com/a/2829036'''
    save_stdout = sys.stdout
    sys.stdout = io.BytesIO()
    yield
    sys.stdout = save_stdout

with nostdout(): # pycbc prints a very annyoing new line when imported
    try:
        import pycbc.waveform
        import pycbc.filter
        import pycbc.psd
        has_pycbc=True
    except:
        warnings.warn("Can't import pycbc",ImportWarning)
        has_pycbc=False

def plotting():
    ''' Import stuff for plotting'''

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

# Defaults values
defaults={  'directory' : os.path.dirname(__file__),
            'mcn' : int(1e8),
            'mcbins' : int(1e5),
            'approximant' : 'IMRPhenomD',
            'psd_from_string' : True,
            'psd_string' : 'aLIGOZeroDetHighPower',
            'psd_from_path': False,
            'psd_path': None,
            'flow' : 10.,
            'deltaf' : 1./40.,
            'snrthreshold': 8.,
            'massmin' : 1.,
            'massmax' : 100.,
            'zmin' : 1e-4,
            'zmax' : 2.2,
            'mc1d' : int(200)}




class averageangles(object):
    '''
    Compute the detection probability, averaged over all angles (sky location, polarization, inclination, etc), as a function of the projection parameter w. This is defined in arxiv:9301003, but here we follow the notation of arxiv:1405.7016

    Usage:
        p = averageangles(directory=os.path.dirname(__file__), binfile=None, mcn=int(1e8), mcbins=int(1e5))
        det = p(w) # with 0<=w<=1

    Parameters:
        directory: where checkpoints are stored
        binfile: checkpoint file (if None computed from other kwargs)
        mcn: resolution parameter (number of Monte Carlo samples)
        mcbins: resolution parameter (number of interpolated bins)
        w: projection parameter 0<=w<=1, see arxiv:1405.7016 (can be float or array)

    Returns:
        det: GW detectability (float or array)
    '''


    def __init__(self,  directory=defaults['directory'],
                        binfile=None,
                        mcn=defaults['mcn'],
                        mcbins=defaults['mcbins']):

        self._interpolate = None
        self.mcn = mcn # Size of monte carlo sample
        assert isinstance(self.mcn,int)
        self.mcbins = mcbins # Number of bins before interpolating
        assert isinstance(self.mcbins,int)
        self.directory=directory

        if binfile is None:
            #self.binfileonly = self.__class__.__name__+'_'+'_'.join([x+'_'+str(eval(x)) for x in ['mcn','mcbins']]) +'.pkl'
            self.binfileonly = self.__class__.__name__+'_'+'_'.join([k+'_'+str(v) for k,v in zip(['mcn','mcbins'],[self.mcn,self.mcbins])]) +'.pkl'

        if self.directory=='':
            self.directory='.'
        self.binfile=self.directory+'/'+self.binfileonly

        # True if all values are the default ones
        self.is_default=all( [v==defaults[k] for k,v in zip(['mcn','mcbins'],[self.mcn,self.mcbins])] )

    def montecarlo_samples(self,mcn):
        ''' Sample the w parameters over the various angles'''

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
        ''' Compute interpolation. If available, read from file'''

        if self._interpolate is None:

            # Takes some time. Store a pickle...
            if not os.path.isfile(self.binfile):
                if not os.path.exists(self.directory) and self.directory!='':
                    os.makedirs(self.directory)

                if self.is_default:
                    print('\n['+this_module+'] You are using defaults values. You can download this interpolant. Use:\n')
                    print('curl https://raw.githubusercontent.com/dgerosa/gwdet/master/checkpoints/'+self.binfileonly+'.tar.gz -o '+self.binfile+'.tar.gz; tar -xzvf '+self.binfile+'.tar.gz -C '+self.directory+'; rm '+self.binfile+'.tar.gz \n')

                print('['+this_module+'] Storing: '+self.binfile)

                print('['+this_module+'] Interpolating Pw(w)...')
                hist = np.histogram(self.montecarlo_samples(self.mcn),bins=self.mcbins)
                hist_dist = scipy.stats.rv_histogram(hist)

                with open(self.binfile, 'wb') as f: pickle.dump(hist_dist, f)
            with open(self.binfile, 'rb') as f: hist_dist = pickle.load(f)

            self._interpolate = hist_dist.sf # sf give the cdf P(>w) instead of P(<w)

        return self._interpolate

    def eval(self,w):
        ''' Evaluate the interpolant'''

        interpolant = self.interpolate()
        interpolated_values = interpolant(w)

        return interpolated_values# if len(interpolated_values)>1 else interpolated_values[0]

    def __call__(self,w):
        ''' Evaluate the interpolant'''

        return self.eval(w)


# Workaround to make a class method pickable. Compatbile with singletons
# https://stackoverflow.com/a/40749513
#def snr_pickable(x): return detectability()._snr(x)
#def compute_pickable(x): return detectability()._compute(x)

class detectability(object):
    '''
    Compute the detection probability of a non-spinning compact binary. We follow the notation of arxiv:1405.7016.

    Usage:
        p = detectability(directory=os.path.dirname(__file__), binfile=None, binfilepdet=None, approximant='IMRPhenomD', psd='aLIGOZeroDetHighPower', flow=10., deltaf=1./40., snrthreshold=8., massmin=1., massmax=100., zmin=1e-4, zmax=2.2, mc1d=int(200), mcn=int(1e8), mcbins=int(1e5), parallel=True, screen=False)
        det = p(m1,m2,z)

    Parameters:
        directory: where checkpoints are stored
        binfile: checkpoint file (if None computed from other kwargs)
        binfilepdet: checkpoint file (if None computed from other kwargs)
        approximant: waveform appriximant used to compute SNRs. Available list: pycbc.waveform.waveform.print_fd_approximants()
        psd: power spectral density used to compute SNRs. Available list: pycbc.psd.analytical.get_lalsim_psd_list()
        flow: starting freqiency in SNR calculations
        deltaf: resolution parameter (frequency sampling)
        snrthreshold: minimum detectable signal
        massmin,massax: limits on the component masses in Msun. Interpolated inside, extrapolated outside
        zmin,zmax: limits on the redshift. Interpolated inside, extrapolated outside
        mc1d: resolution parameter (number of grid point per dimension)
        mcn: resolution parameter (number of Monte Carlo samples)
        mcbins: resolution parameter (number of interpolated bins)
        parallel: use parallel runs on all CPUs available
        screen: debug option, prints all SNRs computed
        m1: component mass in Msun (can be float or array)
        m2: component mass in Msun (can be float or array)
        z: redshift (can be float or array)

    Returns:
        det: GW detectability (float or array)
    '''

    def __init__(self,  approximant=defaults['approximant'],
                        psd=defaults['psd'],
                        flow=defaults['flow'],
                        deltaf=defaults['deltaf'],
                        snrthreshold=defaults['snrthreshold'],
                        screen=False,
                        parallel=True,
                        massmin=defaults['massmin'],
                        massmax=defaults['massmax'],
                        zmin=defaults['zmin'],
                        zmax=defaults['zmax'],
                        directory=defaults['directory'],
                        binfile=None,
                        binfilepdet=None,
                        mc1d=defaults['mc1d'],
                        mcn=defaults['mcn'],
                        mcbins=defaults['mcbins'],
                        has_pycbc=has_pycbc
                        ):

        self.approximant=approximant
        self.psd=psd
        self.flow=flow
        self.deltaf=deltaf
        self.snrthreshold=snrthreshold
        self.massmin=massmin
        self.massmax=massmax
        self.zmin=zmin
        self.zmax=zmax
        self.mc1d=mc1d # Size of grid where interpolation is performed
        assert isinstance(self.mc1d,int)

        self.directory=directory
        if binfile is None:
            #self.binfileonly = self.__class__.__name__+'_'+'_'.join([x+'_'+str(eval(x)) for x in ['approximant','psd','flow','deltaf','snrthreshold','massmin','massmax','zmin','zmax','mc1d']]) +'.pkl'
            self.binfileonly = self.__class__.__name__+'_'+'_'.join([k+'_'+str(v) for k,v in zip(['approximant','psd','flow','deltaf','snrthreshold','massmin','massmax','zmin','zmax','mc1d'],[self.approximant,self.psd,self.flow,self.deltaf,self.snrthreshold,self.massmin,self.massmax,self.zmin,self.zmax,self.mc1d])]) +'.pkl'


        if self.directory=='':
            self.directory='.'
        self.binfile=self.directory+'/'+self.binfileonly

        #self.tempfile=directory+'/temp.pkl'

        # Flags
        self.screen=screen
        self.parallel=parallel
        # Initialize a pool to run parallel jobs. See https://stackoverflow.com/a/29030149/4481987
        if self.parallel and not os.path.isfile(self.binfile):
            self.map = pathos.multiprocessing.ProcessingPool(multiprocessing.cpu_count()).imap
        self.has_pycbc=has_pycbc

        # Parameters for pdet
        self.mcn = mcn
        self.mcbins = mcbins
        self.binfilepdet=binfilepdet

        # Lazy loading
        self._interpolate = None
        self._snrinterpolant = None
        self._pdetproj = None

        # True if all values are the default ones
        #self.is_default=all( [eval('self.'+x)==defaults[x] for x in ['approximant','psd','flow','deltaf','snrthreshold','massmin','massmax','zmin','zmax','mc1d']])
        self.is_default=all( [v==defaults[k] for k,v in zip(['approximant','psd','flow','deltaf','snrthreshold','massmin','massmax','zmin','zmax','mc1d'],[self.approximant,self.psd,self.flow,self.deltaf,self.snrthreshold,self.massmin,self.massmax,self.zmin,self.zmax,self.mc1d])] )



    def pdetproj(self):
        ''' A single instance of the pdet class'''

        if self._pdetproj is None:
            self._pdetproj = averageangles(directory=self.directory,binfile=self.binfilepdet,mcn=self.mcn,mcbins=self.mcbins)
        return self._pdetproj

    def snr(self,m1_vals,m2_vals,z_vals):
        ''' Compute the SNR from m1,m2,z '''

        if not hasattr(m1_vals, "__len__"): m1_vals=[m1_vals]
        if not hasattr(m2_vals, "__len__"): m2_vals=[m2_vals]
        if not hasattr(z_vals, "__len__"): z_vals=[z_vals]

        snr=[]
        for m1,m2,z in zip(m1_vals,m2_vals,z_vals):
            lum_dist = astropy.cosmology.Planck15.luminosity_distance(z).value # luminosity distance in Mpc

            assert self.has_pycbc, "pycbc is needed"
            hp, hc = pycbc.waveform.get_fd_waveform(approximant=self.approximant,
                                                    mass1=m1*(1.+z),
                                                    mass2=m2*(1.+z),
                                                    delta_f=self.deltaf,
                                                    f_lower=self.flow,
                                                    distance=lum_dist)
            evaluatedpsd = pycbc.psd.analytical.from_string(self.psd,len(hp), self.deltaf, self.flow)
            snr_one=pycbc.filter.sigma(hp, psd=evaluatedpsd, low_frequency_cutoff=self.flow)
            snr.append(snr_one ) # use hp only because I want optimally oriented sources
            if self.screen==True:
               print("  m1="+str(m1)+" m1="+str(m2)+" z="+str(z)+" SNR="+str(snr_one))

        return np.array(snr) if len(snr)>1 else snr[0]

    def _snr(self,redshiftedmasses):
        ''' Utility method '''

        m1z,m2z = redshiftedmasses

        hp, hc = pycbc.waveform.get_fd_waveform(approximant=self.approximant,
                                                mass1=m1z,
                                                mass2=m2z,
                                                delta_f=self.deltaf,
                                                f_lower=self.flow,
                                                distance=1.)
        evaluatedpsd = pycbc.psd.analytical.from_string(self.psd,len(hp), self.deltaf, self.flow)
        snr=pycbc.filter.sigma(hp, psd=evaluatedpsd, low_frequency_cutoff=self.flow)
        if self.screen==True:
           print("  m1="+str(m1z)+" m1="+str(m2z)+" SNR="+str(snr))

        return snr

    def compute(self,m1,m2,z):
        ''' Direct evaluation of the detection probability'''

        snr = self.snr(m1,m2,z)
        return self.pdetproj().eval(self.snrthreshold/snr)

    def _compute(self,data):
        ''' Utility method '''

        m1,m2,z=data
        ld = astropy.cosmology.Planck15.luminosity_distance(z).value
        snrint = self.snrinterpolant()
        snr = snrint([m1*(1+z),m2*(1+z)])/ld
        Pw= self.pdetproj().eval(self.snrthreshold/snr)

        return Pw

    def snrinterpolant(self):
        ''' Build an interpolation for the SNR '''


        if self._snrinterpolant is None:

            assert self.has_pycbc, 'pycbc is needed'

            # Takes some time. Store a pickle...

            # Takes some time. Store a pickle...
            #if not os.path.isfile(self.tempfile):


            print('['+this_module+'] Interpolating SNR...')

            # See https://stackoverflow.com/a/30059599

            m1z_grid = np.linspace(self.massmin*(1.+self.zmin),self.massmax*(1.+self.zmax),self.mc1d) # Redshifted mass 1
            m2z_grid = np.linspace(self.massmin*(1.+self.zmin),self.massmax*(1.+self.zmax),self.mc1d) # Redshifted mass 1
            grids=[m1z_grid,m2z_grid]

            #meshgrid=np.zeros(reduce(lambda x,y:x*y, [len(x) for x in grids]))

            #print(reduce(lambda x,y:x*y, [len(x) for x in grids]))
            meshgrid=[]
            meshcoord=[]
            for i,m1z in enumerate(m1z_grid):
                for j,m2z in enumerate(m2z_grid):
                        meshcoord.append([i,j])
                        meshgrid.append([m1z,m2z])
            meshgrid=np.array(meshgrid)
            meshcoord=np.array(meshcoord)

            if self.parallel:

                # Shuffle the arrays: https://stackoverflow.com/a/4602224/4481987
                # Useful to better ditribute load across processors
                if True:
                    assert len(meshcoord)==len(meshgrid)
                    p = np.random.permutation(len(meshcoord))
                    meshcoord = meshcoord[p]
                    meshgrid = meshgrid[p]

                #pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
                #meshvalues = pool.imap(snr_pickable, meshgrid)
                #pool.close() # No more work
                #meshvalues = pool.imap(self._snr, meshgrid)

                meshvalues= self.map(self._snr, meshgrid)
                while (True):
                    completed = meshvalues._index
                    if (completed == len(meshgrid)): break
                    print('   [multiprocessing] Waiting for', len(meshgrid)-completed, 'tasks...')
                    time.sleep(1)

                #pool.close()
                #pool.join()

            else:
                meshvalues = map(self._snr, meshgrid)

            #print meshvalues

            valuesforinterpolator = np.zeros([len(x) for x in grids])
            for ij,val in zip(meshcoord,meshvalues):
                i,j=ij
                valuesforinterpolator[i,j]=val

            snrinterpolant = scipy.interpolate.RegularGridInterpolator(points=grids,values=valuesforinterpolator,bounds_error=False,fill_value=None)

            self._snrinterpolant = snrinterpolant

            #    with open(self.tempfile, 'wb') as f: pickle.dump(snrinterpolant, f)

            #with open(self.tempfile, 'rb') as f: self._snrinterpolant = pickle.load(f)

        return self._snrinterpolant


    def interpolate(self):
        ''' Build an interpolation for the detection probability as a function of m1,m2,z'''

        if self._interpolate is None:

            # Takes some time. Store a pickle...
            if not os.path.isfile(self.binfile):
                if not os.path.exists(self.directory) and self.directory!='':
                    os.makedirs(self.directory)

                if self.is_default:
                    print('\n['+this_module+'] You are using defaults values. You can download this interpolant. Use:\n')
                    print('curl https://raw.githubusercontent.com/dgerosa/gwdet/master/checkpoints/'+self.binfileonly+'.tar.gz -o '+self.binfile+'.tar.gz; tar -xzvf '+self.binfile+'.tar.gz -C '+self.directory+'; rm '+self.binfile+'.tar.gz \n')

                assert self.has_pycbc, "pycbc is needed"

                print('['+this_module+'] Storing: '+self.binfile)

                # Make sure the other interpolants are available
                dummy=self.snrinterpolant()
                dummy = self.pdetproj()(0.5)

                print('['+this_module+'] Interpolating Pw(SNR)...')

                # See https://stackoverflow.com/a/30059599
                m1_grid = np.linspace(self.massmin,self.massmax,self.mc1d) # Redshifted mass 1
                m2_grid = np.linspace(self.massmin,self.massmax,self.mc1d) # Redshifted mass 1
                z_grid = np.linspace(self.zmin,self.zmax,self.mc1d) # Redshifted mass 1

                grids=[m1_grid,m2_grid,z_grid]

                meshgrid=[]
                meshcoord=[]
                for i,m1 in enumerate(m1_grid):
                    for j,m2 in enumerate(m2_grid):
                            for k,z in enumerate(z_grid):
                                #print i,j,k
                                meshcoord.append([i,j,k])
                                meshgrid.append([m1,m2,z])
                meshgrid=np.array(meshgrid)
                meshcoord=np.array(meshcoord)

                if self.parallel:

                    # Shuffle the arrays: https://stackoverflow.com/a/4602224/4481987
                    # Useful to better ditribute load across processors
                    if True:
                        assert len(meshcoord)==len(meshgrid)
                        p = np.random.permutation(len(meshcoord))
                        meshcoord = meshcoord[p]
                        meshgrid = meshgrid[p]

                    #pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
                    #meshvalues = pool.imap(compute_pickable, meshgrid)
                    #pool.close() # No more work

                    #pool2 = pathos.multiprocessing.ProcessingPool(multiprocessing.cpu_count())
                    #meshvalues = pool2.imap(self._compute, meshgrid)
                    #pool2.close()

                    meshvalues= self.map(self._compute, meshgrid)



                    while (True):
                        completed = meshvalues._index
                        if (completed == len(meshgrid)): break
                        print('   [multiprocessing] Waiting for', len(meshgrid)-completed, 'tasks...')
                        time.sleep(1)
                    #pool.close()

                else:
                    meshvalues = map(self._compute, meshgrid)


                valuesforinterpolator = np.zeros([len(x) for x in grids])
                for ijk,val in zip(meshcoord,meshvalues):
                    i,j,k=ijk
                    valuesforinterpolator[i,j,k]=val

                interpolant = scipy.interpolate.RegularGridInterpolator(points=grids,values=valuesforinterpolator,bounds_error=False,fill_value=None)

                with open(self.binfile, 'wb') as f: pickle.dump(interpolant, f)

            with open(self.binfile, 'rb') as f: interpolant = pickle.load(f)

            self._interpolate = interpolant

        return self._interpolate


    def eval(self,m1,m2,z):
        ''' Evaluate the interpolant'''

        if not hasattr(m1, "__len__"): m1=[m1]
        if not hasattr(m2, "__len__"): m2=[m2]
        if not hasattr(z, "__len__"): z=[z]

        interpolant = self.interpolate()

        interpolated_values = interpolant(np.transpose([m1,m2,z]))

        return interpolated_values if len(interpolated_values)>1 else interpolated_values[0]

    def __call__(self,m1,m2,z):
        ''' Evaluate the interpolant'''

        return self.eval(m1,m2,z)


def compare_Pw():
    ''' Compare performance of the averageangles interpolator against public data from Emanuele Berti's website'''

    plotting()# Initialized plotting stuff

    # Download file from Emanuele Berti's website if it does not exist
    if not os.path.isfile("Pw_single.dat"):
        urllib.urlretrieve('http://www.phy.olemiss.edu/~berti/research/Pw_single.dat', "Pw_single.dat")

    wEm,PwEm=np.loadtxt("Pw_single.dat",unpack=True)
    wmy = np.linspace(-0.1,1.1,1000)
    p=averageangles()
    Pwmy=p(wmy)
    f, ax = plt.subplots(2, sharex=True)

    ax[0].plot(wEm,PwEm,label="public data")
    ax[0].plot(wmy,Pwmy,ls='dashed',label="this code")
    Pwmy = p(wEm)
    ax[1].plot(wEm,np.abs(PwEm-Pwmy),c='C2')#*2./(PwEm+Pwmy))
    ax[1].set_xlabel('$\omega$')
    ax[0].set_ylabel('$P(>\omega)$')
    ax[1].set_ylabel('Residuals')
    ax[0].legend()
    plt.savefig(sys._getframe().f_code.co_name+".pdf",bbox_inches='tight')


def compare_Psnr():
    ''' Evaluate performace of the detectability interpolator against raw SNRs calculations '''

    plotting() # Initialized plotting stuff

    computed=[]
    interpolated=[]

    n=1000

    m1=np.random.uniform(1,100,n)
    m2=np.random.uniform(1,100,n)
    z=np.random.uniform(1e-4,2.5,n)

    p=detectability()
    computed=p.compute(m1,m2,z)
    interpolated=p(m1,m2,z)

    computed=np.array(computed)
    interpolated=np.array(interpolated)
    residuals=np.abs(computed-interpolated)

    residuals_notzero= residuals[np.logical_and(computed!=0,interpolated!=0)]

    f, ax = plt.subplots(2)

    #ax[0].hist(residuals,weights=[1./n for x in residuals],alpha=0.8,bins=100)
    ax[0].hist(residuals_notzero,weights=[1./n for x in residuals_notzero],alpha=0.8,bins=100)
    ax[0].axvline(np.median(residuals_notzero),ls='dashed',c='red')
    #ax[1].hist(residuals,range=[0,0.0001],weights=[1./n for x in residuals],alpha=0.8,bins=100)
    ax[1].hist(residuals_notzero,range=[0,0.0001],weights=[1./n for x in residuals_notzero],alpha=0.8,bins=100)
    ax[1].axvline(np.median(residuals_notzero),ls='dashed',c='red')
    ax[1].set_xlabel('Residuals $P_{\\rm det}$')

    plt.savefig(sys._getframe().f_code.co_name+".pdf",bbox_inches='tight')
