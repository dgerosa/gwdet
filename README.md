# gwdet
#### Detectability of gravitational-wave signals from compact binary coalescences

This is a short python module to compute the detectability of gravitational-wave signal from compact binaries averaging over sky-location and source inclination.

Right now, it can only handle non-spinning systems, will be generalized to spins eventually (any help is welcome!).



## Installation and checkpoints

You can install this using pip

~~~bash
pip install gwdet
~~~

Dependancies include `numpy`, `scipy`,`matplotlib`,`astropy`,`requests` and `pathos`, which will be installed automatically if not present. The LIGO software `lal` and `pycbc` are needed to use this code with values other than the default ones (see [here](https://davidegerosa.com/installlal/) for a short guide I wrote). 

If you limit yourself to the default values, I provide some checkpoints files which let you use the code without installing any LIGO software. In any case, even if you have `lal`, dowloading these checkpoints will save you a lot of computational time. When you use the code with defaults parameters for the first time, a message like the following will be printed out:

~~~bash
[gwdet] You are using defaults values. You can download this interpolant. Use:
    gwdet.download_defaults(.....)
~~~

To download the checkpoints, just execute that line in a python interpreter:

~~~~bash
python
>>>> import gwdet
>>>> gwdet.download_defaults(.....)
~~~~

The two checkpoint files are currently ~5MB and ~200MB respectively.



## Usage

This code has two classes only, `averageangles` and `detectability`. You need to first create an instance of each class and then use it.

### averageangles

Compute the detection probability, averaged over all angles (sky location, polarization, inclination, etc), as a function of the projection parameter w. This is defined in arxiv:9301003, but here we follow the notation of arxiv:1405.7016:

```
p = averageangles(directory='gwdet_data', binfile=None, mcn=int(1e8), mcbins=int(1e5))

det = p(w) # with 0<=w<=1
```

##### **Parameters**:

- ##### `directory`: where checkpoints are stored

- `binfile`: checkpoint file (if `None` computed from other kwargs)

- `mcn`: resolution parameter (number of Monte Carlo samples)

- `mcbins`: resolution parameter (number of interpolated bins)

- ` w`: projection parameter 0<=w<=1, see arxiv:1405.7016 (can be float or array)

##### **Returns**:

- `det`: GW detectability (float or array)



### detectability

​    Compute the detection probability of a non-spinning compact binary. We follow the notation of arxiv:1405.7016.

```
p = detectability('directory'='gwdet_data', binfile=None, binfilepdet=None, approximant='IMRPhenomD', psd='aLIGOZeroDetHighPower', 'flow'=10., 'deltaf'=1./40., 'snrthreshold'=8., 'massmin'=1., 'massmax'=100., 'zmin'=1e-4, 'zmax'=2.2, 'mc1d'=int(200), mcn=int(1e8), mcbins=int(1e5), parallel=True, screen=False)

det = p(m1,m2,z)
```

**Parameters:**

- `directory`: where checkpoints are stored
- `binfile`: checkpoint file (if None computed from other kwargs)
- `binfilepdet`: checkpoint file (if None computed from other kwargs)
- `approximant`: waveform appriximant used to compute SNRs. Available list: `pycbc.waveform.waveform.print_fd_approximants()`
- `psd`: power spectral density used to compute SNRs. Available list: `pycbc.psd.analytical.get_lalsim_psd_list()`
- `flow`: starting freqiency in SNR calculations
- `deltaf`: resolution parameter (frequency sampling)
- `snrthreshold`: minimum detectable signal
- `massmin`,`massax`: limits on the component masses in Msun. Interpolated inside, extrapolated outside
- `zmin`,`zmax`: limits on the redshift. Interpolated inside, extrapolated outside
- `mc1d`: resolution parameter (number of grid point per dimension)
- `mcn`: resolution parameter (number of Monte Carlo samples)
- `mcbins`: resolution parameter (number of interpolated bins)
- `parallel`: use parallel runs on all CPUs available
- `screen`: debug option, prints all SNRs computed
- `m1`: component mass in Msun (can be float or array)
- `m2`: component mass in Msun (can be float or array)
- `z`: redshift (can be float or array)

**Returns:**

- `det`: GW detectability (float or array)







### Credits

The code is developed and maintained by [Davide Gerosa](https://github.com/dgerosa/precession/blob/master/www.davidegerosa.com). Please, report bugs to

```
dgerosa@caltech.edu
```

I am happy to help you out!
