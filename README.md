# gwdet: Detectability of gravitational-wave signals from compact binary coalescences

This is a short python module to compute the detectability of gravitational-wave signal from compact binaries averaging over sky-location and source inclination.

Right now, it can only handle non-spinning systems, will be generalized to spins eventually (any help is welcome!).

The detectability function is defined by Finn and Chernoff in [arxiv:9301003](https://arxiv.org/abs/gr-qc/9301003) as a function of the projection parameter *Theta* (their Eq. 3.31). Here however, we follow the notation of Dominik+ i n [arxiv:1405.7016](https://arxiv.org/abs/1405.7016), which used the projection parameter *w= Theta/4* (such that *0<=w<=1*). The detectability function is the cumilative distribution *P(>w)* of the projection parameter *w* (see their Eq. A.). Basically this gives the probabilty that an elliptically polarized gravitational-wave signal (like that of a binary) will be detected taking into account the antenna pattern of the detector and the inclination of the source. Here we consider a single detector, but note that in the real world we deal with networks. 

If you have a (non-spinning) binary with masses m1 and m2 at redshift z, you first need to compute its signal-to-noise ratio assuming optimal orientation and location *snr_opt* (i.e. the source is face-on, overhead the detector). Then, specify a threshold, say 8, above which you consider the signal detectable. The probabilty that a specific binary will be detected is just *P(w=8/snr_opt)*.



## Installation and checkpoints

You can install this using pip

~~~bash
pip install gwdet
~~~

Dependancies include `numpy`, `scipy`,`matplotlib`,`astropy`,`requests` and `pathos`, which will be installed automatically if not present. The LIGO software `lal` and `pycbc` are needed to use this code with values other than the default ones (see [here](https://davidegerosa.com/installlal/) for a short guide I wrote). 

If you limit yourself to the default values, I provide some checkpoints files which let you use the code without installing any LIGO software. In any case, even if you have `lal`, dowloading these checkpoints will save you a lot of computational time. When you use the code with defaults parameters for the first time, a message like the following will be printed out:

```
[gwdet] You are using defaults values. You can download this interpolant. Use:
    curl ....
```

To download the checkpoints, just execute that whole command starting with `curl`.  The two checkpoint files are currently ~1MB and ~50MB respectively. 



## Usage

This code has two classes only, `averageangles` and `detectability`. You need to first create an instance of each class and then use it.

### averageangles

Compute the detection probability, averaged over all angles (sky location, polarization, inclination, etc), as a function of the projection parameter w. 

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

​    Compute the detection probability of a non-spinning compact binary.

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

## Cheks and performance

Here I first compare the performance of the P(w) interpolator implemented in `averageangles` against public data from [Emanuele Berti's website](http://www.phy.olemiss.edu/~berti/research/). The agreement is excellent and the residuals are just numerical noise. This plot can be generated with

```
gwdet.compare_Pw()
```



![compare_pw](https://user-images.githubusercontent.com/7237041/30337388-fbd830f6-979c-11e7-9c20-9fde063c9f5e.png)

Seconly, I compare the perfomance of the P(m1,m2,z) interpolator of `detectability` against 1000 bruce force SNR computations from `lal`. Altough occasional mismathces of 3% are found, the median residuals are as small as ~1e-5. This plot can be generated with

```
gwdet.compare_Psnr()
```



![compare_psnr](https://user-images.githubusercontent.com/7237041/30341935-c3bd36e8-97ac-11e7-947d-ac06dae3bedb.png)



## Credits

The code is developed and maintained by [Davide Gerosa](https://davidegerosa.com/). Please, report bugs to

```
dgerosa@caltech.edu
```

I am happy to help you out!
