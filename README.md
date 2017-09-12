# gwdet
#### Detectability of gravitational-wave signals from compact binary coalescences

This is a short python module to compute the detectability of gravitational-wave signal from compact binaries averaging over sky-location and source inclination.

#### Installation and checkpoints

You can install this using pip

~~~bash
pip install gwdet
~~~

Dependancies include `numpy`, `scipy`,`matplotlib`,`astropy` and `pathos`, which can all be installed with pip (and will be installed automatically if not present). The LIGO software `lal` and `pycbc` are needed to use this code with values other than the default ones (see [here](https://davidegerosa.com/installlal/) for a short guide). 

If you limit yourself to the default values, I provide some checkpoints files which let you use the code without installing any LIGO software. In any case, even if you have `lal`, dowloading these checkpoints will save you a lot of computational time. When you use the code with defaults parameters for the first time, a message like the following will be printed:

~~~bash
[gwdet] You are using defaults values. You can download this interpolant. Use:
    gwdet.download_defaults(.....)
~~~

To download the checkpoints just execute that line in a python interpreter:

~~~~bash
python
>>>> import gwdet
>>>> gwdet.download_defaults(.....)
~~~~



#### Usage

This code has two classes only, `averageangles` and `detectability`. 

##### Detectability as a function of the projection parameter $\omega$



