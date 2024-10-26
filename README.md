## MicroscPSF-Py ##

This is a Python implementation of the fast microscope PSF generation tool (using the Gibson-Lanni model).
Technical details can be found [here](http://jizhou.li/project/microsc_psf).

[![PyPI version](https://badge.fury.io/py/MicroscPSF-Py.svg)](https://badge.fury.io/py/MicroscPSF-Py)

### Install ###

#### PyPI ####

```
$ python -m pip install MicroscPSF-Py
````

#### Source ####

```
$ git clone https://github.com/MicroscPSF/MicroscPSF-Py.git
$ cd MicroscPSF-Py
$ python setup.py install
```

### Usage ###

Please see the [examples.ipynb](https://github.com/MicroscPSF/MicroscPSF-Py/blob/master/examples.ipynb) Jupyter notebook.

### Acknowledgements ###

- Original algorithm Li, J., Xue, F., & Blu, T. (2017). Fast and accurate three-dimensional point spread function computation for fluorescence microscopy. JOSA A, 34(6), 1029-1034. [link](https://doi.org/10.1364/JOSAA.34.001029).

- Python implementation by Kyle Douglass and Hazen Babcock.
