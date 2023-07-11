# pyTSIL

(C)Python interface to [TSIL (Stephen P. Martin and David G. Robertson)](https://www.niu.edu/spmartin/TSIL/) [CPC 174
(2006) 133-151](https://arxiv.org/abs/hep-ph/0501132)].

## Install

 * download and build [TSIL](https://www.niu.edu/spmartin/TSIL/)
   * edit the Makefile to build with `-fPIC`
 * install pyTSIL: 
   * `export TSILDIR=/home/user/tsil`
   * `python3 -m pip install --break-system-packages --user`
   * `TSILDIR` should be the location where TSIL has been installed
## Usage

Using the TSIL c-interface:
```
from pyTSIL import TSIL
results = TSIL(1,2,3,4,5,10,1) # yields all one- and two-loop functions computed by TSIL
results["Mxyzuv"] # gives M(x=1,y=2,z=3,u=4,v=5,s=10,Q**2=1)
```

Using wrapper functions:
```
from pyTSIL import loopfunctions

loopfunctions.set_renscale(172.5)
x,y,z,u,v,s = 1,2,3,4,5,6
loopfunctions.M(x,y,z,u,v,s) # returns result of the M integral (complex number)

loopfunctions.return_real(True) # from here on all loop functions return their real parts only

loopfunctions.dBds(x,y,s)

loopfunctions.return_real(False) # again return complex numbers 

loopfunctions.clear_cache() # clears cache (automatically done by `set_renscale` and `return_real`)
```
All wrapper functions are cached using `lru_cache`.


The package also provides some loop functions required for sub-loop renormalisation which are not included in the original TSIL package.
