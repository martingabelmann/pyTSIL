# pyTSIL

(C)Python interface to [TSIL (Stephen P. Martin and David G. Robertson)](https://www.niu.edu/spmartin/TSIL/) [CPC 174
(2006) 133-151](https://arxiv.org/abs/hep-ph/0501132)].

## Install

 * download and build [TSIL](https://www.niu.edu/spmartin/TSIL/)
  * edit the Makefile to build with `-fPIC`
 * build: `TSILDIR=/home/user/tsil python3 setup.py build && sudo -E python3 setup.py  install`
   * `TSILDIR` should be the location where TSIL has been installed
## Usage

```
from pyTSIL import TSIL
results = TSIL(1,2,3,4,5,10,1) # yields all one- and two-loop functions computed by TSIL
results["Mxyzuv"] # gives M(x=1,y=2,z=3,u=4,v=5,s=10,Q**2=1)
```
   
