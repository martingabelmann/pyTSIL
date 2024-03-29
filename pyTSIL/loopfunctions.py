from .TSIL import TSIL
from functools import lru_cache
from cmath import log

global q, qq, REAL

REAL = lambda x: x # noqa: E731

def return_real(val=True):
    """ val=True: all loop functions return only their real part.
        val=False: all loop function return complex numbers.
        By default all loop functions return complex numbers.
    """
    global REAL
    if val:
        REAL = lambda x: x.real # noqa: E731
    else:
        REAL = lambda x: x # noqa: E731
    clear_cache()

def clear_cache():
    """ clears cache of all wrapper functions """
    for func in [
            A, B, dBds,dBdx, C,
            S, I, T, TBAR, U, V, M,
            Aeps, Beps, Ceps, dBdxeps,
            S0, I0, T0, U0, V0, TSILc
            ]:
        info = func.cache_info()
        if info.currsize > 0:
            print(f'clearing cache of {func.__name__} ({info})')
            func.cache_clear()

def set_renscale(qren):
    """ sets the renormalization scale for all loop function wrappers (not the squared one!) """
    global q, qq
    q = qren
    qq = qren**2
    clear_cache()

@lru_cache(maxsize=None)
def TSILc(x,y,z,u,v,s):
    """ wrapper for TSIL, uses globally defined renormalisation scale (see `set_renscale`) and caching """
    return TSIL(x,y,z,u,v,s,qq)


# renormalised loop functions w/o eps pieces
# one-loop
@lru_cache(maxsize=None)
def A(x):
    """ one-loop one-point function """
    # return TSILc(x,x,x,x,x,0)['Ax']
    try:
        return x*(log(x/qq)-1)
    except ValueError:
        return 0

@lru_cache(maxsize=None)
def B(x,y,s):
    """ one-loop two-point function """
    return REAL(TSILc(x,x,y,x,x,s)['Bxz'])

@lru_cache(maxsize=None)
def dBds(x,y,s):
    """ derivative of one-loop two-point function w.r.t. squared external momentum """
    return REAL(TSILc(x,x,y,x,x,s)['dBdsxz'])

@lru_cache(maxsize=None)
def dBdx(x,y,s):
    """ derivative of one-loop two-point function w.r.t. first squared mass argument"""
    return REAL(TSILc(x,x,y,x,x,s)['Bpxz'])

@lru_cache(maxsize=None)
def C(x,y,z,s):
    """ one-loop two-point function  (C0 integral)"""
    return REAL(TSILc(x,y,z,x,x,s)['Cfinxyz'])
# two-loop
@lru_cache(maxsize=None)
def S(x,y,z,s):
    return REAL(TSILc(x,x,x,y,z,s)['Sxuv'])

@lru_cache(maxsize=None)
def I(x,y,z): # noqa: E741, E743
    return REAL(S(x,y,z,0))

@lru_cache(maxsize=None)
def T(x,y,z,s):
    return REAL(TSILc(x,x,x,y,z,s)['Txuv'])

@lru_cache(maxsize=None)
def TBAR(x,y,z,s):
    return REAL(TSILc(x,x,x,y,z,s)['TBARxuv'])

@lru_cache(maxsize=None)
def U(x,y,z,u,s):
    return REAL(TSILc(x,y,y,z,u,s)['Uxzuv'])

@lru_cache(maxsize=None)
def V(x,y,z,u,s):
    return REAL(TSILc(x,y,y,z,u,s)['Vxzuv'])

@lru_cache(maxsize=None)
def M(x,y,z,u,v,s):
    return REAL(TSILc(x,y,z,u,v,s)['Mxyzuv'])

# O(eps) pieces of one-loop functions
@lru_cache(maxsize=None)
def Aeps(x):
    return REAL(TSILc(x,x,x,x,x,0)['Aepsx'])

@lru_cache(maxsize=None)
def Beps(x,y,s):
    return REAL(TSILc(x,x,y,x,x,s)['Bepsxz'])

@lru_cache(maxsize=None)
def Ceps(x,y,z,s):
    return REAL(TSILc(x,y,z,x,x,s)['Cepsxyz'])

@lru_cache(maxsize=None)
def dBdxeps(x,y,s):
    return REAL(TSILc(x,x,y,x,x,s)['Bepsprimexz'])


# two-loop functions w/ eps^0 pieces
@lru_cache(maxsize=None)
def S0(x,y,z,s):
    return REAL(TSILc(y,x,x,x,z,s)['SSuxv0'])

@lru_cache(maxsize=None)
def I0(x,y,z):
    return REAL(S0(x,y,z,0))

@lru_cache(maxsize=None)
def T0(x,y,z,s):
    return REAL(TSILc(x,x,x,y,z,s)['TTxuv0'])

@lru_cache(maxsize=None)
def U0(x,y,z,u,s):
    return REAL(TSILc(x,y,y,z,u,s)['UUxzuv0'])

@lru_cache(maxsize=None)
def V0(x,y,z,u,s):
    return REAL(TSILc(x,y,y,z,u,s)['VVxzuv0'])
