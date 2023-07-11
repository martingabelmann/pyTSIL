#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <tsil.h>
#include <tsil_funcs.h>
#include <Python.h>

// TODO catch errors, build exceptions

void addvalue(PyObject *dict, const char *key, TSIL_COMPLEX val)
{
    PyObject *z = PyComplex_FromDoubles(creall(val), cimagl(val));
    PyDict_SetItemString(dict, key, z);
    return; 
}

TSIL_COMPLEX TSIL_Bepsprime (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{
    TSIL_COMPLEX t1, t2, t3, t4, dt1, dt2, dt3, dt4, s;
    TSIL_COMPLEX log1msox,sqrtdelta,sqrtdelta2,sqrtdelta3,dsqrtdelta;
    TSIL_REAL x, y, sqrtx, sqrty, lnbarx, lnbary;
    x = X;
    y = Y;
    s = TSIL_AddIeps(S);
    lnbarx = TSIL_CLOG(X/QQ);
    if (TSIL_CABS(S) < TSIL_TOL) {
        if (TSIL_FABS(Y) < TSIL_TOL) return (lnbarx-1.0L)/X;
        lnbary = TSIL_CLOG(Y/QQ);
        if (TSIL_FABS(X-Y)/(X+Y) < TSIL_TOL) return lnbarx/(2.0L*X);
        return (2.0L*(-X + Y + X*lnbarx) - Y*lnbarx*lnbarx + Y*( -2.0L + lnbary)*lnbary)/(2.0L*(X-Y)*(X-Y));
    }
    

    if (TSIL_FABS(Y) < TSIL_TOL) {
        if (TSIL_CABS(S-X) < TSIL_TOL) return (lnbarx-2.0L)/X;
        log1msox = TSIL_CLOG(1-s/X);
        return (2.0L*log1msox - 2.0L*lnbarx*log1msox-log1msox*log1msox + 2.0L*TSIL_CLOG(-(X/(s-X))) + 2.0L*TSIL_Dilog(s/(s-X)))/(2.0L*s);
    }
    lnbary = TSIL_CLOG(Y/QQ);
    sqrtx = TSIL_SQRT(x);
    sqrty = TSIL_SQRT(y);
 
    if (TSIL_CABS(S-TSIL_POW(sqrtx+sqrty,2))/(x+y) < TSIL_TOL) {
        return (-8.0L*(sqrtx+sqrty)+4.0L*sqrtx*lnbarx+sqrty*lnbarx*lnbarx - sqrty*(-4.0L+lnbary)*lnbary)/(4.0L*sqrtx*TSIL_POW(sqrtx+sqrty,2));
    }
    if (TSIL_CABS(S-TSIL_POW(sqrtx-sqrty,2))/(x+y) < TSIL_TOL) {
        return (-8.0L*sqrtx+8.0L*sqrty+4.0L*sqrtx*lnbarx-sqrty*lnbarx*lnbarx+sqrty*(-4.0L+lnbary)*lnbary)/(4.0L*sqrtx*TSIL_POW(sqrtx-sqrty,2));
    }

    sqrtdelta=TSIL_CSQRT(S*S - 2*S*x + x*x - 2*S*y - 2*x*y + y*y);
    sqrtdelta2=TSIL_CSQRT(sqrtdelta*sqrtdelta);
    sqrtdelta3=TSIL_CPOW(sqrtdelta*sqrtdelta, 3.0L/2.0L);
    t1=(s-x+y+sqrtdelta)/(2.0L*sqrtdelta); 
    t2=(-s-x+y+sqrtdelta)/(2.0L*sqrtdelta);
    t3=(-s+x+y+sqrtdelta)/(2.0L*x);
    t4=(-s+x+y-sqrtdelta)/(2.0L*x); 

    dsqrtdelta=-((s-x+y)/sqrtdelta);
    dt1=((s+sqrtdelta-x+y)*(s-sqrtdelta2-x+y))/(2.0L*sqrtdelta3);
    dt2=(-s*s+s*(sqrtdelta-sqrtdelta2)-(sqrtdelta2+x-y)*(sqrtdelta-x+y))/(2.0L*sqrtdelta3);
    dt3=(s*(sqrtdelta-x)-(sqrtdelta+x)*(sqrtdelta-x+y))/(2.0L*sqrtdelta*x*x);
    dt4=(s*(sqrtdelta+x)+(sqrtdelta-x)*(sqrtdelta+x-y))/(2.0L*sqrtdelta*x*x);
    
    return 1.0L/4.0L*(
            -4.0L/x + (2.0L*lnbarx)/x 
            + (1.0L/s)*(
                ((x-y)*(lnbarx-lnbary))/x 
                + ((x-y)*(-4.0L+lnbarx+lnbary))/x  
                + (lnbarx-lnbary)*(-4.0L+lnbarx+lnbary)
                + sqrtdelta*(
                    ((dt4*t3-dt3*t4)*(-4.0L+lnbarx+lnbary))/(t3*t4)
                    + (4.0L*dt1*TSIL_CLOG(1-t1))/t1 
                    - (4.0L*dt2*TSIL_CLOG(1-t2))/t2
                    + 2.0L*(-(dt1/t1)+dt2/t2)*(TSIL_CLOG(1-t1)+TSIL_CLOG(1-t2))
                    + 2.0L*(dt1/(-1.0L+t1)+dt2/(-1.0L+t2))*(TSIL_CLOG(t2)-TSIL_CLOG(t1))
                    + (TSIL_CLOG(t4)-TSIL_CLOG(t3))/x
                ) + dsqrtdelta*(
                    2.0L*(TSIL_CLOG(1-t1)+TSIL_CLOG(1-t2))*(-TSIL_CLOG(t1)+TSIL_CLOG(t2))
                    - 1.0L*(-4.0L+lnbarx+lnbary)*(TSIL_CLOG(t3)-TSIL_CLOG(t4))
                    -4.0L*TSIL_Dilog(t1) +4.0L*TSIL_Dilog(t2))
            )
        );
}

TSIL_COMPLEX TSIL_Cfin (TSIL_REAL X, TSIL_REAL Y, TSIL_REAL Z, TSIL_COMPLEX S, TSIL_REAL QQ)
{
    if (TSIL_FABS(Y-Z) < TSIL_TOL) {
        return TSIL_Bp(Z,X,S,QQ);
    } else {
        return (TSIL_B(Y, X, S, QQ) - TSIL_B(Z, X, S, QQ))/(Y-Z);
    }
}

TSIL_COMPLEX TSIL_Ceps (TSIL_REAL X, TSIL_REAL Y, TSIL_REAL Z, TSIL_COMPLEX S, TSIL_REAL QQ)
{
    if (TSIL_FABS(Y-Z) < TSIL_TOL)
        return TSIL_Bepsprime(Z,X,S,QQ);
    else
        return (TSIL_Beps(Y, X, S, QQ) - TSIL_Beps(Z, X, S, QQ))/(Y-Z);
}


static PyObject* TSIL(PyObject* self, PyObject* args)
{
#include "tsil_global.h"
#include "tsil_names.h"
    int i,j,k;
    double x, y, z, u, v, s, qq;
    clock_t t0, t1,t2;
    
    t0 = clock();
    if(!PyArg_ParseTuple(args, "ddddddd", &x, &y, &z, &u, &v, &s, &qq))
        return NULL;

    TSIL_DATA result;
    TSIL_SetParameters (&result, x, y, z, u, v, qq);
    TSIL_Evaluate (&result, s);
    t1 = clock();

    // TODO add verbose flag
    //TSIL_PrintInfo ();
    //TSIL_PrintVersion ();
    //TSIL_PrintStatus (&result);

    PyObject *mydict = PyDict_New();
    PyDict_SetItemString(mydict, "x", PyFloat_FromDouble(x));
    PyDict_SetItemString(mydict, "y", PyFloat_FromDouble(y));
    PyDict_SetItemString(mydict, "z", PyFloat_FromDouble(z));
    PyDict_SetItemString(mydict, "u", PyFloat_FromDouble(u));
    PyDict_SetItemString(mydict, "v", PyFloat_FromDouble(v));
    PyDict_SetItemString(mydict, "s", PyFloat_FromDouble(s));
    PyDict_SetItemString(mydict, "qq", PyFloat_FromDouble(qq));
    addvalue(mydict, "Mxyzuv", TSIL_GetFunction(&result, "M"));
    for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
            if(i<NUM_U_FUNCS && j<NUM_U_PERMS)
                addvalue(mydict, uname[i][j], TSIL_GetFunction(&result, uname[i][j]));
            if(i<NUM_T_FUNCS && j<NUM_T_PERMS) {
                addvalue(mydict, tname[i][j], TSIL_GetFunction(&result, tname[i][j]));
                addvalue(mydict, tbarname[i][j], TSIL_GetFunction(&result, tbarname[i][j]));
            }
            if(i<NUM_S_FUNCS && j<NUM_S_PERMS)
                addvalue(mydict, sname[i][j], TSIL_GetFunction(&result, sname[i][j]));
            if(i<NUM_B_FUNCS && j<NUM_B_PERMS)
                addvalue(mydict, bname[i][j], TSIL_GetFunction(&result, bname[i][j]));
            if(i<NUM_V_FUNCS && j<NUM_V_PERMS)
                addvalue(mydict, vname[i][j], TSIL_GetFunction(&result, vname[i][j]));
        }
    }
    // also add one-loop functions and their derivatives for convenience
    addvalue(mydict, "Ax", TSIL_A(x, qq));
    addvalue(mydict, "Ay", TSIL_A(y, qq));
    addvalue(mydict, "Az", TSIL_A(z, qq));
    addvalue(mydict, "Au", TSIL_A(u, qq));
    addvalue(mydict, "Av", TSIL_A(v, qq));
    addvalue(mydict, "Apx", TSIL_Ap(x, qq));
    addvalue(mydict, "Apy", TSIL_Ap(y, qq));
    addvalue(mydict, "Apz", TSIL_Ap(z, qq));
    addvalue(mydict, "Apu", TSIL_Ap(u, qq));
    addvalue(mydict, "Apv", TSIL_Ap(v, qq));
    addvalue(mydict, "Bpxz", TSIL_Bp(x, z, s, qq));
    addvalue(mydict, "Bpyu", TSIL_Bp(y, u, s, qq));
    addvalue(mydict, "dBdsxz", TSIL_dBds(x, z, s, qq));
    addvalue(mydict, "dBdsyu", TSIL_dBds(y, u, s, qq));
    addvalue(mydict, "Bpxz", TSIL_Bp(x, z, s, qq));

    // boldface functions and eps parts
    addvalue(mydict, "Aepsx", TSIL_Aeps(x, qq));
    addvalue(mydict, "Aepsy", TSIL_Aeps(y, qq));
    addvalue(mydict, "Aepsz", TSIL_Aeps(z, qq));
    addvalue(mydict, "Aepsu", TSIL_Aeps(u, qq));
    addvalue(mydict, "Aepsv", TSIL_Aeps(v, qq));
    addvalue(mydict, "Bepsxz", TSIL_Beps(x, z, s, qq));
    addvalue(mydict, "Bepsyu", TSIL_Beps(y, u, s, qq));

    addvalue(mydict, "Bepsprimexy", TSIL_Bepsprime(x, y, s, qq));
    addvalue(mydict, "Bepsprimexz", TSIL_Bepsprime(x, z, s, qq));
    addvalue(mydict, "Bepsprimeyu", TSIL_Bepsprime(y, u, s, qq));

    // C0 function
    addvalue(mydict, "Cfinxyz", TSIL_Cfin(x, y, z, s, qq));
    addvalue(mydict, "Cepsxyz", TSIL_Ceps(x, y, z, s, qq));


    for (k=0; k<3; k++) {
    	for (j=0; j<6; j++) {
            if(j<NUM_U_FUNCS)
    		    addvalue(mydict, uuname[j][k], result.U[j].bold[k]);
            if(j<NUM_V_FUNCS)
    		    addvalue(mydict, vvname[j][k], result.V[j].bold[k]);
            if(j<NUM_S_FUNCS)
    		    addvalue(mydict, ssname[j][k], result.S[j].bold[k]);
            if(j<NUM_T_FUNCS)
    		    addvalue(mydict, ttname[j][k], result.T[j].bold[k]);
    	}
    }
    t2 = clock();

    PyObject *calctime = PyFloat_FromDouble(difftime(t1, t0)/CLOCKS_PER_SEC);
    PyObject *runtime = PyFloat_FromDouble(difftime(t2, t1)/CLOCKS_PER_SEC);

    PyDict_SetItemString(mydict, "CalcTime", calctime);
    PyDict_SetItemString(mydict, "RunTime", runtime);

    return mydict;
}

// function definitions
static PyMethodDef funcs[] = {
    { "TSIL", TSIL, METH_VARARGS, "TSIL(x,y,z,u,v,s,qq) returns a dict containing all one- and two-loop functions." },
    { NULL, NULL, 0, NULL }
};

// module definition
static struct PyModuleDef pyTSIL = {
    PyModuleDef_HEAD_INIT,
    "pyTSIL",
    "Interface to TSIL (hep-ph/0501132)",
    -1,
    funcs
};

PyMODINIT_FUNC PyInit_pyTSIL(void)
{
    return PyModule_Create(&pyTSIL);
}
