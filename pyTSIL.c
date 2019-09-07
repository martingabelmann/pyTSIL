#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <tsil.h>
#include <Python.h>

// TODO catch errors, build exceptions

void addvalue(PyObject *dict, const char *key, TSIL_COMPLEX val)
{
    PyObject *z = PyComplex_FromDoubles(creall(val), cimagl(val));
    PyDict_SetItemString(dict, key, z);
    return; 
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

    //TSIL_PrintInfo ();
    //TSIL_PrintVersion ();
    TSIL_PrintStatus (&result);

    PyObject *mydict = PyDict_New();
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
    { "TSIL", TSIL, METH_VARARGS, "returns a dict containing all one- and two-loop functions" },
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
