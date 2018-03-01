/*
 * CiGRAM - Circular Gaussian Random grAph Model
 * Copyright (C) 2017  James Gilbert
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details:
 * http://www.gnu.org/licenses/gpl.txt
 */
#include <Python.h>

#include "standard_include.h"


struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

// wrapper for python
static PyObject* generate_graph(PyObject* self, PyObject* args)
{

    double  average_k,  tau, tau2, mixing_parameter, clustering;
	int num_nodes, max_degree, overlapping_nodes, overlap_membership, nmin, nmax;
	bool fixed_range, excess, defect;
	long seed;

	// Convert python arguments into c
	if (!PyArg_ParseTuple(args, "idddddiiiiilbbb", &num_nodes, &average_k, &tau, &tau2, &mixing_parameter, &clustering,
	   &max_degree, &overlapping_nodes, &overlap_membership, &nmin, &nmax,  &seed, &excess, &defect, &fixed_range))
		return NULL;

    // Set the seed from python (no need for seed.dat file)
    srand5(seed);

    PyObject *edgeList = PyList_New(0);
	PyObject *communityList = PyList_New(0);

	// build edge list to be returned,
    generate_benchmark(excess, defect, num_nodes, average_k, max_degree, tau, tau2,
	                    mixing_parameter, overlapping_nodes, overlap_membership, nmin, nmax, fixed_range,
	                    clustering, edgeList, communityList);

	// Derefence our very large lists
	PyObject *tup_return = Py_BuildValue("(OO)", edgeList, communityList);
	Py_DECREF(edgeList);
	Py_DECREF(communityList);
	return tup_return;
}

//  define functions in module
static PyMethodDef LFRmodelMethods[] =
{
     {"generate_graph", generate_graph, METH_VARARGS, "generate graph with specified parameters"},
     {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3

static int LFRmodel_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int LFRmodel_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "lfr_model",
        NULL,
        sizeof(struct module_state),
        LFRmodelMethods,
        NULL,
        LFRmodel_traverse,
        LFRmodel_clear,
        NULL
};

#define INITERROR return NULL
// module initialization python 3
PyMODINIT_FUNC
PyInit_lfr_model(void)
{
    PyObject *module = PyModule_Create(&moduledef);
    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("lfr_model.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

    return module;
}
#else
// module initialization python 2
#define INITERROR return
PyMODINIT_FUNC
initlfr_model(void)
{
     (void) Py_InitModule("lfr_model", LFRmodelMethods);
}
#endif