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

#include "generate_graph.h"

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>

using namespace boost;

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
	int N, K, formConnected, minDegree, minCommunitySize, seed;
	double Density, sigmaF, sigmaR, communitySigmaF,  communitySigmaR, a, ek_per, overlapProb;

	// Convert python arguments into
	if (!PyArg_ParseTuple(args, "iiddddddddiiil", &N, &K, &Density, &sigmaF, &sigmaR, &a, &communitySigmaF,
	 &communitySigmaR, &ek_per, &overlapProb, &formConnected, &minDegree, &minCommunitySize, &seed))
		return NULL;
	
	/*
	*	create graph
	* 	create python list object
	* 	iterate through  graph edges
	* 	create python tuple object
	* 	add edge tuple to python list object
	* 	return python list
	*/
	std::vector< std::vector<int> > communities;
	std::vector<double> theta_communities;
	std::vector<double> theta;

	Graph G = generateGraph(N, K, Density, sigmaF, sigmaR, a, communitySigmaF, communitySigmaR, ek_per, overlapProb,
	 theta, communities, theta_communities, formConnected, minDegree, minCommunitySize, seed);
	PyObject *edgeList = PyList_New(0);

	// build edge list to be returned
	graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		PyObject *tup = Py_BuildValue("(ii)", source(*ei, G), target(*ei, G));
		PyList_Append(edgeList, tup);
		Py_DECREF(tup);
	}


	PyObject *communityList = PyList_New(0);
	// build community list of lists to be returned
	for (uint i = 0; i < communities.size(); ++i){
		PyObject *comm = PyList_New(0);
		for (uint j=0; j < communities[i].size(); ++j){
			PyObject *in = Py_BuildValue("i", communities[i][j]);
			PyList_Append(comm, in);
			Py_DECREF(in);
		}
		PyList_Append(communityList, comm);
		Py_DECREF(comm);
	}

	PyObject *node_positions = PyList_New(0);

	for (uint i = 0; i < theta.size(); ++i) {
		PyObject *theta_i = Py_BuildValue("d", theta[i]);
		PyList_Append(node_positions, theta_i);
		Py_DECREF(theta_i);
	}

	PyObject *community_positions = PyList_New(0);

	for (uint i = 0; i < theta_communities.size(); ++i) {
		PyObject *theta_i = Py_BuildValue("d", theta_communities[i]);
		PyList_Append(community_positions, theta_i);
		Py_DECREF(theta_i);
	}

	// Derefence our very large lists
	PyObject *tup_return = Py_BuildValue("(OOOO)", edgeList, communityList, node_positions, community_positions);
	Py_DECREF(edgeList); 
	Py_DECREF(communityList);
	Py_DECREF(node_positions); 
	Py_DECREF(community_positions);
	return tup_return;
}


//  define functions in module
static PyMethodDef CmodelMethods[] =
{
     {"generate_graph", generate_graph, METH_VARARGS, "generate graph with specified parameters"},
     {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3

static int Cmodel_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int Cmodel_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "cmodel",
        NULL,
        sizeof(struct module_state),
        CmodelMethods,
        NULL,
        Cmodel_traverse,
        Cmodel_clear,
        NULL
};

#define INITERROR return NULL
// module initialization python 3
PyMODINIT_FUNC
PyInit_cmodel(void)
{
    PyObject *module = PyModule_Create(&moduledef);
    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("myextension.Error", NULL, NULL);
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
initcmodel(void)
{
     (void) Py_InitModule("cmodel", CmodelMethods);
}
#endif