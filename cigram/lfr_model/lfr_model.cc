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



int generate_benchmark(bool excess, bool defect, int num_nodes, double  average_k, int  max_degree, double  tau, double  tau2,
	double  mixing_parameter, int  overlapping_nodes, int  overlap_membership, int  nmin, int  nmax, bool  fixed_range,
	double ca, PyObject* edgeList, PyObject* communityList) {

	double dmin=solve_dmin(max_degree, average_k, -tau);
	if (dmin==-1)
		return -1;

	int min_degree=int(dmin);

	double media1=integer_average(max_degree, min_degree, tau);
	double media2=integer_average(max_degree, min_degree+1, tau);

	if (fabs(media1-average_k)>fabs(media2-average_k))
		min_degree++;

	if (!fixed_range) {
		nmax=max_degree;
		nmin=max(int(min_degree), 3);
		// cout<<"-----------------------------------------------------------"<<endl;
		// cout<<"community size range automatically set equal to ["<<nmin<<" , "<<nmax<<"]"<<endl;
	}


	deque <int> degree_seq ;		//  degree sequence of the nodes
	deque <double> cumulative;
	powerlaw(max_degree, min_degree, tau, cumulative);

	for (int i=0; i<num_nodes; i++) {

		int nn=lower_bound(cumulative.begin(), cumulative.end(), ran4())-cumulative.begin()+min_degree;
		degree_seq.push_back(nn);

	}

	sort(degree_seq.begin(), degree_seq.end());

	if(deque_int_sum(degree_seq) % 2!=0)
		degree_seq[max_element(degree_seq.begin(), degree_seq.end()) - degree_seq.begin()]--;

	deque<deque<int> >  member_matrix;
	deque<int> num_seq;
	deque<int> internal_degree_seq;

	// ****** internal_degree and membership

	if(internal_degree_and_membership(mixing_parameter, overlapping_nodes, overlap_membership, num_nodes, member_matrix, excess, defect, degree_seq, num_seq, internal_degree_seq, fixed_range, nmin, nmax, tau2)==-1)
		return -1;

	deque<set<int> > E;					// E is the adjacency matrix written in form of list of edges
	deque<deque<int> > member_list;		// row i cointains the memberships of node i
	deque<deque<int> > link_list;		// row i cointains degree of the node i respect to member_list[i][j]; there is one more number that is the external degree

	// cout<<"building communities... "<<endl;
	if(build_subgraphs(E, member_matrix, member_list, link_list, internal_degree_seq, degree_seq, excess, defect)==-1)
		return -1;

	// cout<<"connecting communities... "<<endl;
	connect_all_the_parts(E, member_list, link_list);

	if(erase_links(E, member_list, excess, defect, mixing_parameter)==-1)
		return -1;

	if(ca!=unlikely) {
		// cout<<"trying to approach an average clustering coefficient ... "<<ca<<endl;
		cclu(E, member_list, member_matrix, ca);
	}

	// cout<<"recording network..."<<endl;
    for (int u=0; u < E.size(); u++) {

		set<int>::iterator itb=E[u].begin();

		while (itb!=E[u].end()) {
		    // Creates a python tuple, appends it to list, dereferences the tuple
		    PyObject *tup = Py_BuildValue("(ii)", u+1, *(itb++)+1);
		    PyList_Append(edgeList, tup);
		    Py_DECREF(tup);
        }
	}

	// building communities
	// This differs from the C++ implementation slightly,
	// The list returned is a tuple for community memberships (vertex_id, membership)
	// This is best handled at the python interface
	for (int i = 0; i < member_list.size(); i++){
		for (int j=0; j < member_list[i].size(); j++) {
			PyObject *tup = Py_BuildValue("(ii)", i+1, member_list[i][j]+1);
		    PyList_Append(communityList, tup);
		    Py_DECREF(tup); // Always dereference python objects!
	    }
	}
	return 0;
}

// wrapper for python
static PyObject* lfr_graph(PyObject* self, PyObject* args)
{

    double  average_k,  tau, tau2, mixing_parameter, ca;
	int num_nodes, max_degree, overlapping_nodes, overlap_membership, nmin, nmax;
	bool fixed_range, excess, defect;
	long seed;

    excess = false;
    defect = false;

	// Convert python arguments into c
	if (!PyArg_ParseTuple(args, "iiddddiiiiibl", &num_nodes, &average_k, &tau, &tau2, &mixing_parameter, &ca,
	   &max_degree, &overlapping_nodes, &overlap_membership, &nmin, &nmax, &fixed_range, &seed))
		return NULL;

    // Set the seed from python (no need for seed.dat file)
    srand5(seed);

    PyObject *edgeList = PyList_New(0);
	PyObject *communityList = PyList_New(0);

	// build edge list to be returned
    generate_benchmark(excess, defect, num_nodes, average_k, max_degree, tau, tau2,
	                    mixing_parameter, overlapping_nodes, overlap_membership, nmin, nmax, fixed_range,
	                    ca, edgeList, communityList);

	// Derefence our very large lists
	PyObject *tup_return = Py_BuildValue("(OO)", edgeList, communityList);
	Py_DECREF(edgeList);
	Py_DECREF(communityList);
	return tup_return;
}

//  define functions in module
static PyMethodDef LFRmodelMethods[] =
{
     {"lfr_graph", lfr_graph, METH_VARARGS, "generate graph with specified parameters"},
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
PyInit_lfrmodel(void)
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
initcmodel(void)
{
     (void) Py_InitModule("lfr_model", LFRmodelMethods);
}
#endif