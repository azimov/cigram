#include <Python.h>

#include <math.h>
#include <iostream>
#include <deque>
#include <set>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <ctime>
#include <iterator>
#include <algorithm>

using namespace std;

/* benchm.cpp functions */
bool they_are_mate(int a, int b, const deque<deque<int> > & member_list);
int deque_int_sum(const deque<int> & a);
double integral (double a, double b);
double average_degree(const double &dmax, const double &dmin, const double &gamma);

double solve_dmin(const double& dmax, const double &dmed, const double &gamma);
double integer_average (int n, int min, double tau);

int change_community_size(deque<int> &seq);
int build_bipartite_network(deque<deque<int> >  & member_matrix,
 const deque<int> & member_numbers, const deque<int> &num_seq);

int internal_degree_and_membership (double mixing_parameter, int overlapping_nodes, int max_mem_num,
    int num_nodes, deque<deque<int> >  & member_matrix,
    bool excess, bool defect, deque<int> & degree_seq,
    deque<int> &num_seq, deque<int> &internal_degree_seq,
     bool fixed_range, int nmin, int nmax, double tau2);

int compute_internal_degree_per_node(int d, int m, deque<int> & a);

int build_subgraph(deque<set<int> > & E, const deque<int> & nodes, const deque<int> & degrees);

int build_subgraphs(deque<set<int> > & E, const deque<deque<int> > & member_matrix, deque<deque<int> > & member_list,
	deque<deque<int> > & link_list, const deque<int> & internal_degree_seq, const deque<int> & degree_seq,
	const bool excess, const bool defect);

int connect_all_the_parts(deque<set<int> > & E, const deque<deque<int> > & member_list,
    const deque<deque<int> > & link_list);

int internal_kin(deque<set<int> > & E, const deque<deque<int> > & member_list, int i);

int erase_links(deque<set<int> > & E, const deque<deque<int> > & member_list, const bool excess,
                const bool defect, const double mixing_parameter);


int generate_benchmark(bool excess, bool defect, int num_nodes, double  average_k, int  max_degree, double  tau, double  tau2,
	double  mixing_parameter, int  overlapping_nodes, int  overlap_membership, int  nmin, int  nmax, bool  fixed_range,
	double ca, PyObject* edgeList, PyObject* communityList);
