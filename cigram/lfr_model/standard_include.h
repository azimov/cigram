#if !defined(STANDARD_INCLUDE_INCLUDED)
#define STANDARD_INCLUDE_INCLUDED

#endif

#define unlikely -214741

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

int benchmark(bool excess, bool defect, int num_nodes, double  average_k, int  max_degree, double  tau, double  tau2,
	double  mixing_parameter, int  overlapping_nodes, int  overlap_membership,
	int  nmin, int  nmax, bool  fixed_range, double ca);


/* random.cpp */

double ran4(bool t, long s);
double ran2(long *idum);
double ran4();

void srand4(void);
void srand5(int rank);
int irand(int n);
void srand_file(void);

int configuration_model(deque<set<int> > & en, deque<int> & degrees);


/* print.cpp */
int cherr();
int cherr(double a);

template <typename uno, typename due>
void prints(pair <uno, due> &sq,  ostream &out);

template <typename uno, typename due>
void prints(pair <uno, due> &sq);

template <typename uno, typename due>
void prints(map <uno, due> &sq,  ostream &out);

template <typename uno, typename due>
void prints(multimap <uno, due> &sq,  ostream &out);

template <typename Seq>
void prints(Seq &sq, ostream &out);

template <typename type_>
void prints(type_ *a, int b);

template<typename T, template<typename> class C>
void printm(C<T>& c, ostream &out);

template <typename uno, typename due>
void prints(map <uno, due> &sq);

template <typename uno, typename due>
void prints(multimap <uno, due> &sq);

template <typename Seq>
void prints(Seq &sq);

template <typename type>
void prints(const deque<type> & sq);

template<typename T, template<typename> class C>
void printm(C<T>& c);


/* histograms.cpp */
template <typename type>
int log_histogram(deque<type> &c, ostream & out, int number_of_bins);

template <typename type>
int histogram (vector <type> &c, ostream & out, int number_of_bins, double b1, double b2);

template <typename type>
int not_norm_histogram_correlated (deque<type> &c, deque<type> &d,
                            ostream & out, int number_of_bins, double b1, double b2);

template <typename type>
int histogram (deque <type> &c, ostream & out, int number_of_bins, double b1, double b2);


template <typename type>
int not_norm_histogram (vector<type> &c, ostream & out, int number_of_bins, double b1, double b2);


template <typename type>
int not_norm_histogram (deque<type> &c, ostream & out, int number_of_bins, double b1, double b2);

int int_histogram (vector <int> &c, ostream & out);

int int_histogram (deque <int> &c, ostream & out);

/* combinatorics.cpp */
template <typename Seq>
double average_func(Seq &sq);

template <typename Seq>
double variance_func(Seq &sq);

template <typename Seq>
double average_pf(Seq &sq);

template <typename Seq>
double variance_pf(Seq &sq);

double log_factorial (int num);
double log_combination (int n, int k);
double binomial(int n, int x, double p);
int binomial_cumulative(int n, double p, deque<double> &cum);

int powerlaw (int n, int min, double tau, deque<double> &cumulative);

int distribution_from_cumulative(const deque<double> &cum, deque<double> &distr);

int cumulative_from_distribution (deque<double> &cum, const deque<double> &distr);

double poisson (int x, double mu);
int shuffle_and_set(int *due, int dim);
int shuffle_s(deque<int> & sq);

template <typename type_>
int shuffle_s(type_ *a, int b);

double compute_r(int x, int k, int kout, int m);

int add_factors (deque<double> & num, deque<double> &den, int  n, int k);

double compute_hypergeometric(int i, int k, int kout, int m);

int random_from_set(set<int> & s);

/* cc.cpp */
int common_neighbors(int a, int b, deque<set<int> > & en);
double compute_cc(deque<set<int> > & en, int i);
double compute_cc(deque<set<int> > & en);
double compute_tot_t(deque<set<int> > & en);
int choose_the_least(deque<set<int> > & en, deque<int> & A, int a, int & cn_a_o);

int cclu(deque<set<int> > & en, const deque<deque<int> > & member_list,
            const deque<deque<int> > & member_matrix, double ca);


/* cast.cpp */
bool cast_string_to_double (string &b, double &h);
int cast_int(double u);
int cast_string_to_char(string &file_name, char *b);