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