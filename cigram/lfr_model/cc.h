int common_neighbors(int a, int b, deque<set<int> > & en);
double compute_cc(deque<set<int> > & en, int i);
double compute_cc(deque<set<int> > & en);
double compute_tot_t(deque<set<int> > & en);
int choose_the_least(deque<set<int> > & en, deque<int> & A, int a, int & cn_a_o);

int cclu(deque<set<int> > & en, const deque<deque<int> > & member_list,
            const deque<deque<int> > & member_matrix, double ca);
