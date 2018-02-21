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