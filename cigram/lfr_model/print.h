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