#define unlikely -214741
#include <math.h>
#include <iostream>
#include <deque>
#include <set>
#include <vector>
#include <map>
#include <string>
#include <iterator>
#include <algorithm>
#include <Python.h>

using namespace std;

enum logtypes {INFO, WARNING, ERROR, DEBUG};
void log_msg(int type, char const *msg);

#include "benchm.h"
#include "random.h"
#include "combinatorics.h"
#include "cast.h"
#include "cc.h"

