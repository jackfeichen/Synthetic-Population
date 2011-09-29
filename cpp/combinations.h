#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

using std::copy;
using std::cout;
using std::endl;
using std::vector;

#ifndef COMBINATIONS_H
#define	COMBINATIONS_H

size_t choose(size_t n, size_t k);

void combinations_recursive(const vector<short> &elems, 
                            vector<vector<short> > &ilst,
                            size_t req_len, vector<size_t> &pos, 
                            size_t depth, size_t margin);

void combinations(const vector<short> &elems, 
                  vector<vector<short> > &ilst,
                  size_t req_len);

typedef std::vector<bool> V;
typedef std::vector<V> Vv;

// A vector of iterators
// which iterate over your individual vector<int>s.
struct Digits {
    V::const_iterator begin;
    V::const_iterator end;
    V::const_iterator me;
};

typedef std::vector<Digits> Vd;

void cart_product(
    Vv& out,  // final result
    Vv& in);  // final result

#endif	/* COMBINATIONS_H */

