#include "combinations.h"

size_t choose(size_t n, size_t k) {
    if (k > n) {
        return 0;
    }
    size_t r = 1;
    for (size_t d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
}

void combinations_recursive(const vector<short> &elems, 
                            vector<vector<short> > &ilst,
                            size_t req_len, vector<size_t> &pos, 
                            size_t depth, size_t margin)
{
    // Have we selected the number of required elements?
    if (depth >= req_len) {
        ilst.push_back(vector<short>());
        for (size_t ii = 0; ii < pos.size(); ++ii) {
            ilst.back().push_back(elems[pos[ii]]);
        }
        return;
    }

    // Are there enough remaining elements to be selected?
    // This test isn't required for the function to be correct, but
    // it can save a good amount of futile function calls.
    if ((elems.size() - margin) < (req_len - depth))
        return;

    // Try to select new elements to the right of the last selected one.
    for (size_t ii = margin; ii < elems.size(); ++ii) {
        pos[depth] = ii;
        combinations_recursive(elems, ilst, req_len, pos, depth + 1, ii + 1);
    }
    return;
}

void combinations(const vector<short> &elems, 
                  vector<vector<short> > &ilst,
                  size_t req_len)
{
    assert(req_len > 0 && req_len <= elems.size());
    vector<size_t> positions(req_len, 0);
    combinations_recursive(elems, ilst, req_len, positions, 0, 0);
}


void cart_product(
    Vv& out,  // final result
    Vv& in)  // final result
{
    Vd vd;

    // Start all of the iterators at the beginning.
    for(Vv::const_iterator it = in.begin();
        it != in.end();
        ++it) {
        Digits d = {(*it).begin(), (*it).end(), (*it).begin()};
        vd.push_back(d);
    }


    while(1) {

        // Construct your first product vector by pulling 
        // out the element of each vector via the iterator.
        V result;
        for(Vd::const_iterator it = vd.begin();
            it != vd.end();
            it++) {
            result.push_back(*(it->me));
        }
        out.push_back(result);

        // Increment the rightmost one, and repeat.

        // When you reach the end, reset that one to the beginning and
        // increment the next-to-last one. You can get the "next-to-last"
        // iterator by pulling it out of the neighboring element in your
        // vector of iterators.
        for(Vd::iterator it = vd.begin(); ; ) {
            // okay, I started at the left instead. sue me
            ++(it->me);
            if(it->me == it->end) {
                if(it+1 == vd.end()) {
                    // I'm the last digit, and I'm about to roll
                    return;
                } else {
                    // cascade
                    it->me = it->begin;
                    ++it;
                }
            } else {
                // normal
                break;
            }
        }
    }
}
