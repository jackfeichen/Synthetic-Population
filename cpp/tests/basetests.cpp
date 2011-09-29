/* 
 * File:   algorithmtests.cpp
 * Author: jackchen
 *
 * Created on Sep 19, 2011, 9:40:23 AM
 */

#include <stdlib.h>
#include <iostream>
#include "algorithm.h"
#include "combinations.h"

using namespace std;

/*
 * Print a vector of some type.
 */
template <typename T>
void PrintVector(vector<T> &v, string separator) 
{
    cout.setf(ios::fixed,ios::floatfield);
    for(typename vector<T>::iterator it=v.begin(); it!=v.end(); it++)
        cout << *it << separator;
    cout << endl;
}

void test_combination() {
    cout << "Testing generation of interactions.\n";

    vector<short> elems(10);
    vector<vector<short> > list;
    for(short i=0; i<10; i++) elems[i] = i;
    combinations(elems, list, 2);

    for(size_t i=0; i<list.size(); i++)
        PrintVector(list[i], "");

    cout << "Completed.\n";
}

void test_cartesian() {
    cout << "Testing generation of Cartesian product.\n";
    
    vector<bool> a(2);
    a[1] = true;
    
    vector<vector<bool> > result;
    vector<vector<bool> > as(2, a);
    cart_product(result, as);

    for(size_t i=0; i<result.size(); i++)
        PrintVector(result[i], "");
    
    cout << "Completed.\n";
}

/*
 * Simple C++ Test Suite
 */
int main(int argc, char** argv) {
    clock_t start = clock();
    cout.setf(ios::fixed,ios::floatfield);

    test_combination();
    cout << "Time elapsed(s):" << ((clock() - start)/(double)CLOCKS_PER_SEC) << endl;
    test_cartesian();
    cout << "Time elapsed(s):" << ((clock() - start)/(double)CLOCKS_PER_SEC) << endl;

    return (EXIT_SUCCESS);
}

