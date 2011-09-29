/* 
 * File:   algorithm.h
 * Author: jackchen
 *
 * Created on September 17, 2011, 7:20 PM
 */

#include <iostream>
#include <iterator>
#include <map>
#include <math.h>
#include <numeric>
#include <stdio.h>
#include <string.h>
#include <ctime>
#include <vector>

using std::map;
using std::ostream;
using std::ostream_iterator;
using std::string;
using std::vector;

#ifndef ALGORITHM_H
#define	ALGORITHM_H

class Algorithm {
    static short nwise_;
    static float mweight_, iweight_;
    unsigned int seed_;
    void Init(int n, int m, unsigned int rndSeed);
    void Build();
    Algorithm() { }
    
public:
    class Initializer {
    public:
        Initializer();
    };
    friend class Initializer;
    static Initializer initializer;

    static vector<vector<bool> > ctable_;
    static map<vector<bool>, vector<bool> > converter_;
    static map<vector<bool>, short> matcher_;
    
    int n_, m_;
    vector<vector<bool> > pop_, sample_;
    vector<vector<short> > ilst_;
    vector<float> expected_marginals_;
    map<vector<short>, float> expected_interactions_;
    
    Algorithm(int n, int m, int rndSeed);
    ~Algorithm();
    
    map<vector<bool>, vector<bool> > GetConverter();
    
    static float Extrapolate(const vector<vector<bool> > &a, 
                             short i, short j, vector<int> &cells,
                             map<vector<bool>, vector<bool> > &converter);
    static float GetInteraction(vector<int> &cells);
    
    float Evaluate(ostream& out);
    static float Evaluate(const vector<vector<bool> > &pop,
                          vector<float> &marginals,
                          map<vector<short>, float> &interactions,
                          vector<vector<short> > &ilst,
                          map<vector<bool>, vector<bool> > &converter,
                          ostream& out);

    void PrintMembers(int start, int end);
    void PrintMembers(int start, int end, ostream& out);
    template <typename T> 
    static void PrintVector(const vector<T> &v, string separator, ostream& out);
    static float Logit(float x);
    static float Error(float x, float y);
    static float URand();
};

#endif	/* ALGORITHM_H */

