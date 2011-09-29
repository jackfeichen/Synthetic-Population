/* 
 * File:   main.cpp
 * Author: jackchen
 *
 * Created on September 17, 2011, 6:43 PM
 */

#include <cstdlib>
#include <iostream>
#include <ctime>
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

/*
 * Executes the population generation.
 */
void test(Algorithm *alg)
{
    // print interactions
//    for(size_t i=0; i<alg->ilst_.size(); i++)
//        PrintVector(alg->ilst_[i], "");

//    alg->PrintMembers(0, 10);
//    PrintVector(alg->expected_marginals_, " ");
    
    // check converter stuff
//    for(int i=0; i<4; i++) 
//        PrintVector(Algorithm::ctable_[i], "");
//    cout << endl;
//    map<vector<bool>, vector<bool> >::iterator it;
//    for(it=Algorithm::converter_.begin(); it!=Algorithm::converter_.end(); ++it)
//        PrintVector(it->second, "");
//    cout << endl;
//    map<vector<bool>, short>::iterator it;
//    for(it=Algorithm::matcher_.begin(); it!=Algorithm::matcher_.end(); ++it)
//        cout << (it->second) << endl;
    
    // test extrapolate
//    vector<vector<bool> > a;
//    for(int i=10; i<20; i++) {
//        PrintVector(alg->sample_[i], " ");
//        a.push_back(alg->sample_[i]);
//    }
//    vector<int> cells(4);
//    map<vector<bool>, vector<bool> > converter = alg->GetConverter();
//    float result = Algorithm::Extrapolate(a,1,2,cells,converter);
//    cout << "Interaction score: " << result << endl;
    
    // expected interactions
//    vector<vector<short> >::const_iterator iter;
//    for(iter=alg->ilst_.begin(); iter!=alg->ilst_.end(); ++iter)
//        cout << "for interaction " << (*iter)[0] << (*iter)[1] <<": " << alg->expected_interactions_[*iter] << endl;
    
    // get error:
//    float error = alg->Evaluate(cout);
//    cout << "resulting error is: " << error << endl;
    
    // should be 0 error here
//    cout << "Using the sample population, resulting error should be 0.\n";
//    alg->pop_ = alg->sample_;
//    float error = alg->Evaluate(cout);
//    cout << "resulting error is: " << error << endl;
}

void test_random()
{
    for(int i=0; i<10; i++)
        cout << Algorithm::URand() << endl;
}

void test_random_population(Algorithm *alg)
{
    srand(time(NULL));
    cout << "Original expected marginals:\n";
    Algorithm::PrintVector(alg->expected_marginals_, "|", cout);
    cout << "The last 10 members:\n";
    alg->PrintMembers(990,1000);

    // generate some sample population
    for(int i=0; i<alg->n_; i++)
        // for each member of the population, lookup the attribute pr
        for(short s=0; s<alg->m_; s++)
            alg->pop_[i][s] = 
                (alg->expected_marginals_[s] > Algorithm::URand()) 
                    ? true : false;
    float error = alg->Evaluate(cout);
    cout << "sum square error is: " << error << endl;

    cout << "The last 10 members:\n";
    alg->PrintMembers(990,1000);
}

void generate_choices(vector<vector<bool> > &choices, int count)
{
    vector<bool> a(2);
    a[1] = true;
    vector<vector<bool> > as(count, a);
    cart_product(choices, as);
}

void test_random_assignment(Algorithm *alg)
{
    cout << "Original expected marginals:\n";
    Algorithm::PrintVector(alg->expected_marginals_, "|", cout);

    // generate 1024 possible combinations based on the attributes
    vector<vector<bool> > choices;
    generate_choices(choices, alg->m_);
    int index, count = choices.size();
    
    // for each member of the population, randomly choose a member
     srand(time(NULL));
    for(size_t i=0; i<alg->n_; i++) {
        index = rand() % count;
        alg->pop_[i] = choices[index];
    }
    
    float error = alg->Evaluate(cout);
    cout << "sum square error is: " << error << endl;

    cout << "The last 10 members:\n";
    alg->PrintMembers(990,1000);
}

/*
 * Main entry into the application.
 */
int main(int argc, char** argv) 
{
    clock_t start = clock();
    Algorithm *alg;
    alg = new Algorithm(1000, 10, 10);

//    test_random();
//    test(alg);
//    test_random_population(alg);
    test_random_assignment(alg);
    
    cout.setf(ios::fixed,ios::floatfield);
    cout << "Time elapsed(s):" << ((clock() - start)/(double)CLOCKS_PER_SEC) << endl;
    return 0;
}

