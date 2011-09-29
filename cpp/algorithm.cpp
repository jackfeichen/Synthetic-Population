#include "algorithm.h"
#include "combinations.h"

// Algorithm constructor
Algorithm::Algorithm(int n, int m, int rndSeed)
{
    Init(n, m, rndSeed);
}

Algorithm::~Algorithm()
{
    
}

/*
 * n-wise interactions
 */
short Algorithm::nwise_ = 2;
float Algorithm::mweight_ = 4.0;
float Algorithm::iweight_ = 1.0;
vector<vector<bool> > Algorithm::ctable_;
map<vector<bool>, vector<bool> > Algorithm::converter_;
map<vector<bool>, short> Algorithm::matcher_;
Algorithm::Initializer initializer;

/*
 * Static initialization of Lookup variables
 */
Algorithm::Initializer::Initializer()
{
    // build contingency table
    int size = pow(2, Algorithm::nwise_);
    Algorithm::ctable_.resize(size);
    for(short s=0; s<size; s++) {
        Algorithm::ctable_[s].resize(size);
        Algorithm::ctable_[s][s] = 1;
    }
    
    // Build a Cartesian product
    vector<bool> a(2);
    a[1] = true;
    vector<vector<bool> > keys;
    vector<vector<bool> > as(2, a);
    cart_product(keys, as);

    // build converter and matcher map
    for(short s=0; s<size; s++) {
        Algorithm::converter_[keys[s]] = Algorithm::ctable_[s];
        Algorithm::matcher_[keys[s]] = s;
    }
}

/*
 * Algorithm initialization
 */
void Algorithm::Init(int n, int m, unsigned int rndSeed)
{
    // initialize primary algorithm attributes
    n_ = n;
    m_ = m;
    seed_ = rndSeed;
    pop_.resize(n);
    for(int i=0; i<n; i++)
        pop_[i].resize(m);
    sample_ = pop_;
    
    // initialize interaction combinations
    vector<short> elems(m);
    for(short i=0; i<10; i++) elems[i] = i;
    combinations(elems, ilst_, Algorithm::nwise_);
    
    Build();
}

void Algorithm::Build()
{
    srand(seed_);
    // generate expected marginals
    expected_marginals_.resize(m_);
    for(short s=0; s<m_; s++)
        expected_marginals_[s] = Algorithm::URand();
    // generate some sample population
    for(int i=0; i<n_; i++)
        // for each member of the population, lookup the attribute pr
        for(short s=0; s<m_; s++)
            sample_[i][s] = (expected_marginals_[s] > Algorithm::URand()) 
                            ? true : false;
    // now recalculate the real distribution of the marginals
    vector<int> sum(m_, 0);
    for(short s=0; s<m_; s++) {
        for(int i=0; i<n_; i++)
            sum[s] += sample_[i][s];
        expected_marginals_[s] = sum[s]/((float)n_);
    }
    
    // now extrapolate the interactions
    float interaction;
    vector<int> cells(m_);
    for(vector<vector<short> >::const_iterator it=ilst_.begin();
        it != ilst_.end(); ++it) {
        // do extrapolate (only does pair-wise for now)
        interaction = Extrapolate(sample_, (*it)[0], (*it)[1], cells, Algorithm::converter_);
        expected_interactions_[*it] = interaction;
    }
}

float Algorithm::Extrapolate(const vector<vector<bool> > &a, 
                             short i, short j, vector<int> &cells,
                             map<vector<bool>, vector<bool> > &converter) 
{
    size_t n = a.size();
    vector<bool> m(2);
    for(size_t s=0; s<n; ++s) {
        m[0] = a[s][i];
        m[1] = a[s][j];
        vector<bool> converted = converter[m];
        for(short t=0; t<converted.size(); t++)
            // take the converted value and update cells
            cells[t] += converted[t];
    }
    return GetInteraction(cells);
}

float Algorithm::GetInteraction(vector<int> &cells)
{
    float x, value=0.0;
    // TODO: correctly implement interaction for an NxN table
    for(short s=0; s<4; s++) {
        x = (cells[s] > 0) ? ((float)cells[s]) : 0.001;
        if(s==0 || s==3) value += log(x);
        else if(s==1 || s==2) value -= log(x);
    }
    return Logit(value);
}

map<vector<bool>, vector<bool> > Algorithm::GetConverter()
{
    return Algorithm::converter_;
}

float Algorithm::Evaluate(ostream& out)
{
    return Algorithm::Evaluate(pop_, expected_marginals_, expected_interactions_, 
                               ilst_, Algorithm::converter_, out);
}

/*
 * Evaluates a population given the target marginals and interactions
 * Return the sum-square-error result and outputs to the provided ostream
 */
float Algorithm::Evaluate(const vector<vector<bool> > &pop,
                          vector<float> &marginals,
                          map<vector<short>, float> &interactions,
                          vector<vector<short> > &ilst,
                          map<vector<bool>, vector<bool> > &converter,
                          ostream& out)
{
    size_t n = pop.size();
    size_t m = marginals.size();

    // error results
    vector<float> error;
    float result = 0.0;
    
    // first, evaluate each of the marginals
    vector<float> calculated_marginals(m);
    vector<int> sum(m, 0);
    for(size_t s=0; s<m; s++) {
        for(int i=0; i<n; i++)
            sum[s] += pop[i][s];
        calculated_marginals[s] = sum[s]/((float)n);
    }
    // calculate the marginal error:
    for(size_t s=0; s<m; s++) {
        error.push_back(mweight_*Error(marginals[s], calculated_marginals[s]));
    }
    
    // now extrapolate the interactions
    float interaction;
    vector<int> cells(m);
    map<vector<short>, float> calculated_interactions;
    vector<vector<short> >::const_iterator iter;
    for(iter=ilst.begin(); iter != ilst.end(); ++iter) {
        // do extrapolate (only does pair-wise for now)
        interaction = Extrapolate(pop, (*iter)[0], (*iter)[1], cells, converter);
        calculated_interactions[*iter] = interaction;
    }
    // calculate the interaction error:
    for(iter=ilst.begin(); iter != ilst.end(); ++iter) {
        error.push_back(iweight_*Error(interactions[*iter], 
                                       calculated_interactions[*iter]));
    }
    
    // print stuff out if the ostream is set
    if(out != NULL) {
        out << "expected marginals:\n";
        PrintVector(marginals, "|", out);
        out << "calculated marginals:\n";
        PrintVector(calculated_marginals, "|", out);
        
        for(size_t s=0; s<m; s++) {
            // do marginals first
            out << "error for marginal:" << s << "=" << error[s] << endl;
            result += error[s];
        }
        int index = m;
        for(iter=ilst.begin(); iter != ilst.end(); ++iter) {
            // do interactions next
            out << "error for interaction:" << (*iter)[0] << (*iter)[1] << "=" << error[index] << endl;
            result += error[index++];
        }
    }
    
    // sum of square errors:
    return result;
}

void Algorithm::PrintMembers(int start, int end)
{
    PrintMembers(start, end, cout);
}

void Algorithm::PrintMembers(int start, int end, ostream& out)
{
    if(start < 0 || end > n_) {
        out << "The start and/or end are outside the bounds of the population.\n";
        return;
    }
    for(int i=start; i<end; i++)
        PrintVector(pop_[i], " ", out);
}

template <typename T>
void Algorithm::PrintVector(const vector<T> &v, string separator, ostream& out)
{
//    cout.setf(ios::fixed,ios::floatfield);
    typename vector<T>::const_iterator it;
    for(it=v.begin(); it!=v.end(); it++)
        out << *it << separator;
    out << endl;
}

// Logit function for normalizing the interaction calculation
float Algorithm::Logit(float x)
{
    return 1.0/(1.0 + exp(-x));
}

// Error calculation function for x and y
float Algorithm::Error(float x, float y)
{
    return pow(x-y, 2);
}

float Algorithm::URand()
{
    return ((float) rand() / (RAND_MAX+1)) ;
}