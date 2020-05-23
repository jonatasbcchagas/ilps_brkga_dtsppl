#ifndef BRKGA_H
#define BRKGA_H

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <climits>

#include "data.h"
#include "tsp_solver.h"
#include "non_dominated_set.h"

using namespace std;


// MersenneTwister.h
// Mersenne Twister random number generator -- a C++ class MTRand
// Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
// Richard J. Wagner  v1.1  28 September 2009  wagnerr@umich.edu

// The Mersenne Twister is an algorithm for generating random numbers.  It
// was designed with consideration of the flaws in various other generators.
// The period, 2^19937-1, and the order of equidistribution, 623 dimensions,
// are far greater.  The generator is also fast; it avoids multiplication and
// division, and it benefits from caches and pipelines.  For more information
// see the inventors' web page at
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html

// Reference
// M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
// Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
// Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.

// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// Copyright (C) 2000 - 2009, Richard J. Wagner
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//   3. The names of its contributors may not be used to endorse or promote 
//      products derived from this software without specific prior written 
//      permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// The original code included the following notice:
// 
//     When you use this, send an email to: m-mat@math.sci.hiroshima-u.ac.jp
//     with an appropriate reference to your work.
// 
// It would be nice to CC: wagnerr@umich.edu and Cokus@math.washington.edu
// when you write.

// Not thread safe (unless auto-initialization is avoided and each thread has
// its own MTRand object)

class MTRand {
// Data
public:
    typedef unsigned long uint32;  // unsigned integer type, at least 32 bits
    
    enum { N = 624 };       // length of state vector
    enum { SAVE = N + 1 };  // length of array for save()

protected:
    enum { M = 397 };  // period parameter
    
    uint32 state[N];   // internal state
    uint32 *pNext;     // next value to get from state
    int left;          // number of values left before reload needed

// Methods
public:
    MTRand( const uint32 oneSeed );  // initialize with a simple uint32
    MTRand( uint32 *const bigSeed, uint32 const seedLength = N );  // or array
    MTRand();  // auto-initialize with /dev/urandom or time() and clock()
    MTRand( const MTRand& o );  // copy
    
    // Do NOT use for CRYPTOGRAPHY without securely hashing several returned
    // values together, otherwise the generator state can be learned after
    // reading 624 consecutive values.
    
    // Access to 32-bit random numbers
    uint32 randInt();                     // integer in [0,2^32-1]
    uint32 randInt( const uint32 n );     // integer in [0,n] for n < 2^32
    //double rand();                        // real number in [0,1] -- disabled by rtoso
    //double rand( const double n );        // real number in [0,n]    -- disabled by rtoso
    double randExc();                     // real number in [0,1)
    double randExc( const double n );     // real number in [0,n)
    double randDblExc();                  // real number in (0,1)
    double randDblExc( const double n );  // real number in (0,n)
    double operator()();                  // same as rand53()
    
    // Access to 53-bit random numbers (capacity of IEEE double precision)
    double rand();        // calls rand53() -- modified by rtoso
    double rand53();      // real number in [0,1)
    
    // Access to nonuniform random number distributions
    double randNorm( const double mean = 0.0, const double stddev = 1.0 );
    
    // Re-seeding functions with same behavior as initializers
    void seed( const uint32 oneSeed );
    void seed( uint32 *const bigSeed, const uint32 seedLength = N );
    void seed();
    
    // Saving and loading generator state
    void save( uint32* saveArray ) const;  // to array of size SAVE
    void load( uint32 *const loadArray );  // from such array
    friend std::ostream& operator<<( std::ostream& os, const MTRand& mtrand );
    friend std::istream& operator>>( std::istream& is, MTRand& mtrand );
    MTRand& operator=( const MTRand& o );

protected:
    void initialize( const uint32 oneSeed );
    void reload();
    uint32 hiBit( const uint32 u ) const { return u & 0x80000000UL; }
    uint32 loBit( const uint32 u ) const { return u & 0x00000001UL; }
    uint32 loBits( const uint32 u ) const { return u & 0x7fffffffUL; }
    uint32 mixBits( const uint32 u, const uint32 v ) const
        { return hiBit(u) | loBits(v); }
    uint32 magic( const uint32 u ) const
        { return loBit(u) ? 0x9908b0dfUL : 0x0UL; }
    uint32 twist( const uint32 m, const uint32 s0, const uint32 s1 ) const
        { return m ^ (mixBits(s0,s1)>>1) ^ magic(s1); }
    static uint32 hash( time_t t, clock_t c );
};

// Functions are defined in order of usage to assist inlining

inline MTRand::uint32 MTRand::hash( time_t t, clock_t c )
{
    // Get a uint32 from t and c
    // Better than uint32(x) in case x is floating point in [0,1]
    // Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)
    
    static uint32 differ = 0;  // guarantee time-based seeds will change
    
    uint32 h1 = 0;
    unsigned char *p = (unsigned char *) &t;
    for( size_t i = 0; i < sizeof(t); ++i )
    {
        h1 *= UCHAR_MAX + 2U;
        h1 += p[i];
    }
    uint32 h2 = 0;
    p = (unsigned char *) &c;
    for( size_t j = 0; j < sizeof(c); ++j )
    {
        h2 *= UCHAR_MAX + 2U;
        h2 += p[j];
    }
    return ( h1 + differ++ ) ^ h2;
}

inline void MTRand::initialize( const uint32 seed )
{
    // Initialize generator state with seed
    // See Knuth TAOCP Vol 2, 3rd Ed, p.106 for multiplier.
    // In previous versions, most significant bits (MSBs) of the seed affect
    // only MSBs of the state array.  Modified 9 Jan 2002 by Makoto Matsumoto.
    register uint32 *s = state;
    register uint32 *r = state;
    register int i = 1;
    *s++ = seed & 0xffffffffUL;
    for( ; i < N; ++i )
    {
        *s++ = ( 1812433253UL * ( *r ^ (*r >> 30) ) + i ) & 0xffffffffUL;
        r++;
    }
}

inline void MTRand::reload()
{
    // Generate N new values in state
    // Made clearer and faster by Matthew Bellew (matthew.bellew@home.com)
    static const int MmN = int(M) - int(N);  // in case enums are unsigned
    register uint32 *p = state;
    register int i;
    for( i = N - M; i--; ++p )
        *p = twist( p[M], p[0], p[1] );
    for( i = M; --i; ++p )
        *p = twist( p[MmN], p[0], p[1] );
    *p = twist( p[MmN], p[0], state[0] );
    
    left = N, pNext = state;
}

inline void MTRand::seed( const uint32 oneSeed )
{
    // Seed the generator with a simple uint32
    initialize(oneSeed);
    reload();
}

inline void MTRand::seed( uint32 *const bigSeed, const uint32 seedLength )
{
    // Seed the generator with an array of uint32's
    // There are 2^19937-1 possible initial states.  This function allows
    // all of those to be accessed by providing at least 19937 bits (with a
    // default seed length of N = 624 uint32's).  Any bits above the lower 32
    // in each element are discarded.
    // Just call seed() if you want to get array from /dev/urandom
    initialize(19650218UL);
    register int i = 1;
    register uint32 j = 0;
    register int k = ( N > seedLength ? N : seedLength );
    for( ; k; --k )
    {
        state[i] =
        state[i] ^ ( (state[i-1] ^ (state[i-1] >> 30)) * 1664525UL );
        state[i] += ( bigSeed[j] & 0xffffffffUL ) + j;
        state[i] &= 0xffffffffUL;
        ++i;  ++j;
        if( i >= N ) { state[0] = state[N-1];  i = 1; }
        if( j >= seedLength ) j = 0;
    }
    for( k = N - 1; k; --k )
    {
        state[i] =
        state[i] ^ ( (state[i-1] ^ (state[i-1] >> 30)) * 1566083941UL );
        state[i] -= i;
        state[i] &= 0xffffffffUL;
        ++i;
        if( i >= N ) { state[0] = state[N-1];  i = 1; }
    }
    state[0] = 0x80000000UL;  // MSB is 1, assuring non-zero initial array
    reload();
}

inline void MTRand::seed()
{
    // Seed the generator with an array from /dev/urandom if available
    // Otherwise use a hash of time() and clock() values
    
    // First try getting an array from /dev/urandom
    FILE* urandom = fopen( "/dev/urandom", "rb" );
    if( urandom )
    {
        uint32 bigSeed[N];
        register uint32 *s = bigSeed;
        register int i = N;
        register bool success = true;
        while( success && i-- )
            success = fread( s++, sizeof(uint32), 1, urandom );
        fclose(urandom);
        if( success ) { seed( bigSeed, N );  return; }
    }
    
    // Was not successful, so use time() and clock() instead
    seed( hash( time(NULL), clock() ) );
}

inline MTRand::MTRand( const uint32 oneSeed )
    { seed(oneSeed); }

inline MTRand::MTRand( uint32 *const bigSeed, const uint32 seedLength )
    { seed(bigSeed,seedLength); }

inline MTRand::MTRand()
    { seed(); }

inline MTRand::MTRand( const MTRand& o )
{
    register const uint32 *t = o.state;
    register uint32 *s = state;
    register int i = N;
    for( ; i--; *s++ = *t++ ) {}
    left = o.left;
    pNext = &state[N-left];
}

inline MTRand::uint32 MTRand::randInt()
{
    // Pull a 32-bit integer from the generator state
    // Every other access function simply transforms the numbers extracted here
    
    if( left == 0 ) reload();
    --left;
    
    register uint32 s1;
    s1 = *pNext++;
    s1 ^= (s1 >> 11);
    s1 ^= (s1 <<  7) & 0x9d2c5680UL;
    s1 ^= (s1 << 15) & 0xefc60000UL;
    return ( s1 ^ (s1 >> 18) );
}

inline MTRand::uint32 MTRand::randInt( const uint32 n )
{
    // Find which bits are used in n
    // Optimized by Magnus Jonsson (magnus@smartelectronix.com)
    uint32 used = n;
    used |= used >> 1;
    used |= used >> 2;
    used |= used >> 4;
    used |= used >> 8;
    used |= used >> 16;
    
    // Draw numbers until one is found in [0,n]
    uint32 i;
    do
        i = randInt() & used;  // toss unused bits to shorten search
    while( i > n );
    return i;
}

//inline double MTRand::rand()
//    { return double(randInt()) * (1.0/4294967295.0); }

//inline double MTRand::rand( const double n )
//    { return rand() * n; }

inline double MTRand::randExc()
    { return double(randInt()) * (1.0/4294967296.0); }

inline double MTRand::randExc( const double n )
    { return randExc() * n; }

inline double MTRand::randDblExc()
    { return ( double(randInt()) + 0.5 ) * (1.0/4294967296.0); }

inline double MTRand::randDblExc( const double n )
    { return randDblExc() * n; }

inline double MTRand::rand53()
{
    uint32 a = randInt() >> 5, b = randInt() >> 6;
    return ( a * 67108864.0 + b ) * (1.0/9007199254740992.0);  // by Isaku Wada
}

inline double MTRand::rand()
    { return rand53(); }

inline double MTRand::randNorm( const double mean, const double stddev )
{
    // Return a real number from a normal (Gaussian) distribution with given
    // mean and standard deviation by polar form of Box-Muller transformation
    double x, y, r;
    do
    {
        x = 2.0 * rand53() - 1.0;
        y = 2.0 * rand53() - 1.0;
        r = x * x + y * y;
    }
    while ( r >= 1.0 || r == 0.0 );
    double s = sqrt( -2.0 * log(r) / r );
    return mean + x * s * stddev;
}

inline double MTRand::operator()()
{
    return rand53();
}

inline void MTRand::save( uint32* saveArray ) const
{
    register const uint32 *s = state;
    register uint32 *sa = saveArray;
    register int i = N;
    for( ; i--; *sa++ = *s++ ) {}
    *sa = left;
}

inline void MTRand::load( uint32 *const loadArray )
{
    register uint32 *s = state;
    register uint32 *la = loadArray;
    register int i = N;
    for( ; i--; *s++ = *la++ ) {}
    left = *la;
    pNext = &state[N-left];
}

inline std::ostream& operator<<( std::ostream& os, const MTRand& mtrand )
{
    register const MTRand::uint32 *s = mtrand.state;
    register int i = mtrand.N;
    for( ; i--; os << *s++ << "    " ) {}
    return os << mtrand.left;
}

inline std::istream& operator>>( std::istream& is, MTRand& mtrand )
{
    register MTRand::uint32 *s = mtrand.state;
    register int i = mtrand.N;
    for( ; i--; is >> *s++ ) {}
    is >> mtrand.left;
    mtrand.pNext = &mtrand.state[mtrand.N-mtrand.left];
    return is;
}

inline MTRand& MTRand::operator=( const MTRand& o )
{
    if( this == &o ) return (*this);
    register const uint32 *t = o.state;
    register uint32 *s = state;
    register int i = N;
    for( ; i--; *s++ = *t++ ) {}
    left = o.left;
    pNext = &state[N-left];
    return (*this);
}

// Change log:
//
// v0.1 - First release on 15 May 2000
//      - Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
//      - Translated from C to C++
//      - Made completely ANSI compliant
//      - Designed convenient interface for initialization, seeding, and
//        obtaining numbers in default or user-defined ranges
//      - Added automatic seeding from /dev/urandom or time() and clock()
//      - Provided functions for saving and loading generator state
//
// v0.2 - Fixed bug which reloaded generator one step too late
//
// v0.3 - Switched to clearer, faster reload() code from Matthew Bellew
//
// v0.4 - Removed trailing newline in saved generator format to be consistent
//        with output format of built-in types
//
// v0.5 - Improved portability by replacing static const int's with enum's and
//        clarifying return values in seed(); suggested by Eric Heimburg
//      - Removed MAXINT constant; use 0xffffffffUL instead
//
// v0.6 - Eliminated seed overflow when uint32 is larger than 32 bits
//      - Changed integer [0,n] generator to give better uniformity
//
// v0.7 - Fixed operator precedence ambiguity in reload()
//      - Added access for real numbers in (0,1) and (0,n)
//
// v0.8 - Included time.h header to properly support time_t and clock_t
//
// v1.0 - Revised seeding to match 26 Jan 2002 update of Nishimura and Matsumoto
//      - Allowed for seeding with arrays of any length
//      - Added access for real numbers in [0,1) with 53-bit resolution
//      - Added access for real numbers from normal (Gaussian) distributions
//      - Increased overall speed by optimizing twist()
//      - Doubled speed of integer [0,n] generation
//      - Fixed out-of-range number generation on 64-bit machines
//      - Improved portability by substituting literal constants for long enum's
//      - Changed license from GNU LGPL to BSD
//
// v1.1 - Corrected parameter label in randNorm from "variance" to "stddev"
//      - Changed randNorm algorithm from basic to polar form for efficiency
//      - Updated includes from deprecated <xxxx.h> to standard <cxxxx> forms
//      - Cleaned declarations and definitions to please Intel compiler
//      - Revised twist() operator to work on ones'-complement machines
//      - Fixed reload() function to work when N and M are unsigned
//      - Added copy constructor and copy operator from Salvador Espana

//=====================================================================================================================//

/*
 * Population.h
 *
 * Encapsulates a population of chromosomes represented by a vector of doubles. We don't decode
 * nor deal with random numbers here; instead, we provide private support methods to set the
 * fitness of a specific chromosome as well as access methods to each allele. Note that the BRKGA
 * class must have access to such methods and thus is a friend.
 *
 *  Created on : Jun 21, 2010 by rtoso
 *  Last update: Nov 15, 2010 by rtoso
 *      Authors: Rodrigo Franco Toso <rtoso@cs.rutgers.edu>
 */

class Population {
    template< class Decoder, class RNG >
    friend class BRKGA;

public:
    unsigned getN() const;    // Size of each chromosome
    unsigned getP() const;    // Size of population

    //double operator()(unsigned i, unsigned j) const;    // Direct access to allele j of chromosome i

    // These methods REQUIRE fitness to be sorted, and thus a call to sortFitness() beforehand
    // (this is done by BRKGA, so rest assured: everything will work just fine with BRKGA).
    double getBestFitness() const;            // Returns the best fitness in this population
    double getFitness(unsigned i) const;    // Returns the fitness of chromosome i
    const std::vector< double >& getChromosome(unsigned i) const;    // Returns i-th best chromosome

private:
    Population(const Population& other);
    Population(unsigned n, unsigned p);
    ~Population();

    std::vector< std::vector< double > > population;        // Population as vectors of prob.
    std::vector< std::pair< double, unsigned > > fitness;    // Fitness (double) of a each chromosome

    void sortFitness();                                    // Sorts 'fitness' by its first parameter
    void setFitness(unsigned i, double f);                // Sets the fitness of chromosome i
    std::vector< double >& getChromosome(unsigned i);    // Returns a chromosome

    double& operator()(unsigned i, unsigned j);        // Direct access to allele j of chromosome i
    std::vector< double >& operator()(unsigned i);    // Direct access to chromosome i
};

inline Population::Population(const Population& pop) :
        population(pop.population),
        fitness(pop.fitness) {
}

inline Population::Population(const unsigned n, const unsigned p) :
        population(p, std::vector< double >(n, 0.0)), fitness(p) {
    if(p == 0) { throw std::range_error("Population size p cannot be zero."); }
    if(n == 0) { throw std::range_error("Chromosome size n cannot be zero."); }
}

inline Population::~Population() {
}

inline unsigned Population::getN() const {
    return population[0].size();
}

inline unsigned Population::getP() const {
    return population.size();
}

inline double Population::getBestFitness() const {
    return getFitness(0);
}

inline double Population::getFitness(unsigned i) const {
    return fitness[i].first;
}

inline const std::vector< double >& Population::getChromosome(unsigned i) const {
    return population[ fitness[i].second ];
}

inline std::vector< double >& Population::getChromosome(unsigned i) {
    return population[ fitness[i].second ];
}

inline void Population::setFitness(unsigned i, double f) {
    fitness[i].first = f;
    fitness[i].second = i;
}

inline void Population::sortFitness() {
    sort(fitness.begin(), fitness.end());
}

//inline double Population::operator()(unsigned chromosome, unsigned allele) const {
//    return population[chromosome][allele];
//}

inline double& Population::operator()(unsigned chromosome, unsigned allele) {
    return population[chromosome][allele];
}

inline std::vector< double >& Population::operator()(unsigned chromosome) {
    return population[chromosome];
}

//=====================================================================================================================//

/*
 * BRKGA.h
 *
 * This class encapsulates a Biased Random-key Genetic Algorithm (for minimization problems) with K
 * independent Populations stored in two vectors of Population, current and previous. It supports
 * multi-threading via OpenMP, and implements the following key methods:
 *
 * - BRKGA() constructor: initializes the populations with parameters described below.
 * - evolve() operator: evolve each Population following the BRKGA methodology. This method
 *                      supports OpenMP to evolve up to K independent Populations in parallel.
 *                      Please note that double Decoder::decode(...) MUST be thread-safe.
 *
 * Required hyperparameters:
 * - n: number of genes in each chromosome
 * - p: number of elements in each population
 * - pe: pct of elite items into each population
 * - pm: pct of mutants introduced at each generation into the population
 * - rhoe: probability that an offspring inherits the allele of its elite parent
 *
 * Optional parameters:
 * - K: number of independent Populations
 * - MAX_THREADS: number of threads to perform parallel decoding -- WARNING: Decoder::decode() MUST
 *                be thread-safe!
 *
 * Required templates are:
 * RNG: random number generator that implements the methods below.
 *     - RNG(unsigned long seed) to initialize a new RNG with 'seed'
 *     - double rand() to return a double precision random deviate in range [0,1)
 *     - unsigned long randInt() to return a >=32-bit unsigned random deviate in range [0,2^32-1)
 *     - unsigned long randInt(N) to return a unsigned random deviate in range [0, N] with N < 2^32
 *
 * Decoder: problem-specific decoder that implements any of the decode methods outlined below. When
 *          compiling and linking BRKGA with -fopenmp (i.e., with multithreading support via
 *          OpenMP), the method must be thread-safe.
 *     - double decode(const vector< double >& chromosome) const, if you don't want to change
 *       chromosomes inside the framework, or
 *     - double decode(vector< double >& chromosome) const, if you'd like to update a chromosome
 *
 *  Created on : Jun 22, 2010 by rtoso
 *  Last update: Sep 28, 2010 by rtoso
 *      Authors: Rodrigo Franco Toso <rtoso@cs.rutgers.edu>
 */

template< class Decoder, class RNG >
class BRKGA {
public:
    /*
     * Default constructor
     * Required hyperparameters:
     * - n: number of genes in each chromosome
     * - p: number of elements in each population
     * - pe: pct of elite items into each population
     * - pm: pct of mutants introduced at each generation into the population
     * - rhoe: probability that an offspring inherits the allele of its elite parent
     *
     * Optional parameters:
     * - K: number of independent Populations
     * - MAX_THREADS: number of threads to perform parallel decoding
     *                WARNING: Decoder::decode() MUST be thread-safe; safe if implemented as
     *                + double Decoder::decode(std::vector< double >& chromosome) const
     */
    BRKGA(unsigned n, unsigned p, double pe, double pm, double rhoe, Decoder& refDecoder, RNG& refRNG, unsigned K = 1, unsigned MAX_THREADS = 1);

    /**
     * Destructor
     */
    ~BRKGA();

    /**
     * Resets all populations with brand new keys
     */
    void reset();

    /**
     * Evolve the current populations following the guidelines of BRKGAs
     * @param generations number of generations (must be even and nonzero)
     * @param J interval to exchange elite chromosomes (must be even; 0 ==> no synchronization)
     * @param M number of elite chromosomes to select from each population in order to exchange
     */
    void evolve(unsigned generations = 1);
    
    /**
     * Exchange elite-solutions between the populations
     * @param M number of elite chromosomes to select from each population
     */
    void exchangeElite(unsigned M);

    /**
     * Returns the current population
     */
    const Population& getPopulation(unsigned k = 0) const;

    /**
     * Returns the chromosome with best fitness so far among all populations
     */
    const std::vector< double >& getBestChromosome() const;

    /**
     * Returns the best fitness found so far among all populations
     */
    double getBestFitness() const;

    // Return copies to the internal parameters:
    unsigned getN() const;
    unsigned getP() const;
    unsigned getPe() const;
    unsigned getPm() const;
    unsigned getPo() const;
    double getRhoe() const;
    unsigned getK() const;
    unsigned getMAX_THREADS() const;

private:

    Data const *data;
    // Hyperparameters:
    const unsigned n;    // number of genes in the chromosome
    const unsigned p;    // number of elements in the population
    const unsigned pe;    // number of elite items in the population
    const unsigned pm;    // number of mutants introduced at each generation into the population
    const double rhoe;    // probability that an offspring inherits the allele of its elite parent

    // Templates:
    RNG& refRNG;                // reference to the random number generator
    Decoder& refDecoder;    // reference to the problem-dependent Decoder

    // Parallel populations parameters:
    const unsigned K;                // number of independent parallel populations
    const unsigned MAX_THREADS;        // number of threads for parallel decoding

    // Data:
    std::vector< Population* > previous;    // previous populations
    std::vector< Population* > current;        // current populations

    // Local operations:
    void initialize(const unsigned i);        // initialize current population 'i' with random keys
    void evolution(Population& curr, Population& next);
    bool isRepeated(const std::vector< double >& chrA, const std::vector< double >& chrB) const;
};

template< class Decoder, class RNG >
BRKGA< Decoder, RNG >::BRKGA(unsigned _n, unsigned _p, double _pe, double _pm, double _rhoe,
        Decoder& decoder, RNG& rng, unsigned _K, unsigned MAX) : n(_n), p(_p),
        pe(unsigned(_pe * p)), pm(unsigned(_pm * p)), rhoe(_rhoe),
        refRNG(rng), refDecoder(decoder), K(_K), MAX_THREADS(MAX),
        previous(K, 0), current(K, 0) {

    // Error check:
    using std::range_error;
    if(n == 0) { throw range_error("Chromosome size equals zero."); }
    if(p == 0) { throw range_error("Population size equals zero."); }
    if(pe == 0) { throw range_error("Elite-set size equals zero."); }
    if(pe > p) { throw range_error("Elite-set size greater than population size (pe > p)."); }
    if(pm > p) { throw range_error("Mutant-set size (pm) greater than population size (p)."); }
    if(pe + pm > p) { throw range_error("elite + mutant sets greater than population size (p)."); }
    if(K == 0) { throw range_error("Number of parallel populations cannot be zero."); }

    // Initialize and decode each chromosome of the current population, then copy to previous:
    for(unsigned i = 0; i < K; ++i) {
        // Allocate:
        current[i] = new Population(n, p);

        // Initialize:
        initialize(i);

        // Then just copy to previous:
        previous[i] = new Population(*current[i]);
    }
}

template< class Decoder, class RNG >
BRKGA< Decoder, RNG >::~BRKGA() {
    for(unsigned i = 0; i < K; ++i) { delete current[i]; delete previous[i]; }
}

template< class Decoder, class RNG >
const Population& BRKGA< Decoder, RNG >::getPopulation(unsigned k) const {
    return (*current[k]);
}

template< class Decoder, class RNG >
double BRKGA< Decoder, RNG >::getBestFitness() const {
    double best = current[0]->fitness[0].first;
    for(unsigned i = 1; i < K; ++i) {
        if(current[i]->fitness[0].first < best) { best = current[i]->fitness[0].first; }
    }

    return best;
}

template< class Decoder, class RNG >
const std::vector< double >& BRKGA< Decoder, RNG >::getBestChromosome() const {
    unsigned bestK = 0;
    for(unsigned i = 1; i < K; ++i) {
        if( current[i]->getBestFitness() < current[bestK]->getBestFitness() ) { bestK = i; }
    }

    return current[bestK]->getChromosome(0);    // The top one :-)
}

template< class Decoder, class RNG >
void BRKGA< Decoder, RNG >::reset() {
    for(unsigned i = 0; i < K; ++i) { initialize(i); }
}

template< class Decoder, class RNG >
void BRKGA< Decoder, RNG >::evolve(unsigned generations) {
    if(generations == 0) { throw std::range_error("Cannot evolve for 0 generations."); }

    for(unsigned i = 0; i < generations; ++i) {
        for(unsigned j = 0; j < K; ++j) {
            evolution(*current[j], *previous[j]);    // First evolve the population (curr, next)
            std::swap(current[j], previous[j]);        // Update (prev = curr; curr = prev == next)
        }
    }
}

template< class Decoder, class RNG >
void BRKGA< Decoder, RNG >::exchangeElite(unsigned M) {
    if(M == 0 || M >= p) { throw std::range_error("M cannot be zero or >= p."); }

    for(unsigned i = 0; i < K; ++i) {
        // Population i will receive some elite members from each Population j below:
        unsigned dest = p - 1;    // Last chromosome of i (will be updated below)
        for(unsigned j = 0; j < K; ++j) {
            if(j == i) { continue; }

            // Copy the M best of Population j into Population i:
            for(unsigned m = 0; m < M; ++m) {
                // Copy the m-th best of Population j into the 'dest'-th position of Population i:
                const std::vector< double >& bestOfJ = current[j]->getChromosome(m);

                std::copy(bestOfJ.begin(), bestOfJ.end(), current[i]->getChromosome(dest).begin());

                current[i]->fitness[dest].first = current[j]->fitness[m].first;

                --dest;
            }
        }
    }

    for(int j = 0; j < int(K); ++j) { current[j]->sortFitness(); }
}

template< class Decoder, class RNG >
inline void BRKGA< Decoder, RNG >::initialize(const unsigned i) {

    for(unsigned j = 0; j < p; ++j) {
        for(unsigned k = 0; k < n; ++k) { (*current[i])(j, k) = refRNG.rand(); }
    }

    vector < vector < int > > distance;
    distance.resize(Data::getInstance().numItems+1);

    for(int i = 0; i < Data::getInstance().numItems + 1; ++i) {
        distance[i].resize(Data::getInstance().numItems + 1);
        for(int j = 0; j < Data::getInstance().numItems + 1; ++j) {
            distance[i][j] = Data::getInstance().pickupDistance[i][j] + Data::getInstance().deliveryDistance[i][j];
        }
    }

    TSPSolver tsp;
    pair < int, vector < int > > result = tsp.solve(Data::getInstance().numItems+1, distance);

    int k = 0;
    double allele = 0.0;
    for(; k < Data::getInstance().numItems; ++k) {
        (*current[i])(0, result.second[k+1] - 1) = allele;
        allele += 0.001;
    }

    for(int x = 1; x <= Data::getInstance().numItems; ++x) {
        allele = 0.0;
        for(int y = 0; y < min(x, Data::getInstance().reloadingDepth + 1); ++y) {
            (*current[i])(0, k) = allele;
            allele += 0.001;
            k += 1;
        }
    }
    
    for(int x = 1; x <= Data::getInstance().numItems; ++x) {
        allele = 0.0;
        for(int y = 0; y < min(Data::getInstance().numItems - x + 1, Data::getInstance().reloadingDepth + 1); ++y) {
            (*current[i])(0, k) = allele;
            allele += 0.001;
            k += 1;
        }
    }
    
    // Decode:
    #ifdef _OPENMP
        #pragma omp parallel for num_threads(MAX_THREADS)
    #endif
    for(int j = 0; j < int(p); ++j) {
        current[i]->setFitness(j, refDecoder.decode((*current[i])(j)) );
    }

    // Sort:
    current[i]->sortFitness();
}

template< class Decoder, class RNG >
inline void BRKGA< Decoder, RNG >::evolution(Population& curr, Population& next) {
    // We now will set every chromosome of 'current', iterating with 'i':
    unsigned i = 0;    // Iterate chromosome by chromosome
    unsigned j = 0;    // Iterate allele by allele

    // 2. The 'pe' best chromosomes are maintained, so we just copy these into 'current':
    while(i < pe) {
        for(j = 0 ; j < n; ++j) { next(i,j) = curr(curr.fitness[i].second, j); }

        next.fitness[i].first = curr.fitness[i].first;
        next.fitness[i].second = i;
        ++i;
    }

    // 3. We'll mate 'p - pe - pm' pairs; initially, i = pe, so we need to iterate until i < p - pm:
    while(i < p - pm) {
        // Select an elite parent:
        const unsigned eliteParent = (refRNG.randInt(pe - 1));

        // Select a non-elite parent:
        const unsigned noneliteParent = pe + (refRNG.randInt(p - pe - 1));

        // Mate:
        for(j = 0; j < n; ++j) {
            const unsigned sourceParent = ((refRNG.rand() < rhoe) ? eliteParent : noneliteParent);

            next(i, j) = curr(curr.fitness[sourceParent].second, j);

            //next(i, j) = (refRNG.rand() < rhoe) ? curr(curr.fitness[eliteParent].second, j) :
            //                                      curr(curr.fitness[noneliteParent].second, j);
        }

        ++i;
    }

    // We'll introduce 'pm' mutants:
    while(i < p) {
        for(j = 0; j < n; ++j) { next(i, j) = refRNG.rand(); }
        ++i;
    }

    // Time to compute fitness, in parallel:
    #ifdef _OPENMP
        #pragma omp parallel for num_threads(MAX_THREADS)
    #endif
    for(int i = int(pe); i < int(p); ++i) {
        next.setFitness( i, refDecoder.decode(next.population[i]) );
    }

    // Now we must sort 'current' by fitness, since things might have changed:
    next.sortFitness();
}

template< class Decoder, class RNG >
unsigned BRKGA<Decoder, RNG>::getN() const { return n; }

template< class Decoder, class RNG >
unsigned BRKGA<Decoder, RNG>::getP() const { return p; }

template< class Decoder, class RNG >
unsigned BRKGA<Decoder, RNG>::getPe() const { return pe; }

template< class Decoder, class RNG >
unsigned BRKGA<Decoder, RNG>::getPm() const { return pm; }

template< class Decoder, class RNG >
unsigned BRKGA<Decoder, RNG>::getPo() const { return p - pe - pm; }

template< class Decoder, class RNG >
double BRKGA<Decoder, RNG>::getRhoe() const { return rhoe; }

template< class Decoder, class RNG >
unsigned BRKGA<Decoder, RNG>::getK() const { return K; }

template< class Decoder, class RNG >
unsigned BRKGA<Decoder, RNG>::getMAX_THREADS() const { return MAX_THREADS; }

//=====================================================================================================================//

/*
 * Decoder.h
 *
 * Any decoder must have the format below, i.e., implement the method decode(std::vector< double >&)
 * returning a double corresponding to the fitness of that vector. If parallel decoding is to be
 * used in the BRKGA framework, then the decode() method _must_ be thread-safe; the best way to
 * guarantee this is by adding 'const' to the end of decode() so that the property will be checked
 * at compile time.
 *
 * The chromosome inside the BRKGA framework can be changed if desired. To do so, just use the
 * first signature of decode() which allows for modification. Please use double values in the
 * interval [0,1) when updating, thus obeying the BRKGA guidelines.
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

class Decoder {

public:
    
    double alpha;
    double beta;
    
    NonDominatedSet nds;
    
    Decoder(double _alpha=1.0, double _beta=1.0);
        
    ~Decoder();

    double decode(const std::vector< double >& chromosome, string solutionFileOut = "");
    
};

#endif
