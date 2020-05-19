#include <iostream>
#include <cstdio>
#include <iomanip>
#include <vector>
#include <cstring>
#include <string>
#include <algorithm>
#include <chrono>
#include <fstream>

#include "data.h"
#include "ilp_formulation_1.h"
#include "ilp_formulation_2.h"
#include "brkga.h"

using namespace std;

inline void runBRKGA(const string outputSolutionFileName) {
    
    int chromosomeSize = Data::getInstance().numItems;       
    double _a = 200;
    double _pe = 0.10;
    double _pm = 0.25;
    double _rhoe = 0.70; 
        
    for(int k = 1; k <= Data::getInstance().numItems; ++k) {
        chromosomeSize += min(k, Data::getInstance().reloadingDepth + 1);
    }
    
    for(int k = 1; k <= Data::getInstance().numItems; ++k) {
        chromosomeSize += min(Data::getInstance().numItems - k + 1, Data::getInstance().reloadingDepth + 1);
    }

    const unsigned p = chromosomeSize * _a;  // size of population
    const double pe = _pe;                   // fraction of population to be the elite-set
    const double pm = _pm;                   // fraction of population to be replaced by mutants
    const double rhoe = _rhoe;               // probability that offspring inherit an allele from elite parent
    const unsigned K = 1;                    // number of independent populations
    const unsigned MAXT = 1;                 // number of threads for parallel decoding
    
    // clog << _a << ' ' << p << ' ' << pe << ' ' << pm << ' ' << rhoe << ' ' << stoppingCriteriaFactor << endl;
    
    Decoder decoder;                  // initialize the decoder

    double runtime = 3600.0;
    
    std::vector < std::pair < double, std::vector < double > > > solutions;

    long unsigned rng_seed[] = {
                                    269070,  99470, 126489, 644764, 547617, 642580,  73456, 462018, 858990, 756112, 
                                    701531, 342080, 613485, 131654, 886148, 909040, 146518, 782904,   3075, 974703, 
                                    170425, 531298, 253045, 488197, 394197, 519912, 606939, 480271, 117561, 900952, 
                                    968235, 345118, 750253, 420440, 761205, 130467, 928803, 768798, 640300, 871462, 
                                    639622,  90614, 187822, 594363, 193911, 846042, 680779, 344008, 759862, 661168, 
                                    223420, 959508,  62985, 349296, 910428, 964420, 422964, 384194, 985214,  57575, 
                                    639619,  90505, 435236, 465842, 102567, 189997, 741017, 611828, 699223, 335142, 
                                     52119,  49256, 324523, 348215, 651525, 517999, 830566, 958538, 880422, 390645, 
                                    148265, 807740, 934464, 524847, 408760, 668587, 257030, 751580,  90477, 594476, 
                                    571216, 306614, 308010, 661191, 890429, 425031,  69108, 435783,  17725, 335928
                                };
                                
    int _NUM_EXECUTIONS = 10;

    solutions.resize(_NUM_EXECUTIONS);
    
    ofstream fout(outputSolutionFileName + ".log");
    fout << "        UB            Time(s)" << endl;
    
    using namespace std::chrono;
    high_resolution_clock::time_point startTime = high_resolution_clock::now();
    
    for(unsigned exec = 0; exec < _NUM_EXECUTIONS; ++exec) {
        
        fout << "exec #" << fixed << exec+1 << endl;
        
        MTRand rng(rng_seed[exec]);  // initialize the random number generator
    
        // initialize the BRKGA-based heuristic
        BRKGA < Decoder, MTRand > algorithm(chromosomeSize, p, pe, pm, rhoe, decoder, rng, K, MAXT);
        
        fout << fixed << setw(10) << (int)algorithm.getBestFitness() << "        " << fixed << setw(12) << setprecision(1) << 0.0 << endl;
        int lastUB = (int)algorithm.getBestFitness();
        
        using namespace std::chrono;
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        
        while(1) {
            
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            duration<double> time_span = duration_cast<duration<double> >(t2 - t1);

            if((double)time_span.count() >= runtime) break;

            algorithm.evolve();  // evolve the population for one generation
            
            if(lastUB != (int)algorithm.getBestFitness()) {
                fout << fixed << setw(10) << (int)algorithm.getBestFitness() << "        " << fixed << setw(12) << setprecision(1) << (double)time_span.count() << endl;
                lastUB = (int)algorithm.getBestFitness();
            }
        }
        
        solutions[exec] = std::make_pair(algorithm.getBestFitness(), algorithm.getBestChromosome());
    }
    
    high_resolution_clock::time_point endTime = high_resolution_clock::now();
    duration < double > time_span = duration_cast < duration < double > > (endTime - startTime);
    
    std::pair < double, vector < double > > bestSolution = solutions[0];
    for(unsigned exec = 1; exec < _NUM_EXECUTIONS; ++exec) {
        if(solutions[exec].first < bestSolution.first) {
            bestSolution = solutions[exec];
        }
    }

    fout << endl;
    char tmp[10000];
    sprintf(tmp, "%-20s ", outputSolutionFileName.c_str()); fout << tmp;
    for(unsigned exec = 0; exec < _NUM_EXECUTIONS; ++exec) {
        sprintf(tmp, "%15.0lf ", solutions[exec].first); fout << tmp;
    }    
    sprintf(tmp, "%15.3lf", (double)time_span.count()); fout << tmp << endl;
    
    fout.close();

    // save the best solution found
    decoder.decode(bestSolution.second, outputSolutionFileName + ".sol");
}
 
inline void usage() { 
    clog << "\n       Usage ./dtsppl --approach <approach_name> --pickuparea <pickup_area_file_name> --deliveryarea <delivery_area_file_name> --n <number_of_items> --l <reloading_depth> --h <relocation_cost> --outputsolution <solution_file_name> " << endl;
    exit(0);
}
    
int main(int argc, char **argv) {

    if(argc < 15) usage();
    
    int numItems, reloadingDepth, costForEachRelocate;
    char parameterStr[1000];
    string approachID, pickupAreaFileName, deliveryAreaFileName, outputSolutionFileName;    
    
    int check_parameters = 0;
    for(int i = 1; i < argc; i += 2) { 
        if(strcmp(argv[i], "--approach") == 0) { sscanf(argv[i+1],"%s", parameterStr); approachID = parameterStr; check_parameters += 1; }
        else if(strcmp(argv[i], "--pickuparea") == 0) { sscanf(argv[i+1],"%s", parameterStr); pickupAreaFileName = parameterStr; check_parameters += 1; }
        else if(strcmp(argv[i], "--deliveryarea") == 0) { sscanf(argv[i+1],"%s", parameterStr); deliveryAreaFileName = parameterStr; check_parameters += 1; }
        else if(strcmp(argv[i], "--n") == 0) { sscanf(argv[i+1],"%d", &numItems); check_parameters += 1; }
        else if(strcmp(argv[i], "--l") == 0) { sscanf(argv[i+1],"%d", &reloadingDepth); check_parameters += 1; }
        else if(strcmp(argv[i], "--h") == 0) { sscanf(argv[i+1],"%d", &costForEachRelocate); check_parameters += 1; }
        else if(strcmp(argv[i], "--outputsolution") == 0) { sscanf(argv[i+1],"%s", parameterStr); outputSolutionFileName = parameterStr; check_parameters += 1; }
        else check_parameters = -INF;        
    }
    
    if(check_parameters != 7) usage();
    
    Data::getInstance().readData(pickupAreaFileName, deliveryAreaFileName, numItems, reloadingDepth, costForEachRelocate);
      
    if(approachID == "ILP1") {
        ILPFormulation1 ILP1;
        ILP1.setAnInitialSolution();
        ILP1.solve(outputSolutionFileName);
    }
    else if(approachID == "ILP2") {
        ILPFormulation2 ILP2;
        ILP2.setAnInitialSolution();
        ILP2.solve(outputSolutionFileName);
    }   
    else if(approachID == "BRKGA") {
        runBRKGA(outputSolutionFileName);
    }
    else usage();
    
    return 0;
}
