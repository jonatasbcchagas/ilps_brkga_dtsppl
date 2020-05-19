#include <iomanip>
#include <fstream>

#include "brkga.h"

Decoder::Decoder(double _alpha, double _beta) {
    alpha = _alpha;
    beta = _beta;
}
        
Decoder::~Decoder() {}

double Decoder::decode(const std::vector< double >& chromosome, string solutionFileOut) {

    int n = Data::getInstance().numItems;
    int l = Data::getInstance().reloadingDepth;

    vector < vector < pair < double, int > > > permutations;

    int id = 0;
    
    permutations.push_back(vector < pair < double, int > > ());
    for(int i = 0; i < n; ++i) {
        permutations.back().push_back(make_pair(chromosome[id], i+1));
        id += 1;
    }
    
    for(int k = 1; k <= n; ++k) {
        permutations.push_back(vector < pair < double, int > > ());
        for(int i = 0; i < min(k, l + 1); ++i) {
            permutations.back().push_back(make_pair(chromosome[id], i));
            id += 1;
        }
    }
    
    for(int k = 1; k <= n; ++k) {
        permutations.push_back(vector < pair < double, int > > ());
        for(int i = 0; i < min(n - k + 1, l + 1); ++i) {
            permutations.back().push_back(make_pair(chromosome[id], i));
            id += 1;
        }
    }
    
    for(int i = 0; i < (int)permutations.size(); ++i) {
        sort(permutations[i].begin(), permutations[i].end());
    }
    
    vector < int > pickupTour;
    vector < int > deliveryTour;
    vector < int > stack;
    int numberOfRelocations = 0;
    int distance = 0;
    
    pickupTour.push_back(0);
    for(int i = 0; i < n; ++i) {
        pickupTour.push_back(permutations[0][i].second);
    }
    pickupTour.push_back(0);
    
    for(int k = 1; k <= n; ++k) {
        
        vector < pair < int, int > > vt;
        
        vt.push_back(make_pair(permutations[k][0].second, pickupTour[k]));            
        for(int i = 1; i < (int) permutations[k].size(); ++i) {
            vt.push_back(make_pair(permutations[k][i].second, stack.back()));
            stack.erase(--stack.end());
        }
        
        for(int i = 0; i < (int)vt.size(); ++i) {
            if(vt[i].first != i) {
                numberOfRelocations += (int)vt.size() - i - 1;
                break;
            }
        }
        
        sort(vt.rbegin(), vt.rend());
        for(int i = 0; i < (int)vt.size(); ++i) {
            stack.push_back(vt[i].second);
        }
    }
    
    deliveryTour.push_back(0);
    
    for(int k = n+1; k <= n+n; ++k) {
    
        vector < pair < int, int > > vt;
        
        for(int i = 0; i < (int) permutations[k].size(); ++i) {
            vt.push_back(make_pair(permutations[k][i].second, stack.back()));
            stack.erase(--stack.end());
        }
        
        for(int i = 0; i < (int)vt.size(); ++i) {
            if(vt[i].first != i) {
                numberOfRelocations += (int)vt.size() - i - 1;
                break;
            }
        }
        
        sort(vt.rbegin(), vt.rend());
        for(int i = 0; i < (int)vt.size(); ++i) {
            stack.push_back(vt[i].second);            
        }
        
        deliveryTour.push_back(stack.back());
        
        stack.erase(--stack.end());
    }
    deliveryTour.push_back(0);
    
    distance = 0;
    for(int i = 1; i < (int)pickupTour.size(); ++i) {
        distance += Data::getInstance().pickupDistance[pickupTour[i-1]][pickupTour[i]];
    }
    for(int i = 1; i < (int)deliveryTour.size(); ++i) {
        distance += Data::getInstance().deliveryDistance[deliveryTour[i-1]][deliveryTour[i]];
    }
        
    int totalCost = alpha * distance + beta * Data::getInstance().costForEachRealoading * numberOfRelocations;
    
    if(solutionFileOut != "") {
            
        vector < vector < int > > container;
        stack.clear();
    
        for(int k = 1; k <= n; ++k) {
            
            vector < pair < int, int > > vt;
            
            vt.push_back(make_pair(permutations[k][0].second, pickupTour[k]));            
            for(int i = 1; i < (int) permutations[k].size(); ++i) {
                vt.push_back(make_pair(permutations[k][i].second, stack.back()));
                stack.erase(--stack.end());
            }
            
            sort(vt.rbegin(), vt.rend());
            for(int i = 0; i < (int)vt.size(); ++i) {
                stack.push_back(vt[i].second);
            }
            container.push_back(stack);
        }
        
        container.push_back(stack);

        for(int k = n+1; k <= n+n; ++k) {
            
            vector < pair < int, int > > vt;
            
            for(int i = 0; i < (int) permutations[k].size(); ++i) {
                vt.push_back(make_pair(permutations[k][i].second, stack.back()));
                stack.erase(--stack.end());
            }
            
            sort(vt.rbegin(), vt.rend());
            for(int i = 0; i < (int)vt.size(); ++i) {
                stack.push_back(vt[i].second);            
            }
            
            stack.erase(--stack.end());
            
            container.push_back(stack);
        }
        
        container.erase(--container.end());
        
        for(int i = 0; i < (int)container.size(); ++i) {
            while((int)container[i].size() < n) container[i].push_back(-1);
        }
        
        ofstream fout(solutionFileOut.c_str());
    
        fout << "Total cost: " << totalCost << endl;
        fout << "Distance traveled: " << distance << endl;
        fout << "Number of relocations: " << numberOfRelocations << endl << endl;
        
        fout << "Loading/unloading plan timeline:" << endl << endl;

        for(int j = n-1; j >= 0; --j) {
            for(int i = 0; i < (int)container.size(); ++i) {
                if(container[i][j] == -1) fout << "   ";
                else fout << setfill('0') << setw(2) << container[i][j] << ' ';
            }
            fout << endl;
        }
        
        fout << endl;

        fout << "Pickup tour  : 00";
        for(int i = 1; i < (int)pickupTour.size(); ++i) {
            fout << " --> " << setfill('0') << setw(2) << pickupTour[i];
        }
        fout << endl;
        fout << "Delivery tour: 00";
        for(int i = 1; i < (int)deliveryTour.size(); ++i) {
            fout << " --> " << setfill('0') << setw(2) << deliveryTour[i];
        }
        fout << endl;
        
        fout.close();
    }
       
    return totalCost;
}
