#ifndef NON_DOMINATED_SET_H
#define NON_DOMINATED_SET_H

#include <iostream>
#include <fstream>
#include <list>

using namespace std;

class NonDominatedSet {

    private:
        
            list < pair < int, int > > allNDS;

    public:
        
            NonDominatedSet() {};
            
            int getRelation(pair < int, int > a, pair < int, int > b) {
                
                int val = 0;
            
                if (a.first < b.first) {
                    if (val == -1) return 0;
                    val = 1;
                } else if (a.first > b.first) {
                    if (val == 1) return 0;
                    val = -1;
                }
                
                if (a.second < b.second) {
                    if (val == -1) return 0;
                    val = 1;
                } else if (a.second > b.second) {
                    if (val == 1) return 0;
                    val = -1;
                }

                return val;
            }
            
            void add(pair < int, int > s) {

                bool isAdded = true;

                list < pair < int, int > > :: iterator it = allNDS.begin();
                
                for(; it != allNDS.end(); ++it) {
                    
                    pair < int, int > other = *it;

                    int rel = getRelation(s, other);

                    // if dominated by or equal in design space
                    if (rel == -1 || (rel == 0 && s == other)) {
                        isAdded = false;
                        break;
                    } else if (rel == 1) it = allNDS.erase(it);

                }

                if (isAdded) allNDS.push_back(s);
            }

            void saveSet(const string outputFileName) {
                
                ofstream fout(outputFileName.c_str());
                allNDS.sort();
                fout << fixed << setw(10) << "F1" << ' ' << fixed << setw(10) << "F2" << endl;
                list < pair < int, int > > :: iterator it = allNDS.begin();                
                for(; it != allNDS.end(); ++it) {
                    fout << fixed << setw(10) << it->first << ' ' << setw(10) << it->second << endl;
                }
                fout.close();
            }
};

#endif

