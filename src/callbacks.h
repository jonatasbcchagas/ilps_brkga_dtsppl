#ifndef LOG_CALLBACK_H
#define LOG_CALLBACK_H

#include "gurobi_c++.h"
#include <string>

#include "non_dominated_set.h"

class LogCallback: public GRBCallback {

    public:

            NonDominatedSet nds;
    
            GRBVar f1, f2;
            int lastLB, lastUB;
            vector < string > logLines;
            char logLine[10000];
        
            LogCallback(GRBVar &_f1, GRBVar &_f2) {
                f1 = _f1; f2 = _f2;
                lastLB = lastUB = 0;
                logLines.push_back("             UB              LB           Time(s)");
            }
            
            void saveSummarizedLog(string fileName) {
                ofstream fout(fileName.c_str());
                for(int i = 0; i < (int)logLines.size(); ++i) fout << logLines[i] << endl; 
                fout.close();
            }
    
    protected:
        
            void callback () {
                
                try {                    
                    if (where == GRB_CB_MIP) {
                        double objbst = getDoubleInfo(GRB_CB_MIP_OBJBST);
                        double objbnd = getDoubleInfo(GRB_CB_MIP_OBJBND);
                        double runtime = getDoubleInfo(GRB_CB_RUNTIME);
                        bool updated = false;
                        if(int(objbnd + 0.5) != lastLB) { lastLB = int(objbnd + 0.5); updated = true; }
                        if(lastLB != 0 && int(objbst + 0.5) != lastUB) { lastUB = int(objbst + 0.5); updated = true; }
                        if(updated == true) { sprintf(logLine, "%15d %15d %17.1lf", int(objbst + 0.5), int(objbnd + 0.5), runtime); logLines.push_back(logLine); }
                    }
                    if (where == GRB_CB_MIPSOL) {
                        nds.add(make_pair(getSolution(f1), getSolution(f2)));
                    }
                } 
                catch (GRBException e) {
                    clog << "Error number: " << e.getErrorCode() << endl;
                    clog << e.getMessage() << endl;
                } catch (...) {
                    clog << "Error during callback" << endl;
                }
            }
};

#endif
