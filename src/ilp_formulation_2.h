#ifndef ILP_FORMULATION_2_H
#define ILP_FORMULATION_2_H

#include "gurobi_c++.h"
#include "callbacks.h"

class ILPFormulation2 {
    
    private:

            GRBEnv *env = NULL;
            GRBModel *model = NULL;
            GRBVar ***chi = NULL;
            GRBVar **u = NULL;
            GRBVar ****y = NULL;
            GRBVar **z = NULL;     
            int status;
            
            void createVariables();
            void addObjectiveFunction();
            void addConstraints();
            void setParameters();
            
    public:
    
            ILPFormulation2();
            ~ILPFormulation2();
            
            void exportModel(string);
            void setAnInitialSolution();
            void solve(const string);
            int getTotalCost() const;
            int getTotalDistanceTraveled() const;
            int getTotalNumberOfRelocations() const;
            void saveSolution(const string);
};

#endif
