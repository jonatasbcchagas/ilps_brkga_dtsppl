#ifndef ILP_FORMULATION_1_H
#define ILP_FORMULATION_1_H

#include "gurobi_c++.h"
#include "callbacks.h"

class ILPFormulation1 {
    
    private:

            GRBEnv *env = NULL;
            GRBModel *model = NULL;
            GRBVar ****x = NULL;
            GRBVar ****y = NULL;
            GRBVar **z = NULL;  
            GRBVar objPart1;
            GRBVar objPart2;
            int status;
            
            void createVariables();
            void addObjectiveFunction();
            void addConstraints();
            void setParameters();
            
    public:
    
            ILPFormulation1();
            ~ILPFormulation1();
            
            void exportModel(string);
            void setAnInitialSolution();
            void solve(const string);
            int getTotalCost() const;
            int getTotalDistanceTraveled() const;
            int getTotalNumberOfRelocations() const;
            void saveSolution(const string);
};

#endif
