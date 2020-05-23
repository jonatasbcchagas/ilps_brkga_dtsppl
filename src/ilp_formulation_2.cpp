#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <cstring>

#include "data.h"
#include "ilp_formulation_2.h"
#include "tsp_solver.h"

using namespace std;

ILPFormulation2::ILPFormulation2() {
    
    createVariables();
    addObjectiveFunction();
    addConstraints();
    setParameters();
}

ILPFormulation2::~ILPFormulation2() {
    delete env;
}

void ILPFormulation2::createVariables() {
    
    env = new GRBEnv();
    model = new GRBModel(*env);

    char name[1000];
    
    // chi_{ijr}
    chi = new GRBVar**[Data::getInstance().numItems+1];
    for(int i = 0; i <= Data::getInstance().numItems; ++i) {
        chi[i] = new GRBVar*[Data::getInstance().numItems+1];
        for(int j = 0; j <= Data::getInstance().numItems; ++j) {
            chi[i][j] = new GRBVar[2]; // pickup and delivery
            for(int r = PICKUP; r <= DELIVERY; ++r) {
                sprintf(name, "chi_%02d_%02d_%c", i, j, (r == PICKUP ? 'P' : 'D'));
                chi[i][j][r] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
            }
        }
    }
    
    // u_{jr}
    u = new GRBVar*[Data::getInstance().numItems+1];
    for(int j = 0; j <= Data::getInstance().numItems; ++j) {
        u[j] = new GRBVar[2]; // pickup and delivery
        for(int r = PICKUP; r <= DELIVERY; ++r) {
            sprintf(name, "u_%02d_%c", j, (r == PICKUP ? 'P' : 'D'));
            u[j][r] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, name);
        }
    }

    // y_{rjkl}
    y = new GRBVar***[2]; // pickup and delivery
    for(int r = PICKUP; r <= DELIVERY; ++r) {
        y[r] = new GRBVar**[Data::getInstance().numItems+1];
        for(int j = 1; j <= Data::getInstance().numItems; ++j) {
            y[r][j] = new GRBVar*[Data::getInstance().numItems+2];
            for(int k = 1; k <= Data::getInstance().numItems+1; ++k) {
                y[r][j][k] = new GRBVar[r == PICKUP ? k + 1 : Data::getInstance().numItems - k + 1 + 1];
                for(int l = 1; l <= (r == PICKUP ? k : Data::getInstance().numItems - k + 1); ++l) {
                   sprintf(name, "y_%02d_%02d_%02d_%c", j, l, k, (r == PICKUP ? 'P' : 'D'));
                   y[r][j][k][l] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
               }
            }
        }
    }
    
    // z_{kr}
    z = new GRBVar*[Data::getInstance().numItems];
    for(int k = 1; k <= Data::getInstance().numItems - 1; ++k) {
        z[k] = new GRBVar[2]; // pickup and delivery
        for(int r = PICKUP; r <= DELIVERY; ++r) {
            sprintf(name, "z_%02d_%c", k, (r == PICKUP ? 'P' : 'D'));
            z[k][r] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, name);
        }
    }
}

void ILPFormulation2::addObjectiveFunction() {
    
    // Objective function
    
    GRBLinExpr _objPart1 = 0;
    GRBLinExpr _objPart2 = 0;
                
    for(int r = PICKUP; r <= DELIVERY; ++r) {
        for(int i = 0; i <= Data::getInstance().numItems; ++i) {
            for(int j = 0; j <= Data::getInstance().numItems; ++j) {
                if(j == i) continue;
                _objPart1 += chi[i][j][r] * (r == PICKUP ? Data::getInstance().pickupDistance[i][j] : Data::getInstance().deliveryDistance[i][j]);
            }
        }
    }
    for(int r = PICKUP; r <= DELIVERY; ++r) {
        for(int k = 1; k <= Data::getInstance().numItems - 1; ++k) {
            _objPart2 += z[k][r];
        }
    }
    
    GRBLinExpr obj = _objPart1 + (_objPart2 * Data::getInstance().costForEachRealoading);
    model->setObjective(obj, GRB_MINIMIZE);
    
    objPart1 = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER);
    model->addConstr(objPart1 == _objPart1);
    
    objPart2 = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER);
    model->addConstr(objPart2 == _objPart2); 
}

void ILPFormulation2::addConstraints() {
    
    // Constraints
    
    for(int r = PICKUP; r <= DELIVERY; ++r) {
        for(int i = 0; i <= Data::getInstance().numItems; ++i) {
            GRBLinExpr expr = 0;
            for(int j = 0; j <= Data::getInstance().numItems; ++j) {
                if(j == i) continue;
                expr += chi[i][j][r];
            }
            model->addConstr(expr == 1);
        }
    }
    
    for(int r = PICKUP; r <= DELIVERY; ++r) {
        for(int j = 0; j <= Data::getInstance().numItems; ++j) {
            GRBLinExpr expr = 0;
            for(int i = 0; i <= Data::getInstance().numItems; ++i) {
                if(j == i) continue;
                expr += chi[i][j][r];
            }
            model->addConstr(expr == 1);
        }
    }

    for(int r = PICKUP; r <= DELIVERY; ++r) {
        for(int i = 0; i <= Data::getInstance().numItems; ++i) {
            for(int j = 1; j <= Data::getInstance().numItems; ++j) {
                if(j == i) continue;
                model->addConstr(u[j][r] - (u[i][r] + 1 - Data::getInstance().numItems * (1 - chi[i][j][r])) >= 0);
            }
        }
    }

    for(int r = PICKUP; r <= DELIVERY; ++r) {
        model->addConstr(u[0][r] == 0);
    }

    for(int r = PICKUP; r <= DELIVERY; ++r) {
        for(int j = 0; j <= Data::getInstance().numItems; ++j) {
            model->addConstr(u[j][r] <= Data::getInstance().numItems);
        }
    }

    for(int k = 1; k <= Data::getInstance().numItems; ++k) {
        for(int l = 1; l <= k; ++l) {
            GRBLinExpr expr = 0;
            for(int j = 1; j <= Data::getInstance().numItems; ++j) {
                expr += y[PICKUP][j][k][l];
            }
            model->addConstr(expr == 1);
        }
    }
    
    for(int k = 1; k <= Data::getInstance().numItems; ++k) {
        for(int l = 1; l <= Data::getInstance().numItems - k + 1; ++l) {
            GRBLinExpr expr = 0;
            for(int j = 1; j <= Data::getInstance().numItems; ++j) {
                expr += y[DELIVERY][j][k][l];
            }
            model->addConstr(expr == 1);
        }
    }
        
    for(int k = 1; k <= Data::getInstance().numItems; ++k) { 
        for(int j = 1; j <= Data::getInstance().numItems; ++j) {
            GRBLinExpr expr = 0;
            for(int l = 1; l <= k; ++l) { 
                expr += y[PICKUP][j][k][l];
            }                        
            model->addConstr(expr <= 1);
            model->addConstr(expr * k - (k - u[j][PICKUP] + 1) >= 0);
        }
    }

    for(int k = 1; k <= Data::getInstance().numItems; ++k) {
        for(int j = 1; j <= Data::getInstance().numItems; ++j) {
            GRBLinExpr expr = 0;
            for(int l = 1; l <= Data::getInstance().numItems - k + 1; ++l) {
                expr += y[DELIVERY][j][k][l];
            }                        
            model->addConstr(expr <= 1);
            model->addConstr(expr * (Data::getInstance().numItems - k + 1) - (u[j][DELIVERY] - k + 1) >= 0);
        }
    }    
      
    for(int l = 1; l <= Data::getInstance().numItems; ++l) {
        for(int j = 1; j <= Data::getInstance().numItems; ++j) {
            model->addConstr(y[DELIVERY][j][1][l] - y[PICKUP][j][Data::getInstance().numItems][l] == 0);
        }
    }
    
    for(int k = 1; k <= Data::getInstance().numItems-1; ++k) {
        for(int l = 1; l <= k - Data::getInstance().reloadingDepth; ++l) {
            for(int j = 1; j <= Data::getInstance().numItems; ++j) {
                model->addConstr(y[PICKUP][j][k][l] - y[PICKUP][j][Data::getInstance().numItems][l] == 0);
            }
        }
    }
    
    for(int k = 2; k <= Data::getInstance().numItems; ++k) {
        for(int l = 1; l <= Data::getInstance().numItems - k - Data::getInstance().reloadingDepth + 1; ++l) {
            for(int j = 1; j <= Data::getInstance().numItems; ++j) {
                model->addConstr(y[DELIVERY][j][k][l] - y[DELIVERY][j][1][l] == 0);
            }
        }
    }
  
    for(int k = 1; k <= Data::getInstance().numItems - 1; ++k) {
        for(int l = 1; l <= k; ++l) {
            for(int j = 1; j <= Data::getInstance().numItems; ++j) {
                model->addConstr(z[k][PICKUP] >= (y[PICKUP][j][k][l] - y[PICKUP][j][k+1][l]) * (k - l + 1));
            }
        }
    }

    for(int k = 1; k <= Data::getInstance().numItems - 1; ++k) {
        for(int l = 1; l <= Data::getInstance().numItems - k; ++l) {
            for(int j = 1; j <= Data::getInstance().numItems; ++j) {
                model->addConstr(z[k][DELIVERY] >= (y[DELIVERY][j][k+1][l] - y[DELIVERY][j][k][l]) * (Data::getInstance().numItems - k - l + 1));
            }
        }
    }
}

void ILPFormulation2::setParameters() {

    model->set(GRB_DoubleParam_TimeLimit, 3600.0);
    model->set(GRB_IntParam_Threads, 1);
    model->set(GRB_IntParam_LogToConsole, 0);
}

void ILPFormulation2::exportModel(string fileLP) {
    model->write(fileLP.c_str());
}

void ILPFormulation2::setAnInitialSolution() {

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

    vector < vector < bool > > adjMatrix((int)result.second.size()-1, vector < bool > ((int)result.second.size()-1, false));

    for(int i = 0; i < (int)result.second.size()-1; ++i) {
        adjMatrix[result.second[i]][result.second[i+1]] = true;
    }

    for(int i = 0; i <= Data::getInstance().numItems; ++i) {
        for(int j = 0; j <= Data::getInstance().numItems; ++j) {
            if(i == j) continue;
            if(adjMatrix[i][j]) {
                chi[i][j][PICKUP].set(GRB_DoubleAttr_Start, 1.0);
                chi[j][i][DELIVERY].set(GRB_DoubleAttr_Start, 1.0);
            }
            else {
                chi[i][j][PICKUP].set(GRB_DoubleAttr_Start, 0.0);
                chi[j][i][DELIVERY].set(GRB_DoubleAttr_Start, 0.0);
            }
        }
    }

    for(int i = 0; i < (int)result.second.size()-1; ++i) {
        u[result.second[i]][PICKUP].set(GRB_DoubleAttr_Start, i);
        u[result.second[(int)result.second.size()-1-i]][DELIVERY].set(GRB_DoubleAttr_Start, i);
    }

    for(int k = 1; k <= Data::getInstance().numItems; ++k) {
        for(int l = 1; l <= k; ++l) {
            for(int j = 1; j <= Data::getInstance().numItems; ++j) {
                if(j == result.second[l]) {
                    y[PICKUP][j][k][l].set(GRB_DoubleAttr_Start, 1.0);
                }
                else {
                    y[PICKUP][j][k][l].set(GRB_DoubleAttr_Start, 0.0);
                }
            }
        }
    }

    for(int k = 1; k <= Data::getInstance().numItems; ++k) {
        for(int l = 1; l <= Data::getInstance().numItems - k + 1; ++l) {
            for(int j = 1; j <= Data::getInstance().numItems; ++j) {
                if(j == result.second[l]) {
                    y[DELIVERY][j][k][l].set(GRB_DoubleAttr_Start, 1.0);
                }
                else {
                    y[DELIVERY][j][k][l].set(GRB_DoubleAttr_Start, 0.0);
                }
            }
        }
    }

    for(int k = 1; k <= Data::getInstance().numItems - 1; ++k) {
        z[k][PICKUP].set(GRB_DoubleAttr_Start, 0.0);
        z[k][DELIVERY].set(GRB_DoubleAttr_Start, 0.0);
    }
}

void ILPFormulation2::solve(const string outputSolutionFileName) {

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    model->set(GRB_StringParam_LogFile, outputSolutionFileName + ".gurobilog");
        
    LogCallback cb(objPart1, objPart2);
    
    model->setCallback(&cb);
    
    model->optimize();
    
    cb.saveSummarizedLog(outputSolutionFileName + ".log");

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration < double > time_span = duration_cast < duration < double > > (t2 - t1);

    int status = model->get(GRB_IntAttr_Status);

    if(status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED) {
        clog << "The model cannot be solved because it is infeasible or unbounded" << endl;
        exit(0);
    }

    if(model->get(GRB_IntAttr_SolCount) >= 1) {
        saveSolution(outputSolutionFileName + ".sol");
    }

    double currentLowerBound = (double)model->get(GRB_DoubleAttr_ObjBoundC);
    int lowerBound = (int)(model->get(GRB_DoubleAttr_ObjBound)+0.5);
    int upperBound = (int)(model->get(GRB_DoubleAttr_ObjVal)+0.5);
    int totalDistanceTraveled = getTotalDistanceTraveled();
    int totalNumberOfRelocations = getTotalNumberOfRelocations();
    double gap = (double)model->get(GRB_DoubleAttr_MIPGap);

    char tmp[10000];
    strcpy(tmp, outputSolutionFileName.c_str());
    strcat(tmp, ".log");    
    ofstream fout(tmp, ofstream::app);  
    if(cb.lastLB != lowerBound) {
        sprintf(tmp, "%15d %15d %17.1lf", upperBound, lowerBound, time_span.count());    
        fout << tmp << endl;
    }
    fout << endl;   
    sprintf(tmp, "%-20s %15.3lf %15d %15d %15d %15d %15.5lf %15.3lf %15d\n", outputSolutionFileName.c_str(), currentLowerBound, lowerBound, upperBound, totalDistanceTraveled, totalNumberOfRelocations, gap, time_span.count(), status == GRB_OPTIMAL ? 1 : 0);
    fout << tmp << endl;
    
    cb.nds.saveSet(outputSolutionFileName + ".nds");
}

int ILPFormulation2::getTotalCost() const {
    
    if(model->get(GRB_IntAttr_SolCount) == 0) {
        return INF;
    }
    
    if(((int)(model->get(GRB_DoubleAttr_ObjVal)+0.5)) != getTotalDistanceTraveled() + Data::getInstance().costForEachRealoading * getTotalNumberOfRelocations()) {
        clog << "Error getTotalCost()!" << endl;
        clog << (int)(model->get(GRB_DoubleAttr_ObjVal)+0.5) << ' ' << getTotalDistanceTraveled() + Data::getInstance().costForEachRealoading * getTotalNumberOfRelocations() << endl;
        return INF;
    }
        
    return getTotalDistanceTraveled() + Data::getInstance().costForEachRealoading * getTotalNumberOfRelocations();
}

int ILPFormulation2::getTotalDistanceTraveled() const {
    
    if(model->get(GRB_IntAttr_SolCount) == 0) return INF;
    return (int)(objPart1.get(GRB_DoubleAttr_X)+0.5);
}

int ILPFormulation2::getTotalNumberOfRelocations() const {
    
    if(model->get(GRB_IntAttr_SolCount) == 0) return INF;
    return (int)(objPart2.get(GRB_DoubleAttr_X)+0.5);
}

void ILPFormulation2::saveSolution(const string outputSolutionFileName) {
    
    int totalCost = getTotalCost();
    int totalDistanceTraveled = getTotalDistanceTraveled();
    int totalNumberOfRelocations = getTotalNumberOfRelocations();    

    ofstream fout(outputSolutionFileName.c_str());

    fout << "Total cost: " << totalCost << endl;
    fout << "Distance traveled: " << totalDistanceTraveled << endl;
    fout << "Number of relocations: " << totalNumberOfRelocations << endl << endl;

    vector< int > pickupTour, deliveryTour;
    
    pickupTour.push_back(0);
    deliveryTour.push_back(0);    

    int prev = 0;
    while(1) {
        for(int j = 0; j <= Data::getInstance().numItems; ++j) {
            if(j != prev && chi[prev][j][PICKUP].get(GRB_DoubleAttr_X) >= 0.5) {
                pickupTour.push_back(j);
                prev = j;
                break;
            }
        }
        if(prev == 0) break;
    }

    prev = 0;
    while(1) {
        for(int j = 0; j <= Data::getInstance().numItems; ++j) {
            if(j != prev && chi[prev][j][DELIVERY].get(GRB_DoubleAttr_X) >= 0.5) {
                deliveryTour.push_back(j);
                prev = j;
                break;
            }
        }
        if(prev == 0) break;
    }

    vector < vector < int > > container;
    container.resize(2 * Data::getInstance().numItems + 1);
    for(int i = 0; i < (int)container.size(); ++i) {
        container[i].resize(Data::getInstance().numItems+1);
        for(int j = 0; j < (int)container[i].size(); ++j) {
            container[i][j] = -1;
        }
    }
    
    for(int k = 1; k <= Data::getInstance().numItems; ++k) {
        for(int l = 1; l <= k; ++l) {
            for(int j = 1; j <= Data::getInstance().numItems; ++j) {
                if(y[PICKUP][j][k][l].get(GRB_DoubleAttr_X) > 0.5) {
                    container[k][l] = j;
                }
            }
        }
    }
    
    for(int k = 1; k <= Data::getInstance().numItems; ++k) {
        for(int l = 1; l <= Data::getInstance().numItems - k + 1; ++l) {
            for(int j = 1; j <= Data::getInstance().numItems; ++j) {
                if(y[DELIVERY][j][k][l].get(GRB_DoubleAttr_X) > 0.5) {
                    container[k + Data::getInstance().numItems][l] = j;
                }
            }
        }
    }
    
    fout << "Loading/unloading plan timeline:" << endl << endl;

    for(int j = Data::getInstance().numItems; j >= 1; --j) {
        for(int i = 1; i < (int)container.size(); ++i) {
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
}
