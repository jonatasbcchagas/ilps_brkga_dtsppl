#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include "gurobi_c++.h"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>

using namespace std;

// Subtour elimination callback.  Whenever a feasible solution is found,
// find the smallest subtour, and add a subtour elimination constraint
// if the tour doesn't visit every node.

class Subtourelim: public GRBCallback {
    
    public:
        
        GRBVar** vars;
        int n;
        Subtourelim(GRBVar** xvars, int xn) {
            vars = xvars;
            n = xn;
        }
        
        void findsubtour(int n, double** sol, int* tourlenP, int* tour) {

            bool* seen = new bool[n];
            int bestind, bestlen;
            int i, node, len, start;

            for (i = 0; i < n; i++) {
                seen[i] = false;
            }

            start = 0;
            bestlen = n+1;
            bestind = -1;
            node = 0;
            
            while (start < n) {
                
                for (node = 0; node < n; node++) {
                    if (!seen[node]) break;
                }
                if (node == n) break;
                
                for (len = 0; len < n; len++) {
                    tour[start+len] = node;
                    seen[node] = true;
                    for (i = 0; i < n; i++) {
                        if (sol[node][i] > 0.5 && !seen[i]) {
                            node = i;
                            break;
                        }
                    }
                    if (i == n) {
                        len++;
                        if (len < bestlen) {
                            bestlen = len;
                            bestind = start;
                        }
                        start += len;
                        break;
                    }
                }
            }

            for (i = 0; i < bestlen; i++) {
                tour[i] = tour[bestind+i];
            }
            *tourlenP = bestlen;

            delete[] seen;
        }

    protected:
        
        void callback() {
            
            try {
                
                if (where == GRB_CB_MIPSOL) {
                    
                    // Found an integer feasible solution - does it visit every node?
                    double **x = new double*[n];
                    int *tour = new int[n];
                    int i, j, len;
                    
                    for (i = 0; i < n; i++) {
                        x[i] = getSolution(vars[i], n);
                    }

                    findsubtour(n, x, &len, tour);

                    if (len < n) {
                        // Add subtour elimination constraint
                        GRBLinExpr expr = 0;
                        for (i = 0; i < len; i++) {
                            for (j = i+1; j < len; j++) {
                                expr += vars[tour[i]][tour[j]];
                            }
                            addLazy(expr <= len-1);
                        }

                        for (i = 0; i < n; i++) {
                            delete[] x[i];
                        }
                        
                        delete[] x;
                        delete[] tour;
                    }
                }
            } catch (GRBException e) {
                clog << "Error number: " << e.getErrorCode() << endl;
                clog << e.getMessage() << endl;
            } catch (...) {
                clog << "Error during callback" << endl;
            }
    }
};

class TSPSolver {

    public:
            TSPSolver(){}

            pair < int, vector < int > > solve(int n, const vector < vector < int > > &distance) {
                
                vector < int > bestTour;
                int bestTourCost;

                GRBEnv *env = NULL;
                GRBVar **vars = NULL;

                int i, j;
                
                vars = new GRBVar*[n];
                for (i = 0; i < n; i++) {
                    vars[i] = new GRBVar[n];
                }

                try {
                    
                    env = new GRBEnv();
                    GRBModel model = GRBModel(*env);

                    model.getEnv().set(GRB_IntParam_OutputFlag, 0);
                    
                    // Must set LazyConstraints parameter when using lazy constraints
                    model.set(GRB_IntParam_LazyConstraints, 1);

                    // Create binary decision variables
                    for (i = 0; i < n; i++) {
                        for (j = 0; j <= i; j++) {
                            vars[i][j] = model.addVar(0.0, 1.0, distance[i][j], GRB_BINARY);
                            vars[j][i] = vars[i][j];
                        }
                    }

                    // Degree-2 constraints
                    for (i = 0; i < n; i++) {
                        GRBLinExpr expr = 0;
                        for (j = 0; j < n; j++) {
                            expr += vars[i][j];
                        }
                        model.addConstr(expr == 2);
                    }

                    // Forbid edge from node back to itself
                    for (i = 0; i < n; i++) {
                        vars[i][i].set(GRB_DoubleAttr_UB, 0);
                    }

                    // Set callback function
                    Subtourelim cb = Subtourelim(vars, n);
                    model.setCallback(&cb);

                    // Optimize model
                    model.optimize();

                    // Extract solution

                    if (model.get(GRB_IntAttr_SolCount) > 0) {
                        
                        double **sol = new double*[n];
                        
                        for (i = 0; i < n; i++) {
                            sol[i] = model.get(GRB_DoubleAttr_X, vars[i], n);
                        }

                        int* tour = new int[n];
                        int len;

                        cb.findsubtour(n, sol, &len, tour);
                        
                        
                        if(tour[1] > tour[len-1]) reverse(tour + 1, tour + len);
                            
                        assert(len == n);

                        bestTourCost = 0;
                        for (i = 0; i < len; i++) {
                            bestTour.push_back(tour[i]);
                            bestTourCost += distance[tour[i]][tour[(i+1)%len]];
                        }
                        bestTour.push_back(tour[0]);
                            
                        for (i = 0; i < n; i++) {
                            delete[] sol[i];
                        }
                        delete[] sol;
                        delete[] tour;
                    }

                } catch (GRBException e) {
                    clog << "Error number: " << e.getErrorCode() << endl;
                    clog << e.getMessage() << endl;
                } catch (...) {
                    clog << "Error during optimization" << endl;
                }

                for (i = 0; i < n; i++) {
                    delete[] vars[i];
                }
                delete[] vars;
                delete env;
                
                return make_pair(bestTourCost, bestTour);
            }
            
};

#endif
