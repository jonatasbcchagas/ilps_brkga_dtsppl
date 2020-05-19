#include "data.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <csignal>

void Data::readData(string pickupAreaFileName, string deliveryAreaFileName, int _numItems, int _reloadingDepth, int _costForEachRealoading) {

    numItems = _numItems;
    reloadingDepth = _reloadingDepth;
    costForEachRealoading = _costForEachRealoading;
    
    string line;
    line = pickupAreaFileName;
    for(int j=0;j<(int)line.length(); ++j) {
        if(line[j] < '0' || line[j] > '9') line[j] = ' ';
    }
    stringstream ss;
    ss << line;
    ss >> line;
    pickupAreaName = "R" + line + "p";
    
    line = deliveryAreaFileName;
    for(int j=0;j<(int)line.length(); ++j) {
        if(line[j] < '0' || line[j] > '9') line[j] = ' ';
    }
    ss.clear();
    ss << line;
    ss >> line;
    deliveryAreaName = "R" + line + "d";
    
    ifstream fin(pickupAreaFileName.c_str());

    if(!fin) {
        clog << "ERROR!" << endl;
        std::_Exit(EXIT_FAILURE);
    }

    for(int i = 0; i < 6; ++i) {
        getline(fin, line);
        if(line.find("DIMENSION") != string::npos) {
            for(int j=0;j<(int)line.length(); ++j) {
                if(line[j] < '0' || line[j] > '9') line[j] = ' ';
            }
            ss.clear();
            ss << line;
            ss >> numPoints;
        }
    }

    DEPOT = 0;

    pickupDistance.resize(numItems+1);
    for(int i = 0; i < numItems+1; ++i) {
        pickupDistance[i].resize(numItems+1);
    }

    deliveryDistance.resize(numItems+1);
    for(int i = 0; i < numItems+1; ++i) {
        deliveryDistance[i].resize(numItems+1);
    }

    vector < pair < long double, long double > > points;
    int id;
    long double x, y;
    
    while(fin >> id >> x >> y) {
        if((int)points.size() < numItems+1) {
            points.push_back(make_pair(x, y));
        }
    }

    for(int i = 0; i < (int)points.size(); ++i) {
        for(int j = i; j < (int)points.size(); ++j) {
            pickupDistance[i][j] = pickupDistance[j][i] = (int) (0.5 +
                                                                        sqrt(
                                                                            ((points[i].first-points[j].first)*(points[i].first-points[j].first)) +
                                                                            ((points[i].second-points[j].second)*(points[i].second-points[j].second))
                                                                        ));
        }
    }
    fin.close();

    fin.open(deliveryAreaFileName.c_str());

    if(!fin) {
        clog << "ERROR" << endl;
        exit(0);
    }

    for(int i = 0; i < 6; ++i) {
        getline(fin, line);
    }

    points.clear();

    while(fin >> id >> x >> y) {
        if((int)points.size() < numItems+1) {
            points.push_back(make_pair(x, y));
        }
    }

    for(int i = 0; i < (int)points.size(); ++i) {
        for(int j = i; j < (int)points.size(); ++j) {
            deliveryDistance[i][j] = deliveryDistance[j][i] = (int) (0.5 +
                                                                            sqrt(
                                                                                ((points[i].first-points[j].first)*(points[i].first-points[j].first)) +
                                                                                ((points[i].second-points[j].second)*(points[i].second-points[j].second))
                                                                            ));
        }
    }
    
    fin.close();
}
