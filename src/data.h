#ifndef DATA_H
#define DATA_H

#include <vector>
#include <string>

using namespace std;

const int PICKUP = 0;
const int DELIVERY = 1;
const int INF = 987654321;
const double EPS = 1e-5;

class Data {
    
    public:
    
            static Data& getInstance() {
                static Data instance;                                        
                return instance;
            }
            
            Data(Data const&)            = delete;
            void operator=(Data const&)  = delete;
            
            string pickupAreaName, deliveryAreaName;
            int numPoints;
            int numItems;
            int capacityOfFleet;
            std::vector < std::vector < int > > pickupDistance, deliveryDistance;
            int DEPOT;
            int reloadingDepth;
            int costForEachRealoading;
            
            void readData(std::string, std::string, int, int = 2, int = 10);
            
    private:

            Data() {};
};

#endif
