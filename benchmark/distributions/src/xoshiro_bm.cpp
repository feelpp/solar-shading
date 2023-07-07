#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>


#include <XoshiroCpp.hpp>

using namespace XoshiroCpp;

int main(int argc, char *argv[]){
    std::string flag;
    std::stringstream ss;
    
    if (argc > 1) {
        // Start from the second argument (index 1) since the first argument (index 0) is the program name
        for(int i = 1; i < argc; ++i) {
            // If it's not the first argument, add a space before
            if(i != 1) {
                ss << " ";
            }
            ss << argv[i];
        }

        flag = ss.str();
    }
    else {
        flag = "default";
    }
    
    
    const int numVectors = 100000000;
    std::vector<int> xoshiroCppVector(numVectors);


    const std::uint64_t seed = 12345;

    Xoshiro256PlusPlus rg(seed);

    std::uniform_int_distribution<int> DIST(0,100);

    auto startXoshiroCpp = std::chrono::steady_clock::now();

    for (int i = 0; i < numVectors; i++) {
        xoshiroCppVector[i] =  DIST(rg);
    }

    auto endXoshiroCpp = std::chrono::steady_clock::now();
    auto durationXoshiroCpp = std::chrono::duration_cast<std::chrono::milliseconds>(endXoshiroCpp - startXoshiroCpp);

    std::cout << flag << ",XoshiroCpp," << durationXoshiroCpp.count() << std::endl;
}