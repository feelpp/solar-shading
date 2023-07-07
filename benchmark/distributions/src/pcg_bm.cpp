#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>

#include <pcg_random.hpp>


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
    std::vector<int> pcgCppVector(numVectors);

    pcg32_fast rng;  // Create an instance of the pcg32 random number generator
    rng.seed(42);  // Seed the pcg32 random number generator
    std::uniform_int_distribution<int> pcgDist(0, 100);  // Create a uniform distribution over the range [0, 100

    auto startPcgCpp = std::chrono::steady_clock::now();

    for (int i = 0; i < numVectors; ++i) {
        // pcgCppVector[i] = static_cast<int>(rng(101));  // Generate random integers between 0 and 100
        pcgCppVector[i] = pcgDist(rng);
    }

    auto endPcgCpp = std::chrono::steady_clock::now();
    auto durationPcgCpp = std::chrono::duration_cast<std::chrono::milliseconds>(endPcgCpp - startPcgCpp);

    std::cout << flag << ",pcg-cpp," << durationPcgCpp.count() << std::endl;
}