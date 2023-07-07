#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>

#include <EigenRand/EigenRand>

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
    Eigen::VectorXi EigenVector(numVectors);
    
    Eigen::Rand::Vmt19937_64 urng{ 42 };

    auto startEigenRand = std::chrono::steady_clock::now();

    EigenVector = Eigen::Rand::uniformInt<Eigen::VectorXi>(numVectors, 1, urng, 0, 100);

    auto endEigenRand = std::chrono::steady_clock::now();
    auto durationEigenRand = std::chrono::duration_cast<std::chrono::milliseconds>(endEigenRand - startEigenRand);

    std::cout << flag << ",EigenRand," << durationEigenRand.count() << std::endl;
}