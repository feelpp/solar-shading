#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>

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

    std::vector<int> stdVector(numVectors);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 100);

    auto startStdUniform = std::chrono::steady_clock::now();

    for (int i = 0; i < numVectors; ++i) {
        stdVector[i] = dist(gen);
    }

    auto endStdUniform = std::chrono::steady_clock::now();
    auto durationStdUniform = std::chrono::duration_cast<std::chrono::milliseconds>(endStdUniform - startStdUniform);

    std::cout << flag << ",std::uniform_int_distribution," << durationStdUniform.count() << std::endl;
}