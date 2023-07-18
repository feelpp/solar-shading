#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>

namespace FS = std::filesystem;

#ifdef FROM_SHELL_SCRIPT
#include "../../../extlibs/pcg-cpp/include/pcg_random.hpp"

#else
#include <pcg_random.hpp>
#endif

int main(int argc, char *argv[]){
    std::string flag;
    std::stringstream ss;
    std::chrono::microseconds mean= std::chrono::microseconds::zero();
    
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
    
    
    const int numVectors = 10000;

    pcg32_fast rng; 
    rng.seed(42); 
    std::uniform_int_distribution<int> pcgDist(0, 100);
    for (int i = 0 ; i < 1000 ; i++)
    {
        auto startPcgCpp = std::chrono::steady_clock::now();
        std::vector<int> pcgCppVector(numVectors);
        std::generate(pcgCppVector.begin(), pcgCppVector.end(), [&]() { return pcgDist(rng); });
        auto endPcgCpp = std::chrono::steady_clock::now();
        auto durationPcgCpp = std::chrono::duration_cast<std::chrono::microseconds>(endPcgCpp - startPcgCpp);

        mean += durationPcgCpp;
    }
    auto meanPcgCpp = mean.count() / 1000.0;
    std::cout << flag << ",pcg-cpp," << meanPcgCpp << std::endl;

    // write the output to a csv file located at the root of the github repos 
    
    std::ofstream file;
    FS::path scriptPath = FS::current_path();
    FS::path scriptDirectory = scriptPath.parent_path();
    FS::path rootDirectory1 = scriptDirectory.parent_path();
    FS::path rootDirectory = rootDirectory1.parent_path();
    FS::path filePath = rootDirectory / "solar-shading" / "results" / "csv" / "uniform_int.csv";

    if (FS::exists(filePath)) {
        file.open(filePath, std::ios_base::app);
    }
    else {
        file.open(filePath);
        file << "Flag,Method,Duration" << std::endl;
    }

    file << flag << ",pcg-cpp," << meanPcgCpp << std::endl;
    file.close();

    return 0;
}