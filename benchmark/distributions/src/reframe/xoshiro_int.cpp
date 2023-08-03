#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>

#ifdef FROM_SHELL_SCRIPT
#include "../../../extlibs/Xoshiro-cpp/XoshiroCpp.hpp"
#else
#include <XoshiroCpp.hpp>
#endif

namespace FS = std::filesystem;

using namespace XoshiroCpp;

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

    const std::uint64_t seed = 12345;
    Xoshiro256PlusPlus rg(seed);
    std::uniform_int_distribution<int> DIST(0,100);

    for (int i=0 ; i<1000 ; i++)
    {
        auto startXoshiroCpp = std::chrono::steady_clock::now();
        std::vector<int> xoshiroCppVector(numVectors);
        std::generate(xoshiroCppVector.begin(), xoshiroCppVector.end(), [&]() { return DIST(rg); });
        auto endXoshiroCpp = std::chrono::steady_clock::now();
        auto durationXoshiroCpp = std::chrono::duration_cast<std::chrono::microseconds>(endXoshiroCpp - startXoshiroCpp);
        mean += durationXoshiroCpp;
    }
    auto meanXoshiroCpp = mean.count() / 1000.0;
    std::cout << flag << ",xoshiro-cpp," << meanXoshiroCpp << std::endl;

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

    file << flag << ",xoshiro-cpp," << meanXoshiroCpp << std::endl;
    file.close();
    return 0;
}