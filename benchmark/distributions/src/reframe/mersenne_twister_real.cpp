#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <filesystem>

namespace FS = std::filesystem;

int main(int argc, char *argv[]){
    std::string flag;
    std::stringstream ss;
    std::chrono::microseconds mean = std::chrono::microseconds::zero();
    
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
    
    std::random_device rd;
    std::mt19937 mersenne(rd());
    std::uniform_real_distribution<double> mersenneDist(0, 100);

    for(int i = 0; i < 1000; ++i)
    {   
        auto startMersenne = std::chrono::steady_clock::now();
        std::vector<double> mersenneVector(numVectors);
        std::generate(mersenneVector.begin(), mersenneVector.end(), [&]() { return mersenneDist(mersenne); });
        auto endMersenne = std::chrono::steady_clock::now();
        auto durationMersenne = std::chrono::duration_cast<std::chrono::microseconds>(endMersenne - startMersenne);
        mean += durationMersenne;
    }
    auto meanMicros = mean.count() / 1000.0;
    std::cout << flag << ",MersenneTwister," << meanMicros << std::endl;

    // write the output to a csv file located at the root of the github repos 
    
    std::ofstream file;
    FS::path scriptPath = FS::current_path();
    FS::path scriptDirectory = scriptPath.parent_path();
    FS::path rootDirectory1 = scriptDirectory.parent_path();
    FS::path rootDirectory = rootDirectory1.parent_path();
    FS::path filePath = rootDirectory / "solar-shading" / "results" / "csv" / "uniform_real.csv";

    if (FS::exists(filePath)) {
        file.open(filePath, std::ios_base::app);
    }
    else {
        file.open(filePath);
        file << "Flag,Method,Duration" << std::endl;
    }

    file << flag << ",MersenneTwister," << meanMicros << std::endl;
    file.close();

    return 0;
}
