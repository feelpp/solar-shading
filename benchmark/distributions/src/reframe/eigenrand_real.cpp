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
#include "../../../extlibs/EigenRand/EigenRand/EigenRand"
#else
#include "EigenRand/EigenRand"
#endif

using namespace Eigen;

int main(int argc, char *argv[]){
    std::chrono::microseconds meanEigenRand = std::chrono::microseconds::zero();
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
    const int numVectors = 10000;

    Rand::P8_mt19937_64 urng{ 42 };
    Rand::UniformRealGen<double> unif_gen{ 0, 100 };
    for (int i = 0; i < 1000; i++)
    {
        auto startEigenRand = std::chrono::steady_clock::now();
        VectorXd EigenVector = unif_gen.generate<VectorXd>(numVectors,1, urng);
        auto endEigenRand = std::chrono::steady_clock::now();
        auto durationEigenRand = std::chrono::duration_cast<std::chrono::microseconds>(endEigenRand - startEigenRand);
        meanEigenRand += durationEigenRand;
    }
    auto meanEigenRandMicros = meanEigenRand.count() / 1000.0;
    std::cout << flag << ",EigenRand," << meanEigenRandMicros << std::endl;

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

    file << flag << ",EigenRand," << meanEigenRandMicros << std::endl;
    file.close();

    return 0;
}