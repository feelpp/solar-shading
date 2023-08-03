#include <iostream>
#include <chrono>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#ifdef FROM_SHELL_SCRIPT
#include "../../../extlibs/EigenRand/EigenRand/EigenRand"
#else
#include "EigenRand/EigenRand"
#endif

namespace FS = std::filesystem;
using namespace Eigen;

int main(int argc, char *argv[]){
    std::chrono::microseconds meanEigenRand = std::chrono::microseconds::zero();
    std::string flag;
    std::stringstream ss;

    if (argc > 1) {
        for(int i = 1; i < argc; ++i) {
            if(i != 1) ss << " ";
            ss << argv[i];
        }
        flag = ss.str();
    }
    else {
        flag = "default";
    }

    const int numVectors = 10000;
    Rand::P8_mt19937_64 urng{ 42 };
    Rand::UniformRealGen<double> unif_gen_real{ 0, 100 };

    VectorXd EigenVector(numVectors, 1);
    VectorXi EigenVectorInt(numVectors, 1);

    for (int i = 0; i < 1000; i++)
    {
        auto startEigenRand = std::chrono::steady_clock::now();
        EigenVector = unif_gen_real.generate<VectorXd>(numVectors,1, urng);
        EigenVectorInt = EigenVector.cast<int>(); 
        auto endEigenRand = std::chrono::steady_clock::now();
        meanEigenRand += std::chrono::duration_cast<std::chrono::microseconds>(endEigenRand - startEigenRand);
    }

    double meanEigenRandMicros = meanEigenRand.count() / 1000.0;
    std::cout << flag << ",EigenRand_CAST," << meanEigenRandMicros << '\n';

    FS::path filePath = FS::current_path().parent_path().parent_path().parent_path() / "solar-shading" / "results" / "csv" / "uniform_int.csv";

    std::ofstream file;
    if (FS::exists(filePath)) file.open(filePath, std::ios_base::app);
    else {
        file.open(filePath);
        file << "Flag,Method,Duration" << '\n';
    }

    file << flag << ",EigenRand_CAST," << meanEigenRandMicros << '\n';
    file.close();

    return 0;
}
