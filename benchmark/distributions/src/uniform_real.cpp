#ifdef FROM_SHELL_SCRIPT
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>

#include "../../extlibs/EigenRand/EigenRand/EigenRand"
#include "../../extlibs/Xoshiro-cpp/XoshiroCpp.hpp"
#include "../../extlibs/pcg-cpp/include/pcg_random.hpp"
#include "../../extlibs/eigen/Eigen/Dense"
#include "../../extlibs/eigen/Eigen/Core"
// #include <mkl.h>

#else
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>

#include "EigenRand/EigenRand"
#include "XoshiroCpp.hpp"
#include "pcg_random.hpp"
#include "Eigen/Dense"
#include "Eigen/Core"
#endif

using namespace XoshiroCpp;

int main(int argc, char *argv[]) {
    std::string flag;
    std::stringstream ss;

    if (argc > 1) {
        for(int i = 1; i < argc; ++i) {
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

    std::ofstream outputFile("../../results/csv/uniform_real.csv", std::ofstream::app);

    const int numVectors = 100000000;
    
    std::vector<double> stdVector(numVectors);
    std::vector<double> mersenneVector(numVectors);
    Eigen::VectorXd EigenVector(numVectors);
    std::vector<double> pcgCppVector(numVectors);
    std::vector<double> xoshiroCppVector(numVectors);
    Eigen::VectorXd EigenVector2(numVectors);

    auto startEigen = std::chrono::steady_clock::now();
    Eigen::VectorXd eigenVector = Eigen::VectorXd::Random(numVectors).array() / 2.0 + 0.5; // Remap the Eigen random from [-1,1] to [0,1]
    auto endEigen = std::chrono::steady_clock::now();
    auto durationEigen = std::chrono::duration_cast<std::chrono::milliseconds>(endEigen - startEigen);
    outputFile<< "\n" << flag << ",Eigen::VectorXd::Random,"<< durationEigen.count();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    auto startStdUniform = std::chrono::steady_clock::now();
    for (int i = 0; i < numVectors; ++i) {
        stdVector[i] = dist(gen);
    }
    auto endStdUniform = std::chrono::steady_clock::now();
    auto durationStdUniform = std::chrono::duration_cast<std::chrono::milliseconds>(endStdUniform - startStdUniform);
    outputFile<< "\n" << flag << ",std::uniform_real_distribution," << durationStdUniform.count();

    std::mt19937 mersenne(gen());
    std::uniform_real_distribution<double> mersenneDist(0.0, 1.0);

    auto startMersenne = std::chrono::steady_clock::now();
    for (int i = 0; i < numVectors; ++i) {
        mersenneVector[i] = mersenneDist(mersenne);
    }
    auto endMersenne = std::chrono::steady_clock::now();
    auto durationMersenne = std::chrono::duration_cast<std::chrono::milliseconds>(endMersenne - startMersenne);
    outputFile<< "\n" << flag << ",Mersenne Twister," << durationMersenne.count();

    Eigen::Rand::Vmt19937_64 urng{ 42 };
    auto startEigenRand = std::chrono::steady_clock::now();
    EigenVector = Eigen::Rand::uniformRealLike(EigenVector, urng);
    auto endEigenRand = std::chrono::steady_clock::now();
    auto durationEigenRand = std::chrono::duration_cast<std::chrono::milliseconds>(endEigenRand - startEigenRand);
    outputFile << "\n" << flag << ",EigenRand,"<< durationEigenRand.count();

    pcg64_fast rng;  
    rng.seed(42);  
    std::uniform_real_distribution<double> pcgDist(0.0, 1.0); 

    auto startPcgCpp = std::chrono::steady_clock::now();
    for (int i = 0; i < numVectors; ++i) {
        pcgCppVector[i] = pcgDist(rng);
    }
    auto endPcgCpp = std::chrono::steady_clock::now();
    auto durationPcgCpp = std::chrono::duration_cast<std::chrono::milliseconds>(endPcgCpp - startPcgCpp);
    outputFile << "\n" << flag << ",pcg-cpp," << durationPcgCpp.count();

    const std::uint64_t seed = 12345;
    Xoshiro256PlusPlus rg(seed);
    std::uniform_real_distribution<double> DIST(0.0,1.0);

    auto startXoshiroCpp = std::chrono::steady_clock::now();
    for (int i = 0; i < numVectors; i++) {
        xoshiroCppVector[i] =  DIST(rg);
    }
    auto endXoshiroCpp = std::chrono::steady_clock::now();
    auto durationXoshiroCpp = std::chrono::duration_cast<std::chrono::milliseconds>(endXoshiroCpp - startXoshiroCpp);
    outputFile << "\n" << flag << ",XoshiroCpp,"<< durationXoshiroCpp.count();
    outputFile.close();


    #ifdef RFM_USE_REFRAME
    // Write the results to the console
    std::cout << flag << ",Eigen::VectorXi::Random,"<< durationEigen.count() << "\n" ;
    std::cout << flag << ",std::uniform_int_distribution," << durationStdUniform.count() << "\n" ;
    std::cout << flag << ",Mersenne Twister," << durationMersenne.count() << "\n" ;
    std::cout << flag << ",EigenRanda,"<< durationEigenRand.count() << "\n" ;
    // std::cout << flag << ",Intel MKL sequential," << durationMKL.count()<< "\n" ;
    std::cout << flag << ",pcg-cpp," << durationPcgCpp.count() << "\n" ;
    std::cout << flag << ",XoshiroCpp,"<< durationXoshiroCpp.count() << "\n" ;
    #endif

    return 0;
}
