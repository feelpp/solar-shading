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

    // Check if command-line arguments were provided
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


    // Open the CSV file for writing
    std::ofstream outputFile("../../results/csv/uniform_int.csv", std::ofstream::app);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parameters Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    const int numVectors = 100000000;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Vector Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::vector<int> stdVector(numVectors);
    std::vector<int> mersenneVector(numVectors);
    Eigen::VectorXi EigenVector(numVectors);
    std::vector<int> pcgCppVector(numVectors);
    std::vector<int> xoshiroCppVector(numVectors);
    Eigen::VectorXi EigenVector2(numVectors);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Technique 1: Eigen::VectorXi::Random() ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    auto startEigen = std::chrono::steady_clock::now();

    Eigen::VectorXi eigenVector = Eigen::VectorXi::Random(numVectors);

    auto endEigen = std::chrono::steady_clock::now();
    auto durationEigen = std::chrono::duration_cast<std::chrono::milliseconds>(endEigen - startEigen);

    // Write the results to the CSV file
    outputFile<< "\n" << flag << ",Eigen::VectorXi::Random,"<< durationEigen.count();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Technique 2: std::uniform_int_distribution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 100);

    auto startStdUniform = std::chrono::steady_clock::now();

    for (int i = 0; i < numVectors; ++i) {
        stdVector[i] = dist(gen);
    }

    auto endStdUniform = std::chrono::steady_clock::now();
    auto durationStdUniform = std::chrono::duration_cast<std::chrono::milliseconds>(endStdUniform - startStdUniform);

    outputFile<< "\n" << flag << ",std::uniform_int_distribution," << durationStdUniform.count();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Technique 3: Mersenne Twister ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    std::mt19937 mersenne(gen());
    std::uniform_int_distribution<int> mersenneDist(0, 100);

    auto startMersenne = std::chrono::steady_clock::now();

    for (int i = 0; i < numVectors; ++i) {
        mersenneVector[i] = mersenneDist(mersenne);
    }

    auto endMersenne = std::chrono::steady_clock::now();
    auto durationMersenne = std::chrono::duration_cast<std::chrono::milliseconds>(endMersenne - startMersenne);

    outputFile<< "\n" << flag << ",Mersenne Twister," << durationMersenne.count();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Technique 4a: Eigen::Rand ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Eigen::Rand::Vmt19937_64 urng{ 42 };

    auto startEigenRand = std::chrono::steady_clock::now();

    EigenVector = Eigen::Rand::uniformInt<Eigen::VectorXi>(numVectors, 1, urng, 0, 100);

    auto endEigenRand = std::chrono::steady_clock::now();
    auto durationEigenRand = std::chrono::duration_cast<std::chrono::milliseconds>(endEigenRand - startEigenRand);

    outputFile << "\n" << flag << ",EigenRanda,"<< durationEigenRand.count();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Technique 4b: Eigen::Rand ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Eigen::Rand::Vmt19937_64 urng2{ 42 };

    auto startEigenRand2 = std::chrono::steady_clock::now();

    EigenVector2 = Eigen::Rand::uniformIntLike(EigenVector2, urng2, 0, 100);

    auto endEigenRand2 = std::chrono::steady_clock::now();
    auto durationEigenRand2 = std::chrono::duration_cast<std::chrono::milliseconds>(endEigenRand2 - startEigenRand2);

    outputFile << "\n" << flag << ",EigenRandb,"<< durationEigenRand2.count();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Technique 5: Intel MKL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // auto startMKL = std::chrono::steady_clock::now();

    // std::vector<int> mklVector(numVectors);
    // VSLStreamStatePtr stream;
    // vslNewStream(&stream, VSL_BRNG_MCG31, 1);

    // viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, numVectors, mklVector.data(), 0, 100);

    // auto endMKL = std::chrono::steady_clock::now();
    // auto durationMKL = std::chrono::duration_cast<std::chrono::milliseconds>(endMKL - startMKL);

    // outputFile << "\n" << flag << ",Intel MKL sequential," << durationMKL.count();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Technique 6: pcg-cpp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
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

    // Write the results to the CSV file
    outputFile << "\n" << flag << ",pcg-cpp," << durationPcgCpp.count();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Technique 7: Xoshiro-cpp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const std::uint64_t seed = 12345;

    Xoshiro256PlusPlus rg(seed);

    std::uniform_int_distribution<int> DIST(0,100);

    auto startXoshiroCpp = std::chrono::steady_clock::now();

    for (int i = 0; i < numVectors; i++) {
        xoshiroCppVector[i] =  DIST(rg);
    }

    auto endXoshiroCpp = std::chrono::steady_clock::now();
    auto durationXoshiroCpp = std::chrono::duration_cast<std::chrono::milliseconds>(endXoshiroCpp - startXoshiroCpp);

    outputFile << "\n" << flag << ",XoshiroCpp,"<< durationXoshiroCpp.count();

    // Close the CSV file
    outputFile.close();

    // If the script is launched with reframe, write the results to the console
    #ifdef RFM_USE_REFRAME
    std::cout << flag << ",Eigen::VectorXi::Random,"<< durationEigen.count() << "\n" ;
    std::cout << flag << ",std::uniform_int_distribution," << durationStdUniform.count() << "\n" ;
    std::cout << flag << ",Mersenne Twister," << durationMersenne.count() << "\n" ;
    std::cout << flag << ",EigenRanda,"<< durationEigenRand.count() << "\n" ;
    std::cout << flag << ",EigenRandb,"<< durationEigenRand2.count() << "\n" ;
    // std::cout << flag << ",Intel MKL sequential," << durationMKL.count()<< "\n" ;
    std::cout << flag << ",pcg-cpp," << durationPcgCpp.count() << "\n" ;
    std::cout << flag << ",XoshiroCpp,"<< durationXoshiroCpp.count() << "\n" ;
    #endif
    
    return 0;
}
