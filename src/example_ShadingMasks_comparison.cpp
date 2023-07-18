#include "shadingMask_PCG.hpp"
#include "shadingMask_EigenRand_vectorized.hpp"
#include "shadingMask_EigenRand_cast.hpp"
// #undef jump
#include "shadingMask_XOSHIRO.hpp"
#include "shadingMask.hpp"
#include "shadingMask_mix.hpp"
#include "shadingMask_EigenRand.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feelcore/json.hpp>
#include <experimental/filesystem>

namespace FS = std::filesystem;
using namespace Feel;

int main(int argc, char **argv)
{
    using namespace Feel;
    using mesh_t = Mesh<Simplex<3>>;

    po::options_description shadingoptions( "Shading mask options" );
    shadingoptions.add_options()
        ( "json-filename", po::value<std::string>()->default_value(""), "json file containing the buildings names" );

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _desc=shadingoptions
                           );
    // Load the json and recover the building names
    auto jsonfile = removeComments( readFromFile( Environment::expand( soption("json-filename" ) ) ) );  
    std::istringstream astr( jsonfile );
    json json_buildings = json::parse( astr ); 

    // STD RNG

    // Load the mesh
    auto mesh = loadMesh( _mesh = new mesh_t );


    auto startnormal = std::chrono::high_resolution_clock::now();

    // Compute the shading masks
    ShadingMask<mesh_t> sm(mesh,json_buildings);
    sm.computeMasks();

    auto endnormal = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsednormal = endnormal - startnormal;

    // RNG PCG

    auto startrngpcg = std::chrono::high_resolution_clock::now();

    // Compute the shading masks
    ShadingMaskPCG<mesh_t> sm2(mesh,json_buildings);
    sm2.computeMasks();

    auto endrngpcg = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedrngpcg = endrngpcg - startrngpcg;

    // RNG XOSHIRO

    auto startrngxoshiro = std::chrono::high_resolution_clock::now();

    // Compute the shading masks
    ShadingMaskXOSHIRO<mesh_t> sm3(mesh,json_buildings);
    sm3.computeMasks();

    auto endrngxoshiro = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedrngxoshiro = endrngxoshiro - startrngxoshiro;

    // RNG MIX

    auto startrngmix = std::chrono::high_resolution_clock::now();

    // Compute the shading masks
    ShadingMaskMIX<mesh_t> sm4(mesh,json_buildings);
    sm4.computeMasks();

    auto endrngmix = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedrngmix = endrngmix - startrngmix;

    // RNG EIGEN

    auto startrngeigen = std::chrono::high_resolution_clock::now();

    // Compute the shading masks
    ShadingMaskER<mesh_t> sm5(mesh,json_buildings);
    sm5.computeMasks();

    auto endrngeigen = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedrngeigen = endrngeigen - startrngeigen;

    // RNG EIEN CAST

    auto startrngeigencast = std::chrono::high_resolution_clock::now();

    // Compute the shading masks
    ShadingMaskERC<mesh_t> sm6(mesh,json_buildings);
    sm6.computeMasks();

    auto endrngeigencast = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedrngeigencast = endrngeigencast - startrngeigencast;

    // RNG EIGEN VECTORIZED

    auto startrngeigenvect = std::chrono::high_resolution_clock::now();

    // Compute the shading masks
    ShadingMaskERV<mesh_t> sm7(mesh,json_buildings);
    sm7.computeMasks();

    auto endrngeigenvect = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedrngeigenvect = endrngeigenvect - startrngeigenvect;

    std::cout << fmt::format("Elapsed time for shading mask computation with RNG STD:{} s\n", elapsednormal.count());
    std::cout << fmt::format("Elapsed time for shading mask computation with RNG XOSHIRO:{} s\n", elapsedrngxoshiro.count());
    std::cout << fmt::format("Elapsed time for shading mask computation with RNG PCG:{} s\n", elapsedrngpcg.count());
    std::cout << fmt::format("Elapsed time for shading mask computation with RNG MIX:{} s\n", elapsedrngmix.count());
    std::cout << fmt::format("Elapsed time for shading mask computation with RNG EIGEN:{} s\n", elapsedrngeigen.count());
    std::cout << fmt::format("Elapsed time for shading mask computation with RNG EIGEN CAST:{} s\n", elapsedrngeigencast.count());
    std::cout << fmt::format("Elapsed time for shading mask computation with RNG EIGEN VEC:{} s\n", elapsedrngeigenvect.count());

    // Write the results to a csv file and append results: method + time

    std::ofstream file;
    FS::path scriptPath = FS::current_path();
    FS::path scriptDirectory = scriptPath.parent_path();
    FS::path rootDirectory1 = scriptDirectory.parent_path();
    FS::path rootDirectory = rootDirectory1.parent_path();
    FS::path filePath = rootDirectory / "solar-shading" / "results" / "csv" / "results.csv";



    std::cout << fmt::format("Writing results to file: {}\n", filePath.string());

    file.open(filePath.string(), std::ios_base::app);

    if (!file) {
        // Handle file opening failure
        std::cerr << "Error opening file: " << filePath << std::endl;
    } else {
        // File opened successfully
        file << "\n" << "Shading mask STD RNG," << elapsednormal.count();
        file << "\n" << "Shading mask PCG RNG," << elapsedrngpcg.count();
        file << "\n" << "Shading mask XOSHIRO RNG," << elapsedrngxoshiro.count();
        file << "\n" << "Shading mask MIX RNG," << elapsedrngmix.count();
        file << "\n" << "Shading mask EIGEN RNG," << elapsedrngeigen.count();
        file << "\n" << "Shading mask EIGEN CAST RNG," << elapsedrngeigencast.count();
        file << "\n" << "Shading mask EIGEN VECTORIZED RNG," << elapsedrngeigenvect.count();
        file.close();
    }
    
    std::cout << fmt::format("End Shading mask example\n");
    return 0;
}
