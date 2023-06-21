#include "shadingMask_PCG.hpp"
#undef jump
#include "shadingMask_XOSHIRO.hpp"
#include "shadingMask.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feelcore/json.hpp>

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

    // Load the mesh
    auto mesh2 = loadMesh( _mesh = new mesh_t );

    auto startrngpcg = std::chrono::high_resolution_clock::now();

    // Compute the shading masks
    ShadingMaskPCG<mesh_t> sm2(mesh2,json_buildings);
    sm2.computeMasks();

    auto endrngpcg = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedrngpcg = endrngpcg - startrngpcg;

    // RNG XOSHIRO

    // Load the mesh

    auto mesh3 = loadMesh( _mesh = new mesh_t );

    auto startrngxoshiro = std::chrono::high_resolution_clock::now();

    // Compute the shading masks
    ShadingMaskXOSHIRO<mesh_t> sm3(mesh3,json_buildings);
    sm3.computeMasks();

    auto endrngxoshiro = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedrngxoshiro = endrngxoshiro - startrngxoshiro;

    std::cout << fmt::format("Elapsed time for shading mask computation: {} s\n", elapsednormal.count());
    std::cout << fmt::format("Elapsed time for shading mask computation with RNG XOSHIRO: {} s\n", elapsedrngxoshiro.count());
    std::cout << fmt::format("Elapsed time for shading mask computation with RNG PCG: {} s\n", elapsedrngpcg.count());

    // // Write the results to a csv file and append results: method + time

    // std::ofstream file;
    // file.open("../../../results/results.csv", std::ios_base::app);
    // file << "Shading mask STD RNG," << elapsednormal.count() << "\n";
    // file << "Shading mask XOSHIRO RNG," << elapsedrngxoshiro.count() << "\n";
    // file << "Shading mask PCG RNG," << elapsedrngpcg.count() << "\n";
    // file.close();

    
    std::cout << fmt::format("End Shading mask example\n");
    return 0;
}
