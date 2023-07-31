#include "shadingMask.hpp"
#include "perez.hpp"
#include <iostream>
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
        ( "json-filename", po::value<std::string>()->default_value(""), "json file containing the buildings names" )        
        ;

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _desc=shadingoptions
                           );

    // Load the json and recover the building names
    auto jsonfile = removeComments( readFromFile( Environment::expand( soption("json-filename" ) ) ) );  
    std::istringstream astr( jsonfile );
    json json_buildings = json::parse( astr ); 

    // Load the mesh
    auto mesh = loadMesh( _mesh = new mesh_t );

    // Compute the shading masks
    ShadingMask<mesh_t> sm(mesh,json_buildings);
    sm.computeMasks();

    std::cout << fmt::format("End of the Shading Mask Computation\n");

    // Compute the Perez sky model
    std::cout << fmt::format("Start of the Perez Sky Model Computation\n");
    double longitude = 7.75;
    double latitude = 48.57;
    PerezSkyModel psm;
    // Compute the Perez all-weather sky model for Strasbourg
    psm.computePerezSkyModel(longitude, latitude, "2020-07-22", json_buildings, 284, 40);
    std::cout << fmt::format("End of the Perez Sky Model Computation\n");

    return 0;
}