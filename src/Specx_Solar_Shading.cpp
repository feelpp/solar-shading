#include <iostream>

#include <Specx_Shading_mask.hpp>

#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feelcore/json.hpp>
//#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>


#include "SpDataAccessMode.hpp"
#include "Utils/SpUtils.hpp"

#include "Task/SpTask.hpp"
#include "Legacy/SpRuntime.hpp"
#include "Utils/SpTimer.hpp"
#include "Utils/small_vector.hpp"



using namespace Feel;

int main(int argc, char **argv)
{

 const int NumThreads = SpUtils::DefaultNumThreads();

    std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
    std::cout<<"[WELCOME SPECX ]"<<"\n";
    std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n\n";
    std::cout<<"[SPECX INFO] : Nb Thread="<<SpUtils::DefaultNumThreads()<<"\n";

    using namespace Feel;
    using mesh_t = Mesh<Simplex<FEELPP_TOP_DIM,1,3>>;

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

    std::cout<<"[SPECX INFO] : Load Mesh.\n";
    // Load the mesh
    //tic();
    //auto mesh = loadMesh( _mesh = new mesh_t, _update=MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED );
    //toc("Mesh loaded");

    #if FEELPP_TOP_DIM==2
    auto mesh = loadMesh( _mesh = new mesh_t, _update=MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED );
    #else
    auto mesh = loadMesh( _mesh = new mesh_t );
    #endif


    // Compute the shading masks
    std::cout<<"\n\n";
    std::cout<<"[SPECX INFO] : Compute the shading masks.\n";
    ShadingMask<mesh_t> sm(mesh,json_buildings);
    std::cout<<"\n\n";
    std::cout<<"[SPECX INFO] : Compute...\n";
    sm.computeMasks();
    std::cout<<"\n\n";
    std::cout<<"[SPECX INFO] : END !!!\n\n\n";

    //std::cout << fmt::format("End Shading mask example\n");
    return 0;
}
