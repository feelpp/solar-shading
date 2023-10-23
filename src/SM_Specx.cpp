#include <iostream>

#include <shading_mask_specx.hpp>




#include "SpDataAccessMode.hpp"
#include "Utils/SpUtils.hpp"

#include "Task/SpTask.hpp"
#include "Legacy/SpRuntime.hpp"
#include "Utils/SpTimer.hpp"
#include "Utils/small_vector.hpp"




#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feelcore/json.hpp>
//#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>


//~/SpecxProjects/specxpp/build/Benchmark/specx-solar-shading$
//mpirun -np 1 

//./feelpp_ss_SM_Specx --config-file cases/19_buildings_strasbourg/19_buildings_strasbourg.cfg



using namespace Feel;
//using namespace std;

int main(int argc, char **argv)
{


   const int NumThreads = SpUtils::DefaultNumThreads();
   //SpRuntime runtime(NumThreads);
   //std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
   //std::cout<<"Nb Thread="<<runtime.getNbThreads()<<"\n";
   //std::cout<<"Nb CPU Workers="<<runtime.getNbCpuWorkers()<<"\n";
   //std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
 

    


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
    std::cout<<"\n";
    std::cout<<"[INFO] Read envrironment =>"<<"["<<Environment::expand( soption("json-filename" ) )<<"]"<<"\n";
    std::cout<<"\n";

    auto jsonfile = removeComments( readFromFile( Environment::expand( soption("json-filename" ) ) ) );
    std::istringstream astr( jsonfile );
    json json_buildings = json::parse( astr );

    // Load the mesh
    auto mesh = loadMesh( _mesh = new mesh_t );

    // Compute the shading masks
    
    std::cout<<"[INFO] : START SHADING MASK"<<"\n";
    ShadingMask<mesh_t> sm(mesh,json_buildings);

    std::cout<<"[INFO] : BEGIN COMPUTE"<<"\n";

    sm.computeMasks();

    std::cout<<"[INFO] : END COMPUTE"<<"\n";


    std::cout << fmt::format("End Shading mask example\n");

    //runtime.stopAllThreads();

    //runtime.generateDot("RuntimeShadingMask.dot",true);
    //runtime.generateTrace("RuntimeShadingMask.svg");

    
    return 0;
}
