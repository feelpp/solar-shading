
#include <future>
#include <feel/feelmesh/bvh.hpp>



#include "SpDataAccessMode.hpp"
#include "Utils/SpUtils.hpp"

#include "Task/SpTask.hpp"
#include "Legacy/SpRuntime.hpp"
#include "Utils/SpTimer.hpp"
#include "Utils/small_vector.hpp"

int GetCoorTask(int NbObjects,int NbTasks) 
{
    return(round((float(NbObjects)/float(NbTasks)-float(NbObjects/NbTasks))*float(NbTasks)));
}

int GetNbLoop(int NbObjects,int NbTasks)
{
    return(std::max(1,NbObjects/NbTasks));
}



namespace Feel {

template <typename MeshType>
class ShadingMask
{
    using mesh_type = MeshType;
    typedef typename MeshType::ptrtype mesh_ptrtype;

    using tr_mesh_type = typename std::conditional<MeshType::nDim==MeshType::nRealDim,
                                                typename MeshType::trace_mesh_type,
                                                typename MeshType::type >::type ;
    using tr_mesh_ptrtype = typename std::conditional<MeshType::nDim==MeshType::nRealDim,
                                                typename MeshType::trace_mesh_ptrtype,
                                                typename MeshType::ptrtype >::type ;
    // typedef typename MeshType::trace_mesh_ptrtype tr_mesh_ptrtype;
    typedef typename matrix_node<double>::type matrix_node_type;

public:
    using value_type = double;

    

    ShadingMask(mesh_ptrtype mesh, nl::json const& specs, int intervalsAzimuth=72, int intervalsAltitude=10 )
    {
        // Read the number of rays per triangle and the number of threads
        j_ = specs;
        M_Nrays = specs["Nrays"];
        M_Nthreads = specs["Nthreads"].get<int>() ;

        bool QActiveSpecx=true;
 

        // For each building, save the surface mesh and build the corresponding BVH tree for ray search
        if constexpr( MeshType::nDim==MeshType::nRealDim )
        {

            if( specs["/Buildings"_json_pointer].contains("list") ) // the list of volume markers is provided
            {
                auto markersVolume = specs["Buildings"]["list"].get<std::vector<std::string>>();

                
                bool QActiveSpecx=true;
                if (QActiveSpecx) { 

                    SpTimer timer;
                    auto mV = j_["Buildings"]["list"].get<std::vector<std::string>>();
                    const int NbElements=mV.size();
                    std::cout<<"NbElements="<<NbElements<<"\n";
                    std::cout<<"NbElements="<<mV<<"\n";
                    const int NumThreads = SpUtils::DefaultNumThreads();
                    int NbUsed=std::min(NbElements,NumThreads);
                    SpRuntime runtime_init(NbUsed);

                    std::cout<<"[INFO SPECS] : NbThreads="<<NumThreads<<" NbUsed="<<NbUsed<<"\n";

                    

                    for(std::string buildingName : markersVolume)
                    {
                        runtime_init.task(
                        SpRead(buildingName),
                        [&](const std::string & bn)->bool {

                            std::cout << fmt::format("{}\n",bn);
                            auto volumeSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,bn));
                            auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));

                            auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                            M_bvh_tree_vector.insert(std::make_pair(bn, std::move(bvhBuilding) ));
                            M_submeshes.insert(std::make_pair(bn, surfaceSubmesh ));
                            //usleep(1000);
                            return true;
                        }
                        ).setTaskName("Task 0 :"+buildingName);
                    }

                    runtime_init.waitAllTasks();
                    runtime_init.stopAllThreads();
                    timer.stop();
                    std::cout << "[INFO SPECS] : Task Read Time = "<<timer.getElapsed()<< " s" << std::endl; 
                    runtime_init.generateDot("RuntimeReadList.dot",true);
                    runtime_init.generateTrace("RuntimeReadList.svg");

                }
                else 
                {
                    auto mV = j_["Buildings"]["list"].get<std::vector<std::string>>();
                    int NbElement=mV.size();
                    std::cout<<"NbElements="<<NbElement<<"\n";
                    std::cout<<"NbElements="<<mV<<"\n";

                    for(std::string buildingName : markersVolume)
                    {
                        std::cout << fmt::format("{}\n",buildingName);
                        auto volumeSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                        auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));

                        auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                        M_bvh_tree_vector.insert(std::make_pair( buildingName , std::move(bvhBuilding) ));

                        M_submeshes.insert(std::make_pair( buildingName , surfaceSubmesh ));

                    }
                }
            }
            else if( specs["/Buildings"_json_pointer].contains("fileVolumes")) // a csv containing the volume markers is provided
            {
                std::string buildingName;
                // open the file

                std::ifstream fileVolumes(Environment::expand(specs["Buildings"]["fileVolumes"].get<std::string>()));

                auto mV = j_["Buildings"]["fileVolumes"].get<std::string>();
                std::cout<<"Name File"<<mV<<"\n";
                
                std::ifstream FICH(Environment::expand(specs["Buildings"]["fileVolumes"].get<std::string>()));
                std::string lineA;
                int NbObjects=0;
                std::vector<std::string> ListObjects;
                while (getline(FICH, lineA)) {
                    //std::cout << lineA << endl;
                    ListObjects.push_back(lineA);
                    NbObjects++;
                }
                FICH.close();
                std::cout<<"NbObjects="<<NbObjects<<"\n";

                for(int k = 0 ; k < NbObjects ; ++k){
                    std::cout<<ListObjects[k]<<"\n";
                }
                std::cout<<"******************************\n";
                
                bool QActiveSpecx=true;
                if (QActiveSpecx) { 
                    SpTimer timer;
                    
                    int NumThreads = std::min(5,std::min(NbObjects,SpUtils::DefaultNumThreads()));
                    SpRuntime runtime0(NumThreads);
                    int NbTasksToSubmit = runtime0.getNbThreads();

                    int NbLoops=GetNbLoop(NbObjects,NbTasksToSubmit);
                    std::cout<<"Nb Loops F="<<float(NbObjects)/float(NbTasksToSubmit)<<"\n";
                    std::cout<<"Nb Loops="<<NbLoops<<"\n";
                    int NbCoor=GetCoorTask(NbObjects,NbTasksToSubmit);
                    std::cout<<"Nb Coor="<<std::fixed<<std::setprecision(9)<<NbCoor<<"\n";

                    int Nbidx=NbTasksToSubmit;
                    std::cout<<"Nb idx Th="<<Nbidx<<"\n";

                    //small_vector<SpAbstractTaskWithReturn<double>::SpTaskViewer> elapsed;
                    //elapsed.reserve(NbTasksToSubmit*NbLoops);

                    int initVal = 0;
                    int Nb=0;
                    int NbidxN=Nbidx;
                    for(int idxLoop = 0 ; idxLoop < NbLoops ; ++idxLoop){
                            if (idxLoop==NbLoops-1) { 
                                Nbidx=Nbidx+NbCoor;
                            }
                            for(int idx = 0 ; idx < Nbidx ; ++idx){
                                    int Index=idxLoop*NbidxN+idx;
                                    std::cout <<"Index:"<<Index<<"\n";

                                    //std::promise<int> promise1;

                                    
                                    auto returnValue=runtime0.task(SpRead(initVal),
                                        [&](const int&)-> bool {

                                    //auto returnValue=runtime0.task(SpWrite(Nb),
                                    //    [&](int& NbWrite)-> bool {
                                    //      buildingName=ListObjects[NbWrite];
                                    //      std::cout <<"Building Name (V):"<<NbWrite<< fmt::format("{}\n",buildingName);
                                    //      NbWrite=NbWrite+1;

                                          
                                          buildingName=ListObjects[Index];
                                          //std::cout <<"Building Name (V):"<<Index<< fmt::format("{}\n",buildingName);
 
                                          auto volumeSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                                          auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));

                                          auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                                          M_bvh_tree_vector.insert(std::make_pair(buildingName,std::move(bvhBuilding) ));
                                          M_submeshes.insert(std::make_pair(buildingName,surfaceSubmesh ));
                                          return true;
                                        }
                                    );   

                                    //promise1.set_value(0);
                                    returnValue.wait();
                                    
                                    
                                Nb++;
                            }
                            runtime0.waitAllTasks(); 
                    }
                    std::cout<<"CTRL Nb="<<Nb<<"\n";
                   
                    runtime0.waitAllTasks();
                    //runtime0.stopAllThreads();
                    timer.stop();
                    std::cout << "[INFO SPECS] : Task Read Time = "<<timer.getElapsed()<< " s" << std::endl; 
                    runtime0.generateDot("RuntimeReadVolumes.dot",true);
                    runtime0.generateTrace("RuntimeReadVolumes.svg");

                    //usleep(100000);
                    //getchar();
                    



                    //const int NumThreads = SpUtils::DefaultNumThreads();


                    /*
                    SpRuntime runtime_init2(NumThreads);


                    for(int ido = 0 ; ido < NbObjects ; ++ido)
                    {   
                        getline(fileVolumes,buildingName); 
                    
                        auto returnValue=runtime_init2.task(
                                SpRead(buildingName),
                                [&](const std::string & bn)->bool {
                                    std::cout <<"Building Name (V):"<< fmt::format("{}\n",bn);

                                    auto volumeSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,bn));
                                    auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));

                                    auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                                    M_bvh_tree_vector.insert(std::make_pair(bn,std::move(bvhBuilding) ));
                                    M_submeshes.insert(std::make_pair(bn,surfaceSubmesh ));
                                    
                                    return true;
                                }
                        );
                        //returnValue.setTaskName("Task Read Volume :"+buildingName);


                        //returnValue.wait();
                        runtime_init2.waitAllTasks();
                        //i++;
                    
                    }

                    runtime_init2.waitAllTasks();
                    //runtime_init.waitRemain(19);
                    //runtime_init2.stopAllThreads();
                    timer.stop();
                    std::cout << "[INFO SPECS] : Task Read Time = "<<timer.getElapsed()<< " s" << std::endl; 
                    runtime_init2.generateDot("RuntimeReadVolumes.dot",true);
                    runtime_init2.generateTrace("RuntimeReadVolumes.svg");

                    usleep(100000);
                    */
                    
                }

                else 
                
                
                {
                    auto mV = j_["Buildings"]["fileVolumes"].get<std::string>();
                    int NbElement=mV.size();
                    std::cout<<"NbElement="<<NbElement<<"\n";
                    std::cout<<"NbElement="<<mV<<"\n";
                    
                    while ( getline(fileVolumes,buildingName) )
                    {
                        std::cout <<"Building Name:"<< fmt::format("{}\n",buildingName);
                        auto volumeSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                        auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));
                        auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                        M_bvh_tree_vector.insert(std::make_pair( buildingName , std::move(bvhBuilding) ));

                        M_submeshes.insert(std::make_pair( buildingName , surfaceSubmesh ));

                    }
                    
                }

            }
        }
        else
        {
            if( specs["/Buildings"_json_pointer].contains("fileSurfaces") ) // a csv containing the surface markers is provided
            {
                std::string buildingName;
                std::ifstream fileSurfaces(Environment::expand(specs["Buildings"]["fileSurfaces"].get<std::string>()));
                std::cout << Environment::expand(specs["Buildings"]["fileSurfaces"].get<std::string>()) << std::endl;
                // read, line by line, the building marker

                /*
                bool QActiveSpecx=true;
                if (QActiveSpecx) { 
                    const int NumThreads = SpUtils::DefaultNumThreads();
                    SpRuntime runtime_init(NumThreads);
                    
                    while ( getline(fileSurfaces,buildingName) )
                    {
                        std::promise<int> promise0;

                        runtime_init.task(SpRead(buildingName), [&promise0](const std::string & bn){
                            promise0.get_future().get();
                        }).setTaskName("Start Task");

                        runtime_init.task(
                                SpRead(buildingName),
                                [&](const std::string & bn){
                                    std::cout << fmt::format("{}\n",bn);
                                    auto surfaceSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,bn));
                                    auto listMarkers = surfaceSubmesh->markerNames();
                                    // Delete the marker associated to the building
                                    // to Keep only face markers
                                    auto it = listMarkers.find(bn);
                                    listMarkers.erase(it);
                                    surfaceSubmesh->setMarkerNames(listMarkers);

                                    auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                                    M_bvh_tree_vector.insert(std::make_pair( buildingName , std::move(bvhBuilding) ));

                                    M_submeshes.insert(std::make_pair(bn , surfaceSubmesh ));
                                    usleep(1000);
                                }
                        ).setTaskName("Task 0 :"+buildingName);

                        promise0.set_value(0);
                    }

                    runtime_init.waitAllTasks();
                    runtime_init.stopAllThreads();
                    runtime_init.generateDot("RuntimeInitSurfaces.dot",true);
                    runtime_init.generateTrace("RuntimeInitSurfaces.svg");


                } 
                else */
                {
                    while ( getline(fileSurfaces,buildingName) )
                    {
                        std::cout << fmt::format("{}\n",buildingName);
                        auto surfaceSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                        auto listMarkers = surfaceSubmesh->markerNames();
                        // Delete the marker associated to the building
                        // to Keep only face markers
                        auto it = listMarkers.find(buildingName);
                        listMarkers.erase(it);
                        surfaceSubmesh->setMarkerNames(listMarkers);

                        auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                        M_bvh_tree_vector.insert(std::make_pair( buildingName , std::move(bvhBuilding) ));

                        M_submeshes.insert(std::make_pair( buildingName , surfaceSubmesh ));

                }
                }
            }
        }


        fixAzimuthAltitudeDiscretization(intervalsAzimuth, intervalsAltitude);


        std::random_device rd;
        std::random_device rd2;
        std::mt19937 gen(rd());
        std::mt19937 gen2(rd2());
        gen.seed(std::chrono::high_resolution_clock::now()
                            .time_since_epoch()
                            .count());
        gen2.seed(std::chrono::high_resolution_clock::now()
                            .time_since_epoch()
                            .count());

        M_gen=gen;
        M_gen2=gen2;
    }



    // Subdivide the azimuth angles [0,360]° and altitude angles [0,90]° in subsets for easier computation of the shading masks
    void fixAzimuthAltitudeDiscretization(int intervalsAzimuth=72, int intervalsAltitude=10)
    {
        M_azimuthSize = intervalsAzimuth;
        M_altitudeSize = intervalsAltitude;
        value_type deltaAzimuth = 2* M_PI /intervalsAzimuth;
        M_azimuthAngles.resize(M_azimuthSize);
        for(int i=0; i<intervalsAzimuth; i++)
        {
            M_azimuthAngles[i] = i * deltaAzimuth+1e-6;
        }

        value_type deltaAltitude = 0.5 * M_PI /intervalsAltitude;
        M_altitudeAngles.resize(M_altitudeSize);
        for(int i=0; i<intervalsAltitude; i++)
        {
            M_altitudeAngles[i] = i * deltaAltitude;
        }

    }

    // Choose a random pair of indices in the discretized azimuth and altitude vectors
    void getRandomDirectionSM(std::vector<double> &random_direction, std::mt19937 & M_gen, std::mt19937 & M_gen2, int& index_azimuth, int& index_altitude)
    {
        std::uniform_int_distribution<int> dist_azimuth(0,M_azimuthSize-1);
        std::uniform_int_distribution<int> dist_altitude(0,M_altitudeSize-1);

        int size = random_direction.size();

        if(random_direction.size()==3)
        {
            index_azimuth = dist_azimuth(M_gen);
            index_altitude = dist_altitude(M_gen2);
            double phi = -( M_azimuthAngles[index_azimuth] ) + M_PI*0.5 ; // recover spherical coordinate from azimuth angle
            double theta = M_PI*0.5 - M_altitudeAngles[index_altitude]; // recover spherical coordinate from altitude

            random_direction[0]=math::sin(theta)*math::cos(phi);
            random_direction[1]=math::sin(theta)*math::sin(phi);
            random_direction[2]=math::cos(theta);

        }
        else
        {
            throw std::logic_error( "Wrong dimension " + std::to_string(random_direction.size()) + " for the random direction" );
        }

    }

    Eigen::VectorXd get_random_point(matrix_node_type const& element_points)
    {
        int dimension;

        dimension = column(element_points, 0).size();

        if(dimension==3)
        {
            // Choose points in a parallelogram, uniformly
            Eigen::VectorXd p1(dimension),p2(dimension),p3(dimension),v(dimension),u(dimension),p(dimension);
            for(int i=0;i<3;i++)
            {
                p1(i)=column(element_points, 0)[i];
                p2(i)=column(element_points, 1)[i];
                p3(i)=column(element_points, 2)[i];
            }
            v = p2-p1;
            u = p3-p1;
            while(true)
            {
                unsigned seed2 = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                unsigned seed3 = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                std::default_random_engine generator3(seed2),generator4(seed3);
                std::uniform_real_distribution<double> xi1(0,1),xi2(0,1);
                double s = xi1(generator3);
                double t = xi2(generator4);
                // If the point is on the left of the diagonal, keep it, else take the symmetric one
                bool in_triangle = (s + t <= 1);
                if(in_triangle)
                    p = p1 + s * u + t * v;
                else
                    p= p1 + (1 - s) * u + (1 - t) * v;

                if (isOnSurface(p,p1,p2,p3))
                    return p;
                else
                {
                    throw std::logic_error("Point not on triangle, but it must be");
                    return p1;
                }
            }
        }
        else
        {
            Eigen::VectorXd p1(dimension);
            for(int i=0;i<dimension;i++)
            {
                p1(i)=column(element_points, 0)[i];
            }
            throw std::logic_error( "Problem in the computation of the random point" );
            return p1;
        }

    }

    // 3D case
    // Compute the sum of the areas of three subtriangles
    double elementArea(Eigen::VectorXd const& point,Eigen::VectorXd const& el_p1,Eigen::VectorXd const& el_p2,Eigen::VectorXd const& el_p3)
    {
        double area = 0.;
        if(point.size()==3)
        {
            Eigen::Vector3d point_3d,el_p1_3d,el_p2_3d,el_p3_3d;
            point_3d << point[0], point[1],point[2];
            el_p1_3d << el_p1[0], el_p1[1],el_p1[2];
            el_p2_3d << el_p2[0], el_p2[1],el_p2[2];
            el_p3_3d << el_p3[0], el_p3[1],el_p3[2];

            auto v1 = (point_3d-el_p1_3d).cross(point_3d-el_p2_3d);
            auto v2 = (point_3d-el_p2_3d).cross(point_3d-el_p3_3d);
            auto v3 = (point_3d-el_p3_3d).cross(point_3d-el_p1_3d);
            area = math::sqrt(v1.dot(v1))/2. + math::sqrt(v2.dot(v2))/2. + math::sqrt(v3.dot(v3))/2. ;
            return area;
        }
        else
        {
            throw std::logic_error( "Wrong area calculation for the random direction" );
            return -1.;
        }
    }

    // 3D case
    // Compare the area of the 2d simplex V1V2V3 (as sum of 3 subtriangles V_iV_jB) and the 2d triangle
    // created by the intersection P of the ray with the plane of  V1V2V3 (as sum of V_iV_jP)
    bool isOnSurface(Eigen::VectorXd const &point,Eigen::VectorXd const &el_p1,Eigen::VectorXd const &el_p2,Eigen::VectorXd const &el_p3)
    {
        auto c = (el_p1+el_p2+el_p3)/3.;

        auto elem_area = elementArea(c, el_p1,el_p2,el_p3);
        auto area = elementArea(point, el_p1,el_p2,el_p3);
        if (math::abs(area-elem_area)/area<1e-5)
            return true;
        else
            return false;
    }

    // Compute shading masks for the buildings in the json file
    void computeMasks()
    {
       
        bool QActiveSpecx=true; 
        bool QViewInfoSpecs=true;
        if (QViewInfoSpecs) {
            std::cout << "[INFO SPECS] : COMPUTE" << std::endl;
            std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
        }
        
        int NbThreadUsed;
        int NbElements;

        const int NumThreads = SpUtils::DefaultNumThreads();
       
        if( j_["/Buildings"_json_pointer].contains("list") ) // the list of volume markers is provided
        {

            auto mV = j_["Buildings"]["list"].get<std::vector<std::string>>();
            std::cout<<"NbElements="<<mV<<"\n";
            NbElements=mV.size();
        } 
        else if( j_["/Buildings"_json_pointer].contains("fileVolumes"))
        {
            auto mV = j_["Buildings"]["fileVolumes"].get<std::string>();
            std::cout<<"NbElements="<<mV<<"\n";
            NbElements=mV.size();
            NbElements=NumThreads;           
        }

        
        std::cout<<"NbElements="<<NbElements<<"\n";
        NbThreadUsed=std::min(NbElements,NumThreads);
        SpRuntime runtime_init(NbThreadUsed);

        std::cout<<"[INFO SPECS] : Nb Thread="<<NumThreads<<" NbUsed="<<NbThreadUsed<<"\n";

        SpRuntime runtime_Building(NbThreadUsed);

        SpTimer timer;
        int SpNbTasksToSubmit=0;
        if (QViewInfoSpecs) {
            std::cout<<"[INFO SPECS] : Nb Thread="<<runtime_Building.getNbThreads()<<"\n";
            std::cout<<"[INFO SPECS] : Nb CPU Workers="<<runtime_Building.getNbCpuWorkers()<<"\n";
            std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
        }
 
        



        if( j_["/Buildings"_json_pointer].contains("list") ) // the list of volume markers is provided
        {


            if (QViewInfoSpecs) { std::cout << "[INFO SPECS] : COMPUTE 1" << std::endl; }
            auto markersVolume = j_["Buildings"]["list"].get<std::vector<std::string>>();

            if (QActiveSpecx) {  
                
                if (QViewInfoSpecs) { std::cout << "[INFO SPECS] : START LIST" << std::endl; }

                for(std::string building_name : markersVolume)
                {
                    SpNbTasksToSubmit++;
                    if (QViewInfoSpecs) { std::cout << "[INFO SPECS] : NAME =>"; }
                    runtime_Building.task(
                        SpPriority(1),
                        SpRead(building_name),
                        [&](const std::string & bn){
                        //computeMasksOneBuilding(building_name);
                        computeMasksOneBuilding(bn);
                        //computeMasksOneBuildingSpecx(bn);
                        if (QViewInfoSpecs) { std::cout<<bn<<":"<<building_name<<"\n"; }
                        usleep(1000);
                    }
                    ).setTaskName("Task 1S Compute:"+building_name);
                }
                runtime_Building.waitAllTasks();
             }
             else
            {
                for(std::string building_name : markersVolume)
                {
                    std::cout << "[INFO NO SPECS] : Name =>"<<building_name << std::endl;
                    computeMasksOneBuilding(building_name);
                }
            }

             
        }
        else if( j_["/Buildings"_json_pointer].contains("fileVolumes")) // a csv containing the volume markers is provided
        {
            if (QViewInfoSpecs) { std::cout << "[INFO SPECS] : COMPUTE 2" << std::endl; }
            std::string building_name;
            std::ifstream fileVolumes(Environment::expand(j_["Buildings"]["fileVolumes"].get<std::string>()));


            if (QActiveSpecx) {  
                if (QViewInfoSpecs) { std::cout << "[INFO SPECS] : START VOLUMES" << std::endl; }

                while ( getline(fileVolumes,building_name) )
                {
                    SpNbTasksToSubmit++;
                    if (QViewInfoSpecs) { std::cout << "[INFO SPECS] : NAME =>"<<building_name<<"\n"; }

                    std::promise<int> promise1;

                    runtime_Building.task(SpRead(building_name), [&promise1](const std::string & bn){
                        promise1.get_future().get();
                    }).setTaskName("Start Task");

                    runtime_Building.task(
                        SpPriority(1),
                        SpRead(building_name),
                        [&](const std::string & bn) -> bool {
                        computeMasksOneBuilding(bn);
                        //computeMasksOneBuildingSpecx(bn);
                        //usleep(1000);
                        return true;
                        }
                    ).setTaskName("Task 2V Compute :"+building_name);

                    promise1.set_value(0);
                    if (QViewInfoSpecs) { std::cout << "[INFO SPECS] : OK END"<<"\n"; }
                }
                runtime_Building.waitAllTasks();
                
            }
            else
            {
                while ( getline(fileVolumes,building_name) )
                {
                    //std::cout << "[INFO NO SPECS] : BINGGGGGGGOOOOO =>"<<building_name << std::endl;
                    computeMasksOneBuilding(building_name);
                }
            }



        }
            // read, line by line, the building marker
        else if( j_["/Buildings"_json_pointer].contains("fileSurfaces") ) // a csv containing the surface markers is provided
        {
            if (QViewInfoSpecs) { std::cout << "[INFO SPECS] : COMPUTE 3" << std::endl; }
            std::string building_name;
            std::ifstream fileSurfaces(Environment::expand(j_["Buildings"]["fileSurfaces"].get<std::string>()));

            // read, line by line, the building marker
            while ( getline(fileSurfaces,building_name) )
            {
                if (QViewInfoSpecs) { std::cout << "[INFO NO SPECS] : BINGGGGGGGOOOOO =>"<<building_name << std::endl; }
                computeMasksOneBuilding(building_name);
            }
        }


       
       
        if (QActiveSpecx) {  
            runtime_Building.waitAllTasks();
            runtime_Building.stopAllThreads();
            timer.stop();
            std::cout << "[INFO SPECS] : Task Time = "<<timer.getElapsed()<< " s" << std::endl; 
            std::cout << "[INFO SPECS] : Task Average Time = " << timer.getElapsed()/double(SpNbTasksToSubmit)<<" s"<<std::endl; 
            std::cout << "[INFO SPECS] : Nb Tasks = " << SpNbTasksToSubmit<<std::endl; 
            runtime_Building.generateDot("RuntimeShadingMask2.dot",true);
            runtime_Building.generateTrace("RuntimeShadingMask2.svg");
            std::cout << "[INFO SPECS] : COMPUTE END" << std::endl;
            usleep(1000);
        }

    }



    // Compute shading masks for one building only
    void computeMasksOneBuilding(std::string building_name)//, BVHTree<MeshType::nDim> bvh_tree)
    {

        int dim = M_submeshes[building_name]->realDimension();
        std::vector<double> random_direction(dim);

        //std::cout << "Submeshes markers" << M_submeshes[building_name]->markerNames() << std::endl;
        //std::cout <<"Name "<<building_name<< " Submeshes markers" << M_submeshes[building_name]->markerNames() << std::endl;

        // Loop over the markers of the building
        Eigen::MatrixXd SM_table_marker(M_azimuthSize,M_altitudeSize);
        Eigen::MatrixXd Angle_table_marker(M_azimuthSize,M_altitudeSize);

        for(auto  [marker,marker_id] : M_submeshes[building_name]->markerNames())
        {

            SM_table_marker.setZero();
            Angle_table_marker.setZero();
            auto ray_submesh = createSubmesh(_mesh=M_submeshes[building_name],_range=markedelements(M_submeshes[building_name],marker));

            // Launch Nrays from each triangle of each marker
            for(auto const &el : ray_submesh->elements() ) // from each element of the submesh, launch M_Nrays randomly oriented
            {
                    auto rays_from_element = [&,marker=marker](int n_rays_thread){

                        Eigen::MatrixXd SM_table(M_azimuthSize,M_altitudeSize);
                        SM_table.setZero();

                        Eigen::MatrixXd Angle_table(M_azimuthSize,M_altitudeSize);
                        Angle_table.setZero();

                        int index_altitude;
                        int index_azimuth;
                        for(int i=0;i<n_rays_thread;i++)
                        {

                            // Construct the ray emitting from a random point of the element
                            auto random_origin = get_random_point(el.second.vertices());

                            Eigen::VectorXd rand_dir(dim);
                            Eigen::VectorXd p1(dim),p2(dim),p3(dim),origin(3);
                            bool inward_ray=false;
                            if(dim==3)
                            {
                                for(int i=0;i<dim;i++)
                                {
                                    p1(i)=column(el.second.vertices(), 0)[i];
                                    p2(i)=column(el.second.vertices(), 1)[i];
                                    p3(i)=column(el.second.vertices(), 2)[i];
                                    origin(i) = random_origin[i];
                                    rand_dir(i) = random_direction[i];
                                }
                                auto element_normal = ((p3-p1).head<3>()).cross((p2-p1).head<3>());
                                element_normal.normalize();

                                // Choose the direction randomly among the latitude and azimuth
                                getRandomDirectionSM(random_direction,M_gen,M_gen2,index_azimuth,index_altitude);
                                for(int i=0;i<dim;i++)
                                {
                                    rand_dir(i) = random_direction[i];
                                }
                                if(rand_dir.dot(element_normal)>=0)
                                {
                                    inward_ray=true;
                                }

                            }

                            BVHRay<mesh_type::nRealDim> ray( origin, rand_dir, 1e-8 );

                            int closer_intersection_element = -1;
                            if(inward_ray)
                            {
                                closer_intersection_element = 1;
                            }
                            else
                            {
                                for(auto& [building_name,bvh_building_tree] : M_bvh_tree_vector)
                                {
                                    auto rayIntersectionResult =  bvh_building_tree->intersect(ray) ;
                                    if ( !rayIntersectionResult.empty() )
                                        closer_intersection_element = 1;
                                    if (closer_intersection_element >=0 )
                                        break;
                                }
                            }
                            // If there is an intersection, increase the shading mask table entry by 1 and augment the angle table by 1 as well
                            if ( closer_intersection_element >=0 )
                            {
                                SM_table(index_azimuth,index_altitude)++;
                                Angle_table(index_azimuth,index_altitude)++;
                            }
                            else
                            {
                                Angle_table(index_azimuth,index_altitude)++;
                            }
                        }

                        return std::make_pair(SM_table,Angle_table);
                    };

                // Execute the lambda function on multiple threads using
                // std::async and std::future to collect the results
                std::vector<int> n_rays_thread;
                n_rays_thread.push_back(M_Nrays - (M_Nthreads-1) * (int)(M_Nrays / M_Nthreads));
                for(int t= 1; t < M_Nthreads; ++t){
                   n_rays_thread.push_back( M_Nrays / M_Nthreads);
                }

                // Used to store the future results
                std::vector< std::future< std::pair<Eigen::MatrixXd, Eigen::MatrixXd > > > futures;

                for(int t = 0; t < M_Nthreads; ++t){

                    // Start a new asynchronous task
                    futures.emplace_back(std::async(std::launch::async, rays_from_element, n_rays_thread[t]));
                }

                for( auto& f : futures){
                    // Wait for the result to be ready
                    auto two_tables =  f.get();

                    // Add the tables obtained in threads
                    SM_table_marker +=two_tables.first;
                    Angle_table_marker += two_tables.second;

                }
            }
            // Divide the shading mask by the corresponding value of the angle table
            // If an angle combination has not been selected, suppose there is no shadow
            auto shadingMatrix = SM_table_marker.array().binaryExpr( Angle_table_marker.array() , [](auto x, auto y) { return y==0 ? 0 : x/y; });

            // Shading mask value 0 means that the surface is not shadowed, value 1 it is fully shadowed
            // Save the shading mask table to a csv file
            saveShadingMask(building_name,marker,shadingMatrix.matrix());

        }
    }


    // Compute shading masks for one building only
    void computeMasksOneBuildingSpecx(std::string building_name)//, BVHTree<MeshType::nDim> bvh_tree)
    {

            int dim = M_submeshes[building_name]->realDimension();
        std::vector<double> random_direction(dim);

        //std::cout << "Submeshes markers" << M_submeshes[building_name]->markerNames() << std::endl;
        //std::cout <<"Name "<<building_name<< " Submeshes markers" << M_submeshes[building_name]->markerNames() << std::endl;

        // Loop over the markers of the building
        Eigen::MatrixXd SM_table_marker(M_azimuthSize,M_altitudeSize);
        Eigen::MatrixXd Angle_table_marker(M_azimuthSize,M_altitudeSize);

        SpRuntime runtime(100);

        for(auto  [marker,marker_id] : M_submeshes[building_name]->markerNames())
        {

            SM_table_marker.setZero();
            Angle_table_marker.setZero();
            auto ray_submesh = createSubmesh(_mesh=M_submeshes[building_name],_range=markedelements(M_submeshes[building_name],marker));

            // Launch Nrays from each triangle of each marker
            for(auto const &el : ray_submesh->elements() ) // from each element of the submesh, launch M_Nrays randomly oriented
            {
                    auto rays_from_element = [&,marker=marker](int n_rays_thread){

                        Eigen::MatrixXd SM_table(M_azimuthSize,M_altitudeSize);
                        SM_table.setZero();

                        Eigen::MatrixXd Angle_table(M_azimuthSize,M_altitudeSize);
                        Angle_table.setZero();

                        int index_altitude;
                        int index_azimuth;
                        for(int i=0;i<n_rays_thread;i++)
                        {

                            int val;
                            runtime.task(
                            SpRead(val),
                            [&](const int &) {

                                // Construct the ray emitting from a random point of the element
                                auto random_origin = get_random_point(el.second.vertices());

                                Eigen::VectorXd rand_dir(dim);
                                Eigen::VectorXd p1(dim),p2(dim),p3(dim),origin(3);
                                bool inward_ray=false;
                                if(dim==3)
                                {
                                    for(int i=0;i<dim;i++)
                                    {
                                        p1(i)=column(el.second.vertices(), 0)[i];
                                        p2(i)=column(el.second.vertices(), 1)[i];
                                        p3(i)=column(el.second.vertices(), 2)[i];
                                        origin(i) = random_origin[i];
                                        rand_dir(i) = random_direction[i];
                                    }
                                    auto element_normal = ((p3-p1).head<3>()).cross((p2-p1).head<3>());
                                    element_normal.normalize();

                                    // Choose the direction randomly among the latitude and azimuth
                                    getRandomDirectionSM(random_direction,M_gen,M_gen2,index_azimuth,index_altitude);
                                    for(int i=0;i<dim;i++)
                                    {
                                        rand_dir(i) = random_direction[i];
                                    }
                                    if(rand_dir.dot(element_normal)>=0)
                                    {
                                        inward_ray=true;
                                    }

                                }

                                BVHRay<mesh_type::nRealDim> ray( origin, rand_dir, 1e-8 );

                                int closer_intersection_element = -1;
                                if(inward_ray)
                                {
                                    closer_intersection_element = 1;
                                }
                                else
                                {
                                    for(auto& [building_name,bvh_building_tree] : M_bvh_tree_vector)
                                    {
                                        auto rayIntersectionResult =  bvh_building_tree->intersect(ray) ;
                                        if ( !rayIntersectionResult.empty() )
                                            closer_intersection_element = 1;
                                        if (closer_intersection_element >=0 )
                                            break;
                                    }
                                }
                                // If there is an intersection, increase the shading mask table entry by 1 and augment the angle table by 1 as well
                                if ( closer_intersection_element >=0 )
                                {
                                    SM_table(index_azimuth,index_altitude)++;
                                    Angle_table(index_azimuth,index_altitude)++;
                                }
                                else
                                {
                                    Angle_table(index_azimuth,index_altitude)++;
                                }

                            } 
                            ).setTaskName("Task");
                        }//End for

                        runtime.waitAllTasks();
                        //runtime.stopAllThreads();
                        //timer.stop();
           
                        //runtime.generateDot("RuntimeRaysThread.dot",true);
                        //runtime.generateTrace("RuntimeRaysThread.svg");




                        return std::make_pair(SM_table,Angle_table);
                    };

                // Execute the lambda function on multiple threads using
                // std::async and std::future to collect the results
                std::vector<int> n_rays_thread;
                n_rays_thread.push_back(M_Nrays - (M_Nthreads-1) * (int)(M_Nrays / M_Nthreads));
                for(int t= 1; t < M_Nthreads; ++t){
                   n_rays_thread.push_back( M_Nrays / M_Nthreads);
                }

                // Used to store the future results
                std::vector< std::future< std::pair<Eigen::MatrixXd, Eigen::MatrixXd > > > futures;

                for(int t = 0; t < M_Nthreads; ++t){

                    // Start a new asynchronous task
                    futures.emplace_back(std::async(std::launch::async, rays_from_element, n_rays_thread[t]));
                }

                for( auto& f : futures){
                    // Wait for the result to be ready
                    auto two_tables =  f.get();

                    // Add the tables obtained in threads
                    SM_table_marker +=two_tables.first;
                    Angle_table_marker += two_tables.second;

                }
            }
            // Divide the shading mask by the corresponding value of the angle table
            // If an angle combination has not been selected, suppose there is no shadow
            auto shadingMatrix = SM_table_marker.array().binaryExpr( Angle_table_marker.array() , [](auto x, auto y) { return y==0 ? 0 : x/y; });

            // Shading mask value 0 means that the surface is not shadowed, value 1 it is fully shadowed
            // Save the shading mask table to a csv file
            saveShadingMask(building_name,marker,shadingMatrix.matrix());

        }   
    }

    // Save the shading mask table to a CSV file
    void saveShadingMask(std::string building_name, std::string marker_name, const Eigen::Ref<const Eigen::MatrixXd>& M)
    {
        int i_long=0;
        int j_alt=0;

        Eigen::MatrixXd matrix_sm(M_azimuthSize,M_altitudeSize);

        // Save the matrix into a CSV file, inside the shadingMasks subfolder of the results folder
        std::ofstream matrix_file;
        std::string shadingMaskFolder = (boost::filesystem::path(Environment::appRepository())/("shadingMasks")).string();
        if (!boost::filesystem::exists(shadingMaskFolder))
            boost::filesystem::create_directory(shadingMaskFolder);

        std::string matrix_filename = shadingMaskFolder+"/SM_Matrix_"+building_name+"_"+marker_name+".csv";
        matrix_file.open(matrix_filename,std::ios_base::out);
        for(int i=0; i<M_azimuthSize; i++)
        {
            for(int j=0; j<M_altitudeSize-1; j++)
            {
                matrix_file << M(i,j) << ",";
            }
            matrix_file << M(i,M_altitudeSize-1) << "\n";
        }
        matrix_file.close();
    }


    //std::map<std::string,BVHTree<MeshType::nDim>> M_bvh_tree_vector;
    std::map<std::string,std::unique_ptr<BVH<typename tr_mesh_type::element_type>>> M_bvh_tree_vector;
    std::map<std::string,tr_mesh_ptrtype> M_submeshes;
    std::map<int,node_type> M_faces_to_normals;

    Eigen::VectorXd M_azimuthAngles, M_altitudeAngles;

    int M_azimuthSize;
    int M_altitudeSize;
    int M_Nrays;
    int M_Nthreads;

    nl::json j_;

    std::random_device M_rd;
    std::random_device M_rd2;
    std::mt19937 M_gen;
    std::mt19937 M_gen2;
};
} // namespace Feel
