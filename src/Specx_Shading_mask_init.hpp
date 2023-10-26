
namespace Feel
{
template <typename MeshType>
ShadingMask<MeshType>::ShadingMask(mesh_ptrtype mesh, nl::json const& specs, int intervalsAzimuth, int intervalsAltitude )
{
    auto start_computation = std::chrono::system_clock::now();
    std::time_t beginning_time = std::chrono::system_clock::to_time_t(start_computation);
    M_metadataJson["shadingMask"]["Timestamp"]["Beginning"] = strtok(std::ctime(&beginning_time),"\n");;

    // Read the number of rays per triangle and the number of threads
    j_ = specs;
    M_Nrays = specs["Nrays"];
    M_Nthreads = specs["Multithread"]["Nthreads"].get<int>();
    M_mthreadtype = specs["Multithread"]["type"].get<std::string>();
    M_saveMasks = specs["SaveMasks"];


    SpecxNbThreadDesired     = SpUtils::DefaultNumThreads();
    SpecxNbThreadDesired     = 5;
    SpecxNbThreadDesired     = M_Nthreads;
    SpecxSaveNbThreadDesired = SpUtils::DefaultNumThreads();;
    SpecxSaveNbThreadDesired = 80;
    QSaveSpecxGenerateReports= true;
    QSpecxLockConfigWaitSave = true;
    //QSpecxLockConfigWaitSave=false;
    QSaveWithSpecx           = false;
    QViewInfoSpecx           = true;

    // Fix the size of the shading mask matrix
    fixAzimuthAltitudeDiscretization(intervalsAzimuth, intervalsAltitude);

    // Create the random number generators
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

    // Create and store the directions of the M_Nrays
    int index_azimuth, index_altitude;
    M_raysdirections.resize(M_Nrays);
    std::vector<double> random_direction(3);
    std::map< std::pair<int,int>, std::vector<int> > check_directions;
    for(int i=0; i<M_Nrays; i++)
    {
        getRandomDirectionSM(random_direction,M_gen,M_gen2,index_azimuth,index_altitude);                    
        M_raysdirections[i] = std::make_tuple(random_direction,index_azimuth,index_altitude);
        check_directions[std::make_pair(index_azimuth,index_altitude)].push_back(i);
    }
    // Check if all the possible combinations of [0,intervalsAzimuth] x [0,intervalsAltitude] have at least one associated ray 
        std::map< std::pair<int,int>, std::vector<int> >::iterator it;
    for(int i=0; i<intervalsAzimuth; i++)
    {
        for(int j=0; j<intervalsAltitude; j++)
        {
            it = check_directions.find(std::make_pair(i,j));
            if(it == check_directions.end())
            {
                std::cout << fmt::format("Direction associated with indices ({},{}) is missing. Replacing one direction with it \n",i,j);

                std::pair<int,int> kl_pair;
                bool leave_loop = false;

                for(int k=0;k<intervalsAzimuth; k++)
                {
                    for(int l=0; l<intervalsAltitude; l++)
                    {
                        kl_pair=std::make_pair(k,l);
                        if( check_directions[kl_pair].size() > 1 )
                        {
                            std::vector<int> list_indices = check_directions[kl_pair];
                            int index_to_substitute = list_indices.back();
                            check_directions[kl_pair].pop_back();
                            std::cout << fmt::format("Inserting direction associated with indices ({},{}) and deleting one direction associated with indices ({},{}) \n",i,j,k,l);

                            // Compute the direction associated with the indices (i,j)
                            double phi = -( M_azimuthAngles[i] ) + M_PI*0.5 ; // recover spherical coordinate from azimuth angle
                            double theta = M_PI*0.5 - M_altitudeAngles[j]; // recover spherical coordinate from altitude

                            random_direction[0]=math::sin(theta)*math::cos(phi);
                            random_direction[1]=math::sin(theta)*math::sin(phi);
                            random_direction[2]=math::cos(theta);

                            M_raysdirections[index_to_substitute] = std::make_tuple(random_direction,i,j);
                            check_directions[std::make_pair(i,j)].push_back(index_to_substitute);

                            leave_loop = true;
                        }
                        if(leave_loop)
                            break;
                    }
                    if(leave_loop)
                            break;
                }
            }
        }
    }
    
    if constexpr( MeshType::nDim==MeshType::nRealDim )
    {
        // For each building, save the surface mesh and build the corresponding BVH tree for ray search
        if( specs["/Buildings"_json_pointer].contains("list") ) // the list of volume markers is provided
        {

            std::cout<<"[SPECX INFO] : *Section List\n";

            auto markersVolume = specs["Buildings"]["list"].get<std::vector<std::string>>();
            int nBuildings = 0;
            int nFaces = 0;

            NbObjects=markersVolume.size();
            int NbLoop=0; int NbTh=0; int NbCoor=0;
            GetSpecxPreprocessingParameters(NbObjects,SpecxNbThreadDesired,NbLoop,NbTh,NbCoor);
            std::cout<<"[SPECX INFO] : SpecxNbThreadDesired="<<SpecxNbThreadDesired<<"\n";
            std::cout<<"[SPECX INFO] : Nb Objects="<<NbObjects<<"\n";
            std::cout<<"[SPECX INFO] : Nb Loops="<<NbLoop<<" Th="<<NbTh<<" Coor="<<NbCoor<<"\n";
            std::cout<<"[SPECX INFO] : NbThreads="<<NbTh<<" used\n";
            SpRuntime runtime0(NbTh);
            tic();
            std::string buildingName;
            std::atomic<int> Index(0);

            int NbidxN=NbTh;
                for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                    if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }
                    for(int idx = 0 ; idx < NbTh ; ++idx){
                            auto returnValue=runtime0.task(SpWrite(Index),
                                [&](std::atomic<int>&  IndexA)-> bool {
                                    std::cout <<"Index:"<<IndexA<<" ";
                                    buildingName=markersVolume[Index];
                                    std::cout << fmt::format("{}\t",buildingName);
                                    auto volumeSubmesh  = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                                    auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));
                                    auto bvhBuilding    = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));

                                    M_bvh_tree_vector.insert(std::make_pair( buildingName, std::move(bvhBuilding) ));
                                    M_submeshes.insert(std::make_pair( buildingName, surfaceSubmesh ));

                                    nBuildings +=1;
                                    nFaces += nelements(elements(surfaceSubmesh));
                                    usleep(10);
                                    IndexA++;
                                    std::cout <<"\n";

                                    return(true);
                                }
                            );                         
                           

                    }//End for 

                    runtime0.waitAllTasks();
                    

                }//End For idxLoop
            

            /*
            for(std::string buildingName : markersVolume)
            {
                std::cout << fmt::format("{}\n",buildingName);
                auto volumeSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));

                auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                M_bvh_tree_vector.insert(std::make_pair( buildingName , std::move(bvhBuilding) ));

                M_submeshes.insert(std::make_pair( buildingName , surfaceSubmesh ));

                nBuildings +=1;
                nFaces += nelements(elements(surfaceSubmesh));
            }
            */
        


           
            auto bvhBuildingTime = toc("BVH built");


            runtime0.waitAllTasks();
            runtime0.stopAllThreads();
            std::cout<<"[SPECX INFO] : STOP ALL SPECX THREAD\n";
            runtime0.generateDot("Runtime0.dot",true);
            runtime0.generateTrace("Runtime0.svg");



            LOG(INFO) << "BVHs construction: end";

            M_metadataJson["shadingMask"]["Timer"]["BVHs_total_building_time"] = bvhBuildingTime;

            M_metadataJson["shadingMask"]["Method"] = "listVolumes";
            M_metadataJson["shadingMask"]["nBuildings"] = nBuildings;
            M_metadataJson["shadingMask"]["nFaces"] = nFaces;
        }

        else if ((true) && ( specs["/Buildings"_json_pointer].contains("fileVolumes"))) // a csv containing the volume markers is provided
        {
            std::cout<<"[SPECX INFO] : *Section FileVolumes\n";

            std::string buildingName;
            // open the file

            std::ifstream fileVolumes(Environment::expand(specs["Buildings"]["fileVolumes"].get<std::string>()));

            int nBuildings = 0;
            int nFaces = 0;

            std::cout<<"[SPECX INFO] : List\n";
            std::vector<std::string> ListObjects;
            ListObjects=GetListNameObjects(Environment::expand(specs["Buildings"]["fileVolumes"].get<std::string>()));
            NbObjects=ListObjects.size();
            //for(int k = 0 ; k < NbObjects ; ++k){ std::cout<<ListObjects[k]<<"\n";} 
            int NbLoop=0; int NbTh=0; int NbCoor=0; 
            GetSpecxPreprocessingParameters(NbObjects,SpecxNbThreadDesired,NbLoop,NbTh,NbCoor);
            std::cout<<"[SPECX INFO] : SpecxNbThreadDesired="<<SpecxNbThreadDesired<<"\n";
            std::cout<<"[SPECX INFO] : Nb Objects="<<NbObjects<<"\n";
            std::cout<<"[SPECX INFO] : Nb Loops="<<NbLoop<<" Th="<<NbTh<<" Coor="<<NbCoor<<"\n";
            std::cout<<"[SPECX INFO] : NbThreads="<<NbTh<<" used\n";
            SpRuntime runtime0(NbTh);
            //SpRuntime<SpSpeculativeModel::SP_MODEL_2> runtime0;
            
            tic();

            int Index=0;

            //std::atomic<int> Index(0);



            int NumModel=1;
           

            
            if (NumModel==2) {
                int NbidxN=NbTh;
                for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                    if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }

                    std::promise<int> promise0;

                    runtime0.task(SpWrite(Index), [&promise0](int& valParam){
                                promise0.get_future().get();
                            }).setTaskName("ST");

                    for(int idx = 0 ; idx < NbTh ; ++idx){
                            auto returnValue=runtime0.task(SpPotentialWrite(Index),
                                [&](int& IndexA)-> bool {
                                    std::cout <<"Index:"<<IndexA<<" ";
                                    buildingName=ListObjects[IndexA];
                                    std::cout << fmt::format("{}\t",buildingName);
                                    auto volumeSubmesh  = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                                    auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));
                                    auto bvhBuilding    = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));

                                    M_bvh_tree_vector.insert(std::make_pair( buildingName, std::move(bvhBuilding) ));
                                    M_submeshes.insert(std::make_pair( buildingName, surfaceSubmesh ));

                                    nBuildings +=1;
                                    nFaces += nelements(elements(surfaceSubmesh));
                                    usleep(10);
                                    IndexA++;
                                    std::cout <<"\n";
                                    return(true);
                                }
                            );  
                    }//End for 
                    promise0.set_value(0);
                    runtime0.waitAllTasks();
                }//End For idxLoop
            } //END MODEL 2


            if (NumModel==3) {
                int NbidxN=NbTh;


                std::promise<int> promise0;
                    runtime0.task(SpWrite(Index), [&promise0](int& valParam){
                                promise0.get_future().get();
                            }).setTaskName("ST");

                for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                    if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }
                    for(int idx = 0 ; idx < NbTh ; ++idx){
                            auto returnValue=runtime0.task(SpPotentialWrite(Index),
                                [&](int& IndexA)-> bool {
                                    std::cout <<"Index:"<<IndexA<<" ";
                                    buildingName=ListObjects[IndexA];
                                    std::cout << fmt::format("{}\t",buildingName);
                                    auto volumeSubmesh  = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                                    auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));
                                    auto bvhBuilding    = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));

                                    M_bvh_tree_vector.insert(std::make_pair( buildingName, std::move(bvhBuilding) ));
                                    M_submeshes.insert(std::make_pair( buildingName, surfaceSubmesh ));

                                    nBuildings +=1;
                                    nFaces += nelements(elements(surfaceSubmesh));
                                    usleep(10);
                                    IndexA++;
                                    std::cout <<"\n";
                                    return(true);
                                }
                            );  
                    }//End for 
                }//End For idxLoop
                promise0.set_value(0);
                runtime0.waitAllTasks();
            } //END MODEL 3
        

            

            
            if (NumModel==1) {
                int NbidxN=NbTh;
                for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                    if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }
                    for(int idx = 0 ; idx < NbTh ; ++idx){
                            auto returnValue=runtime0.task(SpWrite(Index),
                                //[&](std::atomic<int>&  IndexA)-> bool {
                                [&](int& IndexA)-> bool {
                                    std::cout <<"Index:"<<IndexA<<" ";
                                    buildingName=ListObjects[IndexA];
                                    std::cout << fmt::format("{}\t",buildingName);
                                    auto volumeSubmesh  = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                                    auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));
                                    auto bvhBuilding    = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));

                                    M_bvh_tree_vector.insert(std::make_pair( buildingName, std::move(bvhBuilding) ));
                                    M_submeshes.insert(std::make_pair( buildingName, surfaceSubmesh ));

                                    nBuildings +=1;
                                    nFaces += nelements(elements(surfaceSubmesh));
                                    usleep(10);
                                    IndexA++;
                                    std::cout <<"\n";

                                    return(true);
                                }
                            );                         
                           

                    }//End for 
                    runtime0.waitAllTasks();
                }//End For idxLoop
            } //END MODEL 1

            

            nBuildings=NbObjects;


            // read, line by line, the building marker
            /*
            while ( getline(fileVolumes,buildingName) )
            {
                std::cout << fmt::format("{}\n",buildingName);
                auto volumeSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));
                auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                M_bvh_tree_vector.insert(std::make_pair( buildingName , std::move(bvhBuilding) ));

                M_submeshes.insert(std::make_pair( buildingName , surfaceSubmesh ));

                nBuildings +=1;
                nFaces += nelements(elements(surfaceSubmesh));

            }
            */



            std::cout <<"\n";
            auto bvhBuildingTime = toc("BVHs built");

            
            runtime0.waitAllTasks();
            runtime0.stopAllThreads();
            std::cout<<"[SPECX INFO] : STOP ALL SPECX THREAD\n";
            runtime0.generateDot("Runtime0.dot",true);
            runtime0.generateTrace("Runtime0.svg");
            

            LOG(INFO) << "BVHs construction: end";

            M_metadataJson["shadingMask"]["Timer"]["BVHs_total_building_time"] = bvhBuildingTime;

            M_metadataJson["shadingMask"]["Method"] = "fileVolumes";
            M_metadataJson["shadingMask"]["nBuildings"] = nBuildings;
            M_metadataJson["shadingMask"]["nFaces"] = nFaces;

        }
        else if ((false) && ( specs["/Buildings"_json_pointer].contains("fileVolumes"))) // a csv containing the volume markers is provided
        {
            std::cout<<"[SPECX INFO] : *Section FileVolumes\n";

            std::string buildingName;
            // open the file

            std::ifstream fileVolumes(Environment::expand(specs["Buildings"]["fileVolumes"].get<std::string>()));

            int nBuildings = 0;
            int nFaces = 0;

            std::cout<<"[SPECX INFO] : List\n";
            std::vector<std::string> ListObjects;
            ListObjects=GetListNameObjects(Environment::expand(specs["Buildings"]["fileVolumes"].get<std::string>()));
            NbObjects=ListObjects.size();
            //for(int k = 0 ; k < NbObjects ; ++k){ std::cout<<ListObjects[k]<<"\n";} 
            int NbLoop=0; int NbTh=0; int NbCoor=0; 
            GetSpecxPreprocessingParameters(NbObjects,SpecxNbThreadDesired,NbLoop,NbTh,NbCoor);
            std::cout<<"[SPECX INFO] : SpecxNbThreadDesired="<<SpecxNbThreadDesired<<"\n";
            std::cout<<"[SPECX INFO] : Nb Objects="<<NbObjects<<"\n";
            std::cout<<"[SPECX INFO] : Nb Loops="<<NbLoop<<" Th="<<NbTh<<" Coor="<<NbCoor<<"\n";
            std::cout<<"[SPECX INFO] : NbThreads="<<NbTh<<" used\n";

            
            tic();
            int Index=0;
            int NumModel=2;    

            if (NumModel==2) {
                SpConsumerThread cs;
                int NbidxN=NbTh;
                for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                    if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }
                    for(int idx = 0 ; idx < NbTh ; ++idx){
                            cs.submitJob([&](){
                                std::cout <<"Index:"<<Index<<" ";
                                buildingName=ListObjects[Index];
                                std::cout << fmt::format("{}\t",buildingName);
                                auto volumeSubmesh  = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                                auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));
                                auto bvhBuilding    = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));

                                M_bvh_tree_vector.insert(std::make_pair( buildingName, std::move(bvhBuilding) ));
                                M_submeshes.insert(std::make_pair( buildingName, surfaceSubmesh ));

                                nBuildings +=1;
                                nFaces += nelements(elements(surfaceSubmesh));
                                usleep(10);
                                Index++;
                                std::cout <<"\n";
                                }
                            );                         
                           

                    }//End for 

                    cs.submitJobAndWait([&](){
               
                    });
                    
                }//End For idxLoop
                nBuildings=NbObjects;
            }

            

            if (NumModel==1) {
                while ( getline(fileVolumes,buildingName) )
                {
                    std::cout << fmt::format("{}\n",buildingName);
                    auto volumeSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName));
                    auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh));
                    auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                    M_bvh_tree_vector.insert(std::make_pair( buildingName , std::move(bvhBuilding) ));

                    M_submeshes.insert(std::make_pair( buildingName , surfaceSubmesh ));

                    nBuildings +=1;
                    nFaces += nelements(elements(surfaceSubmesh));

                }
            }

            std::cout <<"\n";
            auto bvhBuildingTime = toc("BVHs built");

            LOG(INFO) << "BVHs construction: end";

            M_metadataJson["shadingMask"]["Timer"]["BVHs_total_building_time"] = bvhBuildingTime;

            M_metadataJson["shadingMask"]["Method"] = "fileVolumes";
            M_metadataJson["shadingMask"]["nBuildings"] = nBuildings;
            M_metadataJson["shadingMask"]["nFaces"] = nFaces;

        }
    }
    else
    {
        // For each building, save the surface mesh and build the corresponding BVH tree for ray search
        if( specs["/Buildings"_json_pointer].contains("fileSurfaces") ) // a csv containing the surface markers is provided
        {
            std::cout<<"[SPECX INFO] : *Section FileSurfaces\n";


            std::string buildingName;
            std::ifstream fileSurfaces(Environment::expand(specs["Buildings"]["fileSurfaces"].get<std::string>()));
            //std::cout << Environment::expand(specs["Buildings"]["fileSurfaces"].get<std::string>()) << std::endl;

            int nBuildings = 0;
            int nFaces = 0;


            std::cout<<"[SPECX INFO] : List\n";
            std::vector<std::string> ListObjects;
            ListObjects=GetListNameObjects(Environment::expand(specs["Buildings"]["fileSurfaces"].get<std::string>()));
            NbObjects=ListObjects.size();
            //for(int k = 0 ; k < NbObjects ; ++k){ std::cout<<ListObjects[k]<<"\n";} //CTRL
            
            std::cout<<"[SPECX INFO] : Nb Objects="<<NbObjects<<"\n";
            std::cout<<"[SPECX INFO] : NbThreads Desired="<<M_Nthreads<<"\n";
                
            int NbLoop=0; int NbTh=0; int NbCoor=0;
            GetSpecxPreprocessingParameters(NbObjects,SpecxNbThreadDesired,NbLoop,NbTh,NbCoor);
            std::cout<<"[SPECX INFO] : Nb Loops="<<NbLoop<<" Th="<<NbTh<<" Coor="<<NbCoor<<"\n";

            SpRuntime runtime0(NbTh);

            tic();
            
            std::cout<<"[SPECX INFO] : >>>";
            int initVal = 0;


            int NbidxN=NbTh;
                for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                    if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }
                    for(int idx = 0 ; idx < NbTh ; ++idx){
                            int Index=idxLoop*NbidxN+idx;
                            //std::cout <<"Index:"<<Index<<"\n";
                            buildingName=ListObjects[Index];
                            
                            auto returnValue=runtime0.task(SpRead(buildingName),
                                [&](const std::string & buildingNameS)-> bool {
                                    std::cout << fmt::format("{}\t",buildingNameS);
                                    
                                    auto surfaceSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingNameS));
                                    auto listMarkers = surfaceSubmesh->markerNames();
               
                                    auto it = listMarkers.find(buildingNameS);
                                    listMarkers.erase(it);
                                    surfaceSubmesh->setMarkerNames(listMarkers);

                                    auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));


                                    M_bvh_tree_vector.insert(std::make_pair( buildingNameS , std::move(bvhBuilding) ));

                                    M_submeshes.insert(std::make_pair( buildingNameS , surfaceSubmesh ));

                                    nBuildings +=1;
                                    nFaces += nelements(elements(surfaceSubmesh));
                                   
                                    return(true);
                                }
                            );
                                                        
                            returnValue.wait();
                            
                            

                    }//End for 
                    runtime0.waitAllTasks();

                }//End For idxLoop
  
                nBuildings=NbObjects;



            
            // read, line by line, the building marker

            /*
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

                nBuildings +=1;
                nFaces += nelements(elements(surfaceSubmesh));

            }
            */

            std::cout <<"\n";
            auto bvhBuildingTime = toc("BVHs built");
            
            runtime0.waitAllTasks();
            runtime0.stopAllThreads();
            std::cout<<"[SPECX INFO] : STOP ALL SPECX THREAD\n";
            runtime0.generateDot("Runtime0.dot",true);
            runtime0.generateTrace("Runtime0.svg");
            //runtime0.finish(); SpRuntime<SpSpeculativeModel::SP_MODEL_1>

            LOG(INFO) << "BVHs construction: end";

            NbObjects=nBuildings;

            M_metadataJson["shadingMask"]["Timer"]["BVHs_total_building_time"] = bvhBuildingTime;

            M_metadataJson["shadingMask"]["Method"] = "fileSurfaces";
            M_metadataJson["shadingMask"]["nBuildings"] = nBuildings;
            M_metadataJson["shadingMask"]["nFaces"] = nFaces;
        }
        // Store only the view on the surface mesh faces
        else if( specs["/Buildings"_json_pointer].contains("fileFaces") ) // a csv containing the face markers is provided
        {


            
            std::cout<<"[SPECX INFO] : *Section FileFaces\n";
            std::string faceName;
            std::ifstream fileFaces(Environment::expand(specs["Buildings"]["fileFaces"].get<std::string>()));

            int nMarkers = 0;

            while ( getline(fileFaces,faceName) )
            {
                M_listFaceMarkers.push_back(faceName);
                M_listMarkerFaceEntity[faceName];

                nMarkers += 1;
            }
            M_rangeFaces = markedelements(mesh,M_listFaceMarkers);

            int nFaces = 0;
            for( auto const& face :  M_rangeFaces)
            {
                // Create a map connecting face_id element to marker name (which must contain the string "_face_")
                auto f = boost::unwrap_ref( face );
                for( auto m : f.marker() )
                {                        
                    if (mesh->markerName(m).find("_face_") != std::string::npos)
                    {                     
                        M_mapEntityToBuildingFace.insert( std::make_pair( f.id(), mesh->markerName(m) ) );
                        M_listMarkerFaceEntity[mesh->markerName(m)].push_back(std::ref(f));
                    }
                }     

                nFaces += 1;               
            }

            M_metadataJson["shadingMask"]["Method"] = "fileFaces";
            M_metadataJson["shadingMask"]["nMarkers"] = nMarkers;
            M_metadataJson["shadingMask"]["nFaces"] = nFaces;

            // Create a BVH containing all the faces of the buildings
            tic();
            M_bvh = boundingVolumeHierarchy( _range=M_rangeFaces );
            auto bvhBuildingTime = toc("BVH built");
            LOG(INFO) << "BVH construction: end";

            M_metadataJson["shadingMask"]["Timer"]["BVH_building_time"] = bvhBuildingTime;
        }
        else if( specs["/Buildings"_json_pointer].contains("aggregatedMarkers"))
        {

            std::cout<<"[SPECX INFO] : *Section AggregatedMarkers\n";

            
            M_rangeFaces = markedelements(mesh,"building"); // it contains all the faces of all buildings
            tic();
            int nFaces = 0;
            int nMarkers = 0;
            for( auto const& face :  M_rangeFaces)
            {
                auto f = boost::unwrap_ref( face );
                std::vector<std::string> composite_marker; // collects all necessary marker substrings to compose the face marker using buildingId and faceId
                for( auto m : f.marker() )
                {
                    auto markerName = mesh->markerName(m);

                    if(markerName.find("buildingVerticalFace_") != std::string::npos)
                    {
                        auto pos = markerName.find_last_of('_');
                        composite_marker.push_back("_face" + markerName.substr(pos, std::string::npos));
                    }
                    else if(markerName.find("buildingRoof") != std::string::npos)
                    {
                        composite_marker.push_back( "_face_roof");
                    }
                    else if((markerName.find("buildingId_") != std::string::npos))
                    {
                        auto pos = markerName.find_last_of('_');
                        // insert the building name at the beginning of the vector
                        composite_marker.insert(composite_marker.begin(), "building" + markerName.substr(pos+1));
                    }
                    else
                    {
                        // markers "building" and "terrain" are not useful
                    }
                }

                std::string faceName;
                for(auto marker_ : composite_marker )
                    faceName += marker_;


                if( M_listMarkerFaceEntity[faceName].empty() )
                {
                    M_listFaceMarkers.push_back(faceName);
                    nMarkers +=1;
                }
                M_listMarkerFaceEntity[faceName].push_back(std::ref(f));
                M_mapEntityToBuildingFace.insert( std::make_pair( f.id(), faceName ) );

                nFaces += 1;

            }



            auto dataStructureBuildingTime = toc("Building markers and associated data structures");
            M_metadataJson["shadingMask"]["Method"] = "aggregatedMarkers";
            M_metadataJson["shadingMask"]["nBuildingFaces"] = nFaces;
            M_metadataJson["shadingMask"]["nMarkers"] = nMarkers;
            M_metadataJson["shadingMask"]["Timer"]["DataStructures_building_time"] = dataStructureBuildingTime;

            // Create a BVH containing all the faces of the buildings
            LOG(INFO) << "BVH construction: beginning";
            tic();


            M_bvh = boundingVolumeHierarchy( _range=M_rangeFaces );
            auto bvhBuildingTime = toc("BVH built");
            LOG(INFO) << "BVH construction: end";

            M_metadataJson["shadingMask"]["Timer"]["BVH_building_time"] = bvhBuildingTime;
        }
    }        
}


template <typename MeshType>
void
ShadingMask<MeshType>::fixAzimuthAltitudeDiscretization(int intervalsAzimuth, int intervalsAltitude)
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

}