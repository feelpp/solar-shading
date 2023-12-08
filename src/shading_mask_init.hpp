
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
            auto markersVolume = specs["Buildings"]["list"].get<std::vector<std::string>>();
            int nBuildings = 0;
            int nFaces = 0;

            tic();
            for(std::string buildingName : markersVolume)
            {
                std::cout << fmt::format("{}\n",buildingName);
                auto volumeSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName),_update=0,_view=1);
                auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh),_update=0,_view=1);
                auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                M_bvh_tree_vector.insert(std::make_pair( buildingName , std::move(bvhBuilding) ));

                M_submeshes.insert(std::make_pair( buildingName , surfaceSubmesh ));

                nBuildings +=1;
                nFaces += nelements(elements(surfaceSubmesh));
            }
            auto bvhBuildingTime = toc("BVH built");
            LOG(INFO) << "BVHs construction: end";

            M_metadataJson["shadingMask"]["Timer"]["BVHs_total_building_time"] = bvhBuildingTime;

            M_metadataJson["shadingMask"]["Method"] = "listVolumes";
            M_metadataJson["shadingMask"]["nBuildings"] = nBuildings;
            M_metadataJson["shadingMask"]["nFaces"] = nFaces;
        }
        else if( specs["/Buildings"_json_pointer].contains("fileVolumes")) // a csv containing the volume markers is provided
        {
            std::string buildingName;
            // open the file

            std::ifstream fileVolumes(Environment::expand(specs["Buildings"]["fileVolumes"].get<std::string>()));

            int nBuildings = 0;
            int nFaces = 0;

            tic();
            // read, line by line, the building marker
            while ( getline(fileVolumes,buildingName) )
            {
                std::cout << fmt::format("{}\n",buildingName);
                auto volumeSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName),_update=0,_view=1);
                auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh),_update=0,_view=1);
                auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
                M_bvh_tree_vector.insert(std::make_pair( buildingName , std::move(bvhBuilding) ));

                M_submeshes.insert(std::make_pair( buildingName , surfaceSubmesh ));

                nBuildings +=1;
                nFaces += nelements(elements(surfaceSubmesh));

            }
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
            std::string buildingName;
            std::ifstream fileSurfaces(Environment::expand(specs["Buildings"]["fileSurfaces"].get<std::string>()));
            std::cout << Environment::expand(specs["Buildings"]["fileSurfaces"].get<std::string>()) << std::endl;

            int nBuildings = 0;
            int nFaces = 0;
            tic();
            // read, line by line, the building marker
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
            auto bvhBuildingTime = toc("BVHs built");
            LOG(INFO) << "BVHs construction: end";

            M_metadataJson["shadingMask"]["Timer"]["BVHs_total_building_time"] = bvhBuildingTime;

            M_metadataJson["shadingMask"]["Method"] = "fileSurfaces";
            M_metadataJson["shadingMask"]["nBuildings"] = nBuildings;
            M_metadataJson["shadingMask"]["nFaces"] = nFaces;
        }
        // Store only the view on the surface mesh faces
        else if( specs["/Buildings"_json_pointer].contains("fileFaces") ) // a csv containing the face markers is provided
        {
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





/*
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

    // Fix the size of the shading mask matrix
    fixAzimuthAltitudeDiscretization(intervalsAzimuth, intervalsAltitude);

    // Create the random number generators
    makeRandomNumberGeneratorsSeed();
    makeCreateM_NraysMatrix(intervalsAzimuth,intervalsAltitude);

    loadMeshDataSubPart2(mesh,1);
    

    if constexpr( MeshType::nDim==MeshType::nRealDim )
    {
       
    }
    else
    {
        // For each building, save the surface mesh and build the corresponding BVH tree for ray search
        if( specs["/Buildings"_json_pointer].contains("fileSurfaces") ) // a csv containing the surface markers is provided
        {
            std::string buildingName;
            std::ifstream fileSurfaces(Environment::expand(specs["Buildings"]["fileSurfaces"].get<std::string>()));
            std::cout << Environment::expand(specs["Buildings"]["fileSurfaces"].get<std::string>()) << std::endl;

            int nBuildings = 0;
            int nFaces = 0;
            tic();
            // read, line by line, the building marker
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
            auto bvhBuildingTime = toc("BVHs built");
            LOG(INFO) << "BVHs construction: end";

            M_metadataJson["shadingMask"]["Timer"]["BVHs_total_building_time"] = bvhBuildingTime;

            M_metadataJson["shadingMask"]["Method"] = "fileSurfaces";
            M_metadataJson["shadingMask"]["nBuildings"] = nBuildings;
            M_metadataJson["shadingMask"]["nFaces"] = nFaces;
        }
        // Store only the view on the surface mesh faces
        else if( specs["/Buildings"_json_pointer].contains("fileFaces") ) // a csv containing the face markers is provided
        {
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
*/





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

template <typename MeshType>
void
ShadingMask<MeshType>::makeRandomNumberGeneratorsSeed()
    {
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
    }




template <typename MeshType>
void
ShadingMask<MeshType>::makeCreateM_NraysMatrix(int intervalsAzimuth, int intervalsAltitude)
    {

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
    }







template <typename MeshType>
void
ShadingMask<MeshType>::loadMeshDataSubPart1(mesh_ptrtype mesh,int numOp)
    {   
        // [INFO]: refactoring OK for this parts
        // For each building, save the surface mesh and build the corresponding BVH tree for ray search
        tic();
        int nBuildings = 0;
        int nFaces = 0; 
        int nbObjects = 0;
        std::string buildingName;
        std::string nameFile;
        std::vector<std::string> listObjects;
        nl::json const& markersVolume; 

        if(numOp==1) // the list of volume markers is provided
        {
            markersVolume = j_["Buildings"]["list"].get<std::vector<std::string>>(); //To review if it possible to...
            nbObjects=markersVolume.size(); 
        }
        if (numOp==2) // a csv containing the volume markers is provided
        {
            nameFile=Environment::expand(j_["Buildings"]["fileVolumes"].get<std::string>());
            listObjects=GetListNameObjects(nameFile);
            nbObjects=listObjects.size();       
        }
            
        for(int idx = 0 ; idx <nbObjects; ++idx)
        {
            if (numOp==1) { buildingName=markersVolume[idx]; }
            if (numOp==2) { buildingName=listObjects[idx]; }
            std::cout << fmt::format("{}\n",buildingName);
            auto volumeSubmesh  = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName),_update=0,_view=1);
            auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh),_update=0,_view=1);
            auto bvhBuilding = boundingVolumeHierarchy(_range=elements(surfaceSubmesh));
            M_bvh_tree_vector.insert(std::make_pair( buildingName , std::move(bvhBuilding) ));
            M_submeshes.insert(std::make_pair( buildingName , surfaceSubmesh ));
            nBuildings +=1;
            nFaces += nelements(elements(surfaceSubmesh));
        }

        //BEGIN:SAVE META INFO
        auto bvhBuildingTime = toc("BVHs built");
        LOG(INFO) << "BVHs construction: end";
        M_metadataJson["shadingMask"]["Timer"]["BVHs_total_building_time"] = bvhBuildingTime;
        M_metadataJson["shadingMask"]["Method"] = "fileVolumes";
        M_metadataJson["shadingMask"]["nBuildings"] = nBuildings;
        M_metadataJson["shadingMask"]["nFaces"] = nFaces;
        //END:SAVE META INFO
    }




template <typename MeshType>
void
ShadingMask<MeshType>::loadMeshDataSubPart2(mesh_ptrtype mesh,int numOp)
    {
        // [INFO]: refactoring OK for this parts
        // For each building, save the surface mesh and build the corresponding BVH tree for ray search
        tic();
        int nBuildings = 0;
        int nFaces = 0;
        int nMarkers = 0;
        int nbObjects = 0;
        std::string buildingName;
        std::string faceName;
        std::vector<std::string> listObjects;
        std::string nameFile;

        
        if (numOp==1) { nameFile=Environment::expand(j_["Buildings"]["fileSurfaces"].get<std::string>()); } // a csv containing the surface markers is provided
        if (numOp==2) { nameFile=Environment::expand(j_["Buildings"]["fileFaces"].get<std::string>());} // a csv containing the face markers is provided
        listObjects=GetListNameObjects(nameFile); nbObjects=listObjects.size();

        if (numOp==1) //fileSurfaces
        { 
            // read, line by line, the building marker
            for(int idx = 0 ; idx <nbObjects; ++idx)
            {
                buildingName=listObjects[idx];
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
        }

        if (numOp==2) //fileFaces // Store only the view on the surface mesh faces
        {
            std::ifstream fileFaces(Environment::expand(j_["Buildings"]["fileFaces"].get<std::string>()));
            while ( getline(fileFaces,faceName) )
            {
                M_listFaceMarkers.push_back(faceName);
                M_listMarkerFaceEntity[faceName];
                nMarkers += 1;
            }
            M_rangeFaces = markedelements(mesh,M_listFaceMarkers);
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
        }

    

        //BEGIN:SAVE META INFO
        M_metadataJson["shadingMask"]["nFaces"] = nFaces;     
        if (numOp==1) { 
            M_metadataJson["shadingMask"]["Method"] = "fileSurfaces";
            M_metadataJson["shadingMask"]["nBuildings"] = nBuildings;
        }

        if (numOp==2) { 
            M_metadataJson["shadingMask"]["Method"] = "fileFaces";
            M_metadataJson["shadingMask"]["nMarkers"] = nMarkers;
        }
        //END:SAVE META INFO

        //BEGIN::CREATE BVH
        if (numOp==2) {  // Create a BVH containing all the faces of the buildings
            M_bvh = boundingVolumeHierarchy( _range=M_rangeFaces );  
            auto bvhBuildingTime = toc("BVHs built");
            M_metadataJson["shadingMask"]["Timer"]["BVHs_total_building_time"] = bvhBuildingTime;            
        }    
        //END::CREATE BVH


        LOG(INFO) << "BVHs construction: end";           

    }





template <typename MeshType>
void
ShadingMask<MeshType>::loadMeshDataSubPart3(mesh_ptrtype mesh)
    {
        // [INFO]: refactoring OK for this parts
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

            //BEGIN:SAVE META INFO
            auto dataStructureBuildingTime = toc("Building markers and associated data structures");
            M_metadataJson["shadingMask"]["Method"] = "aggregatedMarkers";
            M_metadataJson["shadingMask"]["nBuildingFaces"] = nFaces;
            M_metadataJson["shadingMask"]["nMarkers"] = nMarkers;
            M_metadataJson["shadingMask"]["Timer"]["DataStructures_building_time"] = dataStructureBuildingTime;
            //BEGIN:SAVE META INFO

            //BEGIN::CREATE BVH
            // Create a BVH containing all the faces of the buildings
            LOG(INFO) << "BVH construction: beginning";
            tic(); 
            M_bvh = boundingVolumeHierarchy( _range=M_rangeFaces );
            auto bvhBuildingTime = toc("BVH built");
            LOG(INFO) << "BVH construction: end";
            M_metadataJson["shadingMask"]["Timer"]["BVH_building_time"] = bvhBuildingTime;
            //END::CREATE BVH
    }



template <typename MeshType>
void
ShadingMask<MeshType>::loadMeshData(mesh_ptrtype mesh, nl::json const& specs)
    {
        // [INFO]: refactoring OK for this parts
        if (specs["/Buildings"_json_pointer].contains("list"))               { loadMeshDataSubPart1(mesh,1); }
        if (specs["/Buildings"_json_pointer].contains("fileVolumes"))        { loadMeshDataSubPart1(mesh,2); }
        if (specs["/Buildings"_json_pointer].contains("fileSurfaces"))       { loadMeshDataSubPart2(mesh,1); } 
        if (specs["/Buildings"_json_pointer].contains("fileFaces"))          { loadMeshDataSubPart2(mesh,2); }
        if (specs["/Buildings"_json_pointer].contains("aggregatedMarkers"))  { loadMeshDataSubPart3(mesh); }
    }
}

