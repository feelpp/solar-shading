
namespace Feel
{
template <typename MeshType>
ShadingMask<MeshType>::ShadingMask(mesh_ptrtype mesh, nl::json const& specs, int intervalsAzimuth, int intervalsAltitude )
{
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
                            list_indices.pop_back();

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
        else if( specs["/Buildings"_json_pointer].contains("fileVolumes")) // a csv containing the volume markers is provided
        {
            std::string buildingName;
            // open the file

            std::ifstream fileVolumes(Environment::expand(specs["Buildings"]["fileVolumes"].get<std::string>()));

            // read, line by line, the building marker
            while ( getline(fileVolumes,buildingName) )
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
    else
    {
        // For each building, save the surface mesh and build the corresponding BVH tree for ray search
        if( specs["/Buildings"_json_pointer].contains("fileSurfaces") ) // a csv containing the surface markers is provided
        {
            std::string buildingName;
            std::ifstream fileSurfaces(Environment::expand(specs["Buildings"]["fileSurfaces"].get<std::string>()));
            std::cout << Environment::expand(specs["Buildings"]["fileSurfaces"].get<std::string>()) << std::endl;
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

            }
        }
        // Store only the view on the surface mesh faces
        else if( specs["/Buildings"_json_pointer].contains("fileFaces") ) // a csv containing the face markers is provided
        {
            std::string faceName;
            std::ifstream fileFaces(Environment::expand(specs["Buildings"]["fileFaces"].get<std::string>()));

            while ( getline(fileFaces,faceName) )
            {
                M_listFaceMarkers.push_back(faceName);
                M_listMarkerFaceEntity[faceName];
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
            }
            // Create a BVH containing all the faces of the buildings
            M_bvh = boundingVolumeHierarchy( _range=M_rangeFaces );
        }
        else if( specs["/Buildings"_json_pointer].contains("aggregatedMarkers"))
        {
            M_rangeFaces = markedelements(mesh,"building"); // it contains all the faces of all buildings
            tic();
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
                        composite_marker.insert(composite_marker.begin(), "building" + markerName.substr(pos, std::string::npos));
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
                    M_listFaceMarkers.push_back(faceName);
                M_listMarkerFaceEntity[faceName].push_back(std::ref(f));
                M_mapEntityToBuildingFace.insert( std::make_pair( f.id(), faceName ) );

            }
            toc("Building markers and associated data structures");
            // Create a BVH containing all the faces of the buildings
            LOG(INFO) << "BVH construction: beginning";
            tic();
            M_bvh = boundingVolumeHierarchy( _range=M_rangeFaces );
            toc("BVH built");
            LOG(INFO) << "BVH construction: end";
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