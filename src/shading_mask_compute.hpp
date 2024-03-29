namespace Feel
{

    // Compute shading masks for one building only
    template <typename MeshType>
    void
    ShadingMask<MeshType>::computeMasksOneBuilding(std::string building_name)
    {
        int dim = M_submeshes[building_name]->realDimension();
        std::vector<double> random_direction(dim);

        std::cout << "Submeshes markers" << M_submeshes[building_name]->markerNames() << std::endl;

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

                    auto rays_from_element = [&,marker=marker](int n_rays_thread, int id_thread ){

                        Eigen::MatrixXd SM_table(M_azimuthSize,M_altitudeSize);
                        SM_table.setZero();

                        Eigen::MatrixXd Angle_table(M_azimuthSize,M_altitudeSize);
                        Angle_table.setZero();

                        int index_altitude;
                        int index_azimuth;
                        int initial_index_rays = n_rays_thread * id_thread ;

                        for(int j=0;j<n_rays_thread;j++)
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
                                }
                                auto element_normal = ((p3-p1).head<3>()).cross((p2-p1).head<3>());
                                element_normal.normalize();

                                // Choose the direction randomly among the latitude and azimuth
                                random_direction = std::get<0>(M_raysdirections[initial_index_rays + j]);
                                index_azimuth = std::get<1>(M_raysdirections[initial_index_rays + j]);
                                index_altitude = std::get<2>(M_raysdirections[initial_index_rays + j]);
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
                                    auto rayIntersectionResult =  bvh_building_tree->intersect(_ray=ray) ;
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
                    futures.emplace_back(std::async(std::launch::async, rays_from_element, n_rays_thread[t],t));
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
            if(M_saveMasks)
                saveShadingMask(building_name,marker,shadingMatrix.matrix());

        }
    }

    // Compute shading masks for the buildings in the json file
    template <typename MeshType>
    void
    ShadingMask<MeshType>::computeMasks()
    {
        if( j_["/Buildings"_json_pointer].contains("list") ) // the list of volume markers is provided
        {
            auto markersVolume = j_["Buildings"]["list"].get<std::vector<std::string>>();
            tic();
            for(std::string building_name : markersVolume)
            {
                computeMasksOneBuilding(building_name);//,M_bvh_tree_vector[building_name]);
            }
            auto timeComputation = toc("Shading masks computed using raytracing");
            M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
            M_metadataJson["shadingMask"]["Nthreads"] = M_Nthreads;
            M_metadataJson["shadingMask"]["NraysPerElement"] = M_Nrays;
        }
        else if( j_["/Buildings"_json_pointer].contains("fileVolumes")) // a csv containing the volume markers is provided
        {
            std::string building_name;
            std::ifstream fileVolumes(Environment::expand(j_["Buildings"]["fileVolumes"].get<std::string>()));
            tic();
            // read, line by line, the building marker
            while ( getline(fileVolumes,building_name) )
            {
                computeMasksOneBuilding(building_name);
            }
            auto timeComputation = toc("Shading masks computed using raytracing");
            M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
            M_metadataJson["shadingMask"]["Nthreads"] = M_Nthreads;
            M_metadataJson["shadingMask"]["NraysPerElement"] = M_Nrays;
        }
            // read, line by line, the building marker
        else if( j_["/Buildings"_json_pointer].contains("fileSurfaces") ) // a csv containing the surface markers is provided
        {
            std::string building_name;
            std::ifstream fileSurfaces(Environment::expand(j_["Buildings"]["fileSurfaces"].get<std::string>()));
            tic();
            // read, line by line, the building marker
            while ( getline(fileSurfaces,building_name) )
            {
                computeMasksOneBuilding(building_name);
            }
            auto timeComputation = toc("Shading masks computed using raytracing");
            M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
            M_metadataJson["shadingMask"]["Nthreads"] = M_Nthreads;
            M_metadataJson["shadingMask"]["NraysPerElement"] = M_Nrays;

        }
        else if( j_["/Buildings"_json_pointer].contains("fileFaces") ||  j_["/Buildings"_json_pointer].contains("aggregatedMarkers") ) // a csv containing the face markers is provided, or they are computed using aggregated markers
        {
            std::vector<double> random_direction(3);

            std::map<std::string, int> markerLineMap;

            int matrixSize = M_azimuthSize * M_altitudeSize;

            // Store large vectors containing all the shading mask matrices whose columns are stacked onto each other
            std::vector<double> SM_tables(M_listFaceMarkers.size() * matrixSize,0);
            std::vector<double> Angle_tables(M_listFaceMarkers.size() * matrixSize,0);

            std::cout << "Allocated SM_tables and Angle_tables of size " << M_listFaceMarkers.size() * matrixSize << std::endl;

            int markerNumber = 0;

            // Multithread over rays
            if (M_mthreadtype == "ray")
            {
                tic();
                for(auto const &eltWrap : M_rangeFaces ) // from each element of the submesh, launch M_Nrays randomly oriented
                {
                    auto const& el = unwrap_ref( eltWrap );

                    auto elMarker = M_mapEntityToBuildingFace.at(el.id());

                    auto multithreading_over_rays = [&](int n_rays_thread, int id_thread){

                            Eigen::VectorXd SM_vector(matrixSize);
                            Eigen::VectorXd Angle_vector(matrixSize);

                            SM_vector.setZero();
                            Angle_vector.setZero();

                            int index_altitude;
                            int index_azimuth;

                            int initial_index_rays = n_rays_thread * id_thread ;
                            for(int j=0;j<n_rays_thread;j++)
                            {

                                // Construct the ray emitting from a random point of the element
                                auto random_origin = get_random_point(el.vertices());

                                Eigen::VectorXd rand_dir(3);
                                Eigen::VectorXd p1(3),p2(3),p3(3),origin(3);
                                bool inward_ray=false;

                                for(int i=0;i<3;i++)
                                {
                                    p1(i)=column(el.vertices(), 0)[i];
                                    p2(i)=column(el.vertices(), 1)[i];
                                    p3(i)=column(el.vertices(), 2)[i];
                                    origin(i) = random_origin[i];
                                }
                                auto element_normal = ((p3-p1).head<3>()).cross((p2-p1).head<3>());
                                element_normal.normalize();

                                random_direction = std::get<0>(M_raysdirections[initial_index_rays + j]);
                                index_azimuth = std::get<1>(M_raysdirections[initial_index_rays + j]);
                                index_altitude = std::get<2>(M_raysdirections[initial_index_rays + j]);
                                for(int i=0;i<3;i++)
                                {
                                    rand_dir(i) = random_direction[i];
                                }
                                if(rand_dir.dot(element_normal)>=0)
                                {
                                    inward_ray=true;
                                }

                                BVHRay<mesh_type::nRealDim> ray( origin, rand_dir, 1e-8 );

                                int closer_intersection_element = -1;
                                if(inward_ray)
                                {
                                    closer_intersection_element = 1;
                                }
                                else
                                {
                                    auto rayIntersectionResult =  M_bvh->intersect(_ray=ray) ;
                                    if ( !rayIntersectionResult.empty() )
                                        closer_intersection_element = 1;
                                }
                                // Compute the index associated to the entry to modify
                                // The vector is constituted of M_altitudeSize blocks of M_azimuthSize stacked onto each other
                                int vector_entry = index_azimuth + M_azimuthSize*index_altitude;

                                // If there is an intersection, increase the shading mask table entry by 1 and augment the angle table by 1 as well
                                if ( closer_intersection_element >=0 )
                                {
                                    SM_vector(vector_entry)++;
                                    Angle_vector(vector_entry)++;
                                }
                                else
                                {
                                    Angle_vector(vector_entry)++;
                                }
                            }
                            return std::make_pair(SM_vector,Angle_vector);
                        };

                    // Execute the lambda function on multiple threads using
                    // std::async and std::future to collect the results
                    std::vector<int> n_rays_thread;
                    n_rays_thread.push_back(M_Nrays - (M_Nthreads-1) * (int)(M_Nrays / M_Nthreads));

                    for(int t= 1; t < M_Nthreads; ++t){
                    n_rays_thread.push_back( M_Nrays / M_Nthreads);
                    }

                    // Used to store the future results
                    std::vector< std::future< std::pair<Eigen::VectorXd, Eigen::VectorXd > > > futures;

                    for(int t = 0; t < M_Nthreads; ++t){

                        // Start a new asynchronous task
                        futures.emplace_back(std::async(std::launch::async, multithreading_over_rays, n_rays_thread[t], t));
                    }

                    if( markerLineMap.find(elMarker) == markerLineMap.end())
                    {
                        markerLineMap.insert(std::make_pair(elMarker,markerNumber));
                        markerNumber += 1;
                    }
                    auto initial_index_SM = SM_tables.begin() +  markerLineMap[elMarker] * matrixSize;
                    auto initial_index_Angles = Angle_tables.begin() +  markerLineMap[elMarker] * matrixSize;

                    // Add the tables obtained in threads
                    auto SM_tables_subset = Eigen::Map<Eigen::VectorXd>( &(*initial_index_SM), matrixSize);
                    auto Angle_tables_subset = Eigen::Map<Eigen::VectorXd>( &(*initial_index_Angles), matrixSize);

                    for( auto& f : futures){
                        // Wait for the result to be ready
                        auto two_vectors =  f.get();

                        SM_tables_subset += two_vectors.first;
                        Angle_tables_subset += two_vectors.second;

                    }
                }
                auto timeComputation = toc("Shading masks computed using raytracing");
                M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
                M_metadataJson["shadingMask"]["Nthreads"] = M_Nthreads;
                M_metadataJson["shadingMask"]["NraysPerElement"] = M_Nrays;
            }
// Multithread over markers
            else if (M_mthreadtype == "markers")
            {

                    auto multithreading_over_markers = [&](std::vector<std::string> marker_list_thread, int id_thread, int start_index){

                            int index_altitude;
                            int index_azimuth;

                            int initial_index_marker;
                            int i_marker = 0;

                            int len_marker_list_thread = marker_list_thread.size();

                            int vector_entry;

                            for( auto const& marker : marker_list_thread)
                            {
                                auto faces_with_marker = M_listMarkerFaceEntity[marker];

                                initial_index_marker = start_index + i_marker;

                                auto initial_index_SM = SM_tables.begin() +  initial_index_marker * matrixSize;
                                auto initial_index_Angles = Angle_tables.begin() +  initial_index_marker * matrixSize;

                                // Extract a view from the vectors SM_tables and Angle_tables
                                auto SM_vector = Eigen::Map<Eigen::VectorXd>( &(*initial_index_SM), matrixSize);
                                auto Angle_vector = Eigen::Map<Eigen::VectorXd>( &(*initial_index_Angles), matrixSize);


                                for(auto const& face : faces_with_marker)
                                {
                                    for(int j=0;j<M_Nrays;j++)
                                    {
                                        // Construct the ray emitting from a random point of the element
                                        auto random_origin = get_random_point(face.vertices());

                                        Eigen::VectorXd rand_dir(3);
                                        Eigen::VectorXd p1(3),p2(3),p3(3),origin(3);
                                        bool inward_ray=false;

                                        for(int i=0;i<3;i++)
                                        {
                                            p1(i)=column(face.vertices(), 0)[i];
                                            p2(i)=column(face.vertices(), 1)[i];
                                            p3(i)=column(face.vertices(), 2)[i];
                                            origin(i) = random_origin[i];
                                        }
                                        auto element_normal = ((p3-p1).head<3>()).cross((p2-p1).head<3>());
                                        element_normal.normalize();

                                        // Choose the direction randomly among the latitude and azimuth

                                        random_direction = std::get<0>(M_raysdirections[j]);
                                        index_azimuth = std::get<1>(M_raysdirections[j]);
                                        index_altitude = std::get<2>(M_raysdirections[j]);
                                        for(int i=0;i<3;i++)
                                        {
                                            rand_dir(i) = random_direction[i];
                                        }
                                        if(rand_dir.dot(element_normal)>=0)
                                        {
                                            inward_ray=true;
                                        }

                                        BVHRay<mesh_type::nRealDim> ray( origin, rand_dir, 1e-8 );

                                        int closer_intersection_element = -1;
                                        if(inward_ray)
                                        {
                                            closer_intersection_element = 1;
                                        }
                                        else
                                        {
                                            auto rayIntersectionResult =  M_bvh->intersect(_ray=ray) ;
                                            if ( !rayIntersectionResult.empty() )
                                                closer_intersection_element = 1;
                                        }
                                        // Compute the index associated to the entry to modify
                                        // The vector is constituted of M_altitudeSize blocks of M_azimuthSize stacked onto each other
                                        vector_entry = index_azimuth + M_azimuthSize*index_altitude;

                                        // If there is an intersection, increase the shading mask table entry by 1 and augment the angle table by 1 as well
                                        if ( closer_intersection_element >=0 )
                                        {
                                            SM_vector(vector_entry)++;
                                            Angle_vector(vector_entry)++;
                                        }
                                        else
                                        {
                                            Angle_vector(vector_entry)++;
                                        }
                                    }
                                }

                                i_marker += 1;
                                // std::cout << "I_marker " << i_marker << " thread number " << id_thread << " marker " << marker << std::endl;
                            }
                            return true;
                        };
                    tic();
                    // Store the index of M_listFaceMarkers where each thread will stop its computations
                    std::vector<int> marker_threads_list_length;

                    marker_threads_list_length.push_back(M_listFaceMarkers.size() - (M_Nthreads-1) * (int)(M_listFaceMarkers.size() / M_Nthreads));

                    for(int t= 1; t < M_Nthreads; ++t){
                    marker_threads_list_length.push_back( marker_threads_list_length[t-1] + M_listFaceMarkers.size() / M_Nthreads);
                    }

                    int n0 = 0;
                    int t = 0;

                    std::cout << "Size of marker list per thread " << marker_threads_list_length << std::endl;

                    std::vector<std::vector<std::string>> marker_thread_lists(marker_threads_list_length.size());

                    // Store the index of M_listFaceMarkers from where each thread will start its computations
                    std::vector<int> start_index_list(marker_threads_list_length.size());
                    start_index_list[0]=0;

                    for(auto n : marker_threads_list_length)
                    {
                        for(int i=n0; i< n; i++)
                        {
                            marker_thread_lists[t].push_back(M_listFaceMarkers[i]);
                        }
                        n0 = n;
                        t += 1;
                        start_index_list[t]=n;
                    }

                    // Used to store the future results
                    // need to
                    std::vector< std::future< bool > > futures;

                    for(int t= 0; t < M_Nthreads; ++t){
                    // Launch in parallel asynchronously
                    futures.emplace_back(std::async(std::launch::async, multithreading_over_markers, marker_thread_lists[t], t, start_index_list[t]));
                    }

                    // Collect futures just to allow asynchronous executions
                    for( auto& f : futures){
                    // Wait for the result to be ready
                    auto a =  f.get();
                    }
                    auto timeComputation = toc("Shading masks computed using raytracing");
                    M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
                    M_metadataJson["shadingMask"]["Nthreads"] = M_Nthreads;
                    M_metadataJson["shadingMask"]["NraysPerElement"] = M_Nrays;

            }

            // Divide the shading mask by the corresponding value of the angle table
            // If an angle combination has not been selected, suppose there is no shadow
            std::transform(SM_tables.begin(),SM_tables.end(),Angle_tables.begin(),SM_tables.begin(),std::divides<double>());

            // Shading mask value 0 means that the surface is not shadowed, value 1 it is fully shadowed
            // Save the shading mask table to a csv file
            if(M_saveMasks)
            {
                tic();
                if(j_["/Buildings"_json_pointer].contains("fileFaces") )
                for(int i=0; i< M_listFaceMarkers.size(); i++)
                {
                    std::string building_name = std::to_string(i);
                    std::string marker = std::to_string(i);
                    auto initial_index_SM = SM_tables.begin() +  i * matrixSize;
                    auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(&(*initial_index_SM),M_azimuthSize, M_altitudeSize);
                    saveShadingMask(building_name,marker,shadingMatrix.matrix());
                }
                else if(j_["/Buildings"_json_pointer].contains("aggregatedMarkers") )
                {
                    for(int i=0; i< M_listFaceMarkers.size(); i++)
                    {
                        std::string building_name = M_listFaceMarkers[i];
                        std::string marker = "";
                        auto initial_index_SM = SM_tables.begin() +  i * matrixSize;
                        auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(&(*initial_index_SM),M_azimuthSize, M_altitudeSize);
                        saveShadingMask(building_name,marker,shadingMatrix.matrix());
                    }
                }
                auto timeCSVsaving = toc("Mask CSV saved");
                M_metadataJson["shadingMask"]["Timer"]["MaskCSVsaving"] = timeCSVsaving;
            }
        }
        auto end_computation = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end_computation);
        M_metadataJson["shadingMask"]["Timestamp"]["End"] = strtok(std::ctime(&end_time),"\n");
        saveMetadata();
    }

}
