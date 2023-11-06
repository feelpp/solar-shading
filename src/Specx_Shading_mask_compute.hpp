



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
//KILL THIS PART
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
                //std::cout<<"[SPECX INFO] : ******Save Shading Mask in one building.\n";
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

            std::cout<<"[SPECX INFO] : ***Section List\n";
            auto markersVolume = j_["Buildings"]["list"].get<std::vector<std::string>>();
            
            int NbLoop=0; int NbTh=0; int NbCoor=0;
            GetSpecxPreprocessingParameters(NbObjects,SpecxNbThreadDesired,NbLoop,NbTh,NbCoor);
            std::cout<<"[SPECX INFO] : SpecxNbThreadDesired="<<SpecxNbThreadDesired<<"\n";
            std::cout<<"[SPECX INFO] : Nb Objects="<<NbObjects<<"\n";
            std::cout<<"[SPECX INFO] : Nb Loops="<<NbLoop<<" Th="<<NbTh<<" Coor="<<NbCoor<<"\n";
            std::cout<<"[SPECX INFO] : NbThreads="<<NbTh<<" used\n";
            SpRuntime runtime1(NbTh);

            tic();


                std::string buildingName;
                int initVal = 0;

                int NbidxN=NbTh;
                for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                    if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }
                    for(int idx = 0 ; idx < NbTh ; ++idx){
                            int Index=idxLoop*NbidxN+idx;
                            std::cout <<"Index:"<<Index<<"\n";
                            buildingName=markersVolume[Index];
                            std::cout << "[INFO SPECS] : NAME =>"<<buildingName<<"\n";

                            std::promise<int> promise1;

                            runtime1.task(SpRead(buildingName), [&promise1](const std::string & bn){
                                promise1.get_future().get();
                            }).setTaskName("Start Task");

                            runtime1.task(
                                SpPriority(1),
                                SpRead(buildingName),
                                [&](const std::string & bn) -> bool {
                                computeMasksOneBuilding(bn);
                                return true;
                                }
                            ).setTaskName("Task 2V Compute :"+buildingName);

                            promise1.set_value(0);
                            std::cout << "[INFO SPECS] : OK END"<<"\n"; 
                            
                    }//End for 
                    runtime1.waitAllTasks();
                }//End For idxLoop

            /*
            for(std::string building_name : markersVolume)
            {
                computeMasksOneBuilding(building_name);//,M_bvh_tree_vector[building_name]);
            }
            */

            runtime1.waitAllTasks();
            runtime1.stopAllThreads();
            std::cout<<"[SPECX INFO] : STOP ALL SPECX THREAD\n";
            runtime1.generateDot("Runtime1.dot",true);
            runtime1.generateTrace("Runtime1.svg");


            auto timeComputation = toc("Shading masks computed using raytracing");
            M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
            M_metadataJson["shadingMask"]["Nthreads"] = M_Nthreads;
            M_metadataJson["shadingMask"]["NraysPerElement"] = M_Nrays;
        }
        else if( j_["/Buildings"_json_pointer].contains("fileVolumes")) // a csv containing the volume markers is provided
        {
            std::cout<<"[SPECX INFO] : ***Section Volumes\n";
            std::string building_name;
            std::ifstream fileVolumes(Environment::expand(j_["Buildings"]["fileVolumes"].get<std::string>()));
            std::cout<<"[SPECX INFO] : List\n";
            std::vector<std::string> ListObjects;
            ListObjects=GetListNameObjects(Environment::expand(j_["Buildings"]["fileVolumes"].get<std::string>()));
            NbObjects=ListObjects.size();

            //SpecxNbThreadDesired=std::sqrt(NbObjects);
       
            int NbLoop=0; int NbTh=0; int NbCoor=0;
            GetSpecxPreprocessingParameters(NbObjects,SpecxNbThreadDesired,NbLoop,NbTh,NbCoor);
            std::cout<<"[SPECX INFO] : SpecxNbThreadDesired="<<SpecxNbThreadDesired<<"\n";
            std::cout<<"[SPECX INFO] : Nb Objects="<<NbObjects<<"\n";
            std::cout<<"[SPECX INFO] : Nb Loops="<<NbLoop<<" Th="<<NbTh<<" Coor="<<NbCoor<<"\n";
            std::cout<<"[SPECX INFO] : NbThreads="<<NbTh<<" used\n";
            SpRuntime runtime1(NbTh);

            //std::string buildingName;

                int NbidxN=NbTh;
                for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                    if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }
                    for(int idx = 0 ; idx < NbTh ; ++idx){
                            int Index=idxLoop*NbidxN+idx;
                            runtime1.task(
                                SpRead(Index),
                                [&](const int & k) -> bool { 
                                    computeMasksOneBuilding(ListObjects[k]);  
                                    //usleep(10); 
                                    return true; 
                                }
                            ).setTaskName("Op("+std::to_string(Index)+")");  
                            usleep(10);
                            std::atomic_int counter(0);

                    }//End for 
                    runtime1.waitAllTasks();
                    usleep(100);
                }//End For idxLoop


                  


                        
                            
            


            tic();
            // read, line by line, the building marker

            /*
            while ( getline(fileVolumes,building_name) )
            {
                computeMasksOneBuilding(building_name);
            }
            */

            runtime1.stopAllThreads();
            std::cout<<"[SPECX INFO] : STOP ALL SPECX THREAD\n";
            runtime1.generateDot("Runtime1.dot",true);
            runtime1.generateTrace("Runtime1.svg");

            auto timeComputation = toc("Shading masks computed using raytracing");
            M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
            M_metadataJson["shadingMask"]["Nthreads"] = M_Nthreads;
            M_metadataJson["shadingMask"]["NraysPerElement"] = M_Nrays;
        }
            // read, line by line, the building marker
        else if( j_["/Buildings"_json_pointer].contains("fileSurfaces") ) // a csv containing the surface markers is provided
        {
            std::cout<<"[SPECX INFO] : ***Section Surfaces\n";
            std::string building_name;
            std::ifstream fileSurfaces(Environment::expand(j_["Buildings"]["fileSurfaces"].get<std::string>()));


            //std::cout<<Environment::expand(j_["Buildings"]["fileSurfaces"].get<std::string>())<<"\n";

            std::cout<<"[SPECX INFO] : List\n";
            std::vector<std::string> ListObjects;
            ListObjects=GetListNameObjects(Environment::expand(j_["Buildings"]["fileSurfaces"].get<std::string>()));
            NbObjects=ListObjects.size();
            //for(int k = 0 ; k < NbObjects ; ++k){ std::cout<<ListObjects[k]<<"\n";} //CTRL
            
            int NbLoop=0; int NbTh=0; int NbCoor=0;
            GetSpecxPreprocessingParameters(NbObjects,SpecxNbThreadDesired,NbLoop,NbTh,NbCoor);
            std::cout<<"[SPECX INFO] : SpecxNbThreadDesired="<<SpecxNbThreadDesired<<"\n";
            std::cout<<"[SPECX INFO] : Nb Objects="<<NbObjects<<"\n";
            std::cout<<"[SPECX INFO] : Nb Loops="<<NbLoop<<" Th="<<NbTh<<" Coor="<<NbCoor<<"\n";
            std::cout<<"[SPECX INFO] : NbThreads="<<NbTh<<" used\n";
            SpRuntime runtime1(NbTh);

            tic();

            std::string buildingName;
                int initVal = 0;

                int NbidxN=NbTh;
                for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                    if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }
                    for(int idx = 0 ; idx < NbTh ; ++idx){
                            int Index=idxLoop*NbidxN+idx;
                            //std::cout <<"Index:"<<Index<<"\n";
                            buildingName=ListObjects[Index];
                            runtime1.task(
                                SpRead(buildingName),
                                [&](const std::string & bn) -> bool { computeMasksOneBuilding(bn); return true; }
                            ).setTaskName("Task 2V Compute :"+buildingName);
                            
                    }//End for 
                    runtime1.waitAllTasks();
                }//End For idxLoop
            

            // read, line by line, the building marker
            /*
            while ( getline(fileSurfaces,building_name) )
            {
                computeMasksOneBuilding(building_name);
            }
            */
            
            runtime1.stopAllThreads();
            std::cout<<"[SPECX INFO] : STOP ALL SPECX THREAD\n";
            runtime1.generateDot("Runtime1.dot",true);
            runtime1.generateTrace("Runtime1.svg");

            auto timeComputation = toc("Shading masks computed using raytracing");
            M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
            M_metadataJson["shadingMask"]["Nthreads"] = M_Nthreads;
            M_metadataJson["shadingMask"]["NraysPerElement"] = M_Nrays;
            
        }





        else if( j_["/Buildings"_json_pointer].contains("fileFaces") ||  j_["/Buildings"_json_pointer].contains("aggregatedMarkers") ) // a csv containing the face markers is provided, or they are computed using aggregated markers
        {            


            std::cout<<"[SPECX INFO] : ***Section Aggregate Markers\n";

            std::vector<double> random_direction(3);            

            std::map<std::string, int> markerLineMap;

            int matrixSize = M_azimuthSize * M_altitudeSize;

            // Store large vectors containing all the shading mask matrices whose columns are stacked onto each other
            std::vector<double> SM_tables(M_listFaceMarkers.size() * matrixSize,0);
            std::vector<double> Angle_tables(M_listFaceMarkers.size() * matrixSize,0);

            std::cout << "Allocated SM_tables and Angle_tables of size " << M_listFaceMarkers.size() * matrixSize << std::endl;

            int markerNumber = 0;  

            //Only used to control if it is ok or not. 
            if (!QMAKE_WITH_SPECX) {
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
                                                    auto rayIntersectionResult =  M_bvh->intersect(ray) ;
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
                                                            auto rayIntersectionResult =  M_bvh->intersect(ray) ;
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
            }

            if (QMAKE_WITH_SPECX) {
                // Multithread over rays
                if (M_mthreadtype == "ray")
                {

                    std::cout<<"[SPECX INFO] : ***Section Aggregate Markers RAY\n";
                    
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
                                        auto rayIntersectionResult =  M_bvh->intersect(ray) ;
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


                    //BLOCK001:BEGIN
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

                            //std::cout<<start_index_list[t]<<"\n"; 
                        }

                    //BLOCK002:END
                    
                    //QViewInfoSpecx=true;
                    //SpecxNbThreadDesired=20;
                    //SpecxNbThreadDesired=10;
                    //SpecxNbThreadDesired=5;
                    //SpecxNbThreadDesired=M_Nthreads;
                    NbObjects=M_Nthreads;

                    int NbLoop=0; int NbTh=0; int NbCoor=0;
                    GetSpecxPreprocessingParameters(NbObjects,SpecxNbThreadDesired,NbLoop,NbTh,NbCoor);
                    if (QViewInfoSpecx) {
                        std::cout<<"[SPECX INFO] : SpecxNbThreadDesired="<<SpecxNbThreadDesired<<"\n";
                        std::cout<<"[SPECX INFO] : Nb Objects="<<NbObjects<<"\n";
                        std::cout<<"[SPECX INFO] : Nb Loops="<<NbLoop<<" Th="<<NbTh<<" Coor="<<NbCoor<<"\n";
                        std::cout<<"[SPECX INFO] : NbThreads="<<NbTh<<" used\n";
                    }
                    SpRuntime runtime1(NbTh);
                    NbTh=runtime1.getNbThreads(); //CTRL if OK
                    if (QViewInfoSpecx) { std::cout<<"[SPECX INFO] : Num Thread in markers="<<NbTh<<"\n"; }

                    //std::cout<<"\nM_Nthreads="<<M_Nthreads<<"\n";
                    

                    /* BEGIN LAMBDA FUNCTION *********************************************************************/
                                    auto multithreading_over_markers = [&](std::vector<std::string> marker_list_thread, int id_thread, int start_index){

                                                int index_altitude;
                                                int index_azimuth;
                                                
                                                int initial_index_marker;
                                                int i_marker = 0;

                                                int len_marker_list_thread = marker_list_thread.size();

                                                int vector_entry;

                                                //std::cout <<"marker_list_thread:"<<marker_list_thread[id_thread]<<"\n";
                                                //std::cout <<"Index id_thread:"<<id_thread<<"\n";
                                                //std::cout <<"Start Index:"<<start_index<<"\n";
                                                //std::cout <<"marker_list_thread.size="<<marker_list_thread.size()<<"\n";

                                                for( auto const& marker : marker_list_thread)
                                                {
                                                    auto faces_with_marker = M_listMarkerFaceEntity[marker];
                                                    
                                                    initial_index_marker = start_index + i_marker;

                                                    auto initial_index_SM = SM_tables.begin() +  initial_index_marker * matrixSize;
                                                    auto initial_index_Angles = Angle_tables.begin() +  initial_index_marker * matrixSize;

                                                    // Extract a view from the vectors SM_tables and Angle_tables
                                                    auto SM_vector = Eigen::Map<Eigen::VectorXd>( &(*initial_index_SM), matrixSize);
                                                    auto Angle_vector = Eigen::Map<Eigen::VectorXd>( &(*initial_index_Angles), matrixSize);

                                                    //std::cout <<"faces_with_marker="<<faces_with_marker.size()<<"\n";


                                                    //BEGIN:SUB PART
                                                
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
                                                                auto rayIntersectionResult =  M_bvh->intersect(ray) ;
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
                                                    
                                                    //END:SUB PART




                                                    i_marker += 1;
                                                    // std::cout << "I_marker " << i_marker << " thread number " << id_thread << " marker " << marker << std::endl;
                                                }
                                                return true;
                                            };


                    /*END LAMBDA FUNCTION *********************************************************************/

                    usleep(10);
                    //getchar();
                        
                    tic();
                    SpTimer timerTask1;


                        int NbidxN=NbTh;
                        for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                            if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }
                            for(int idx = 0 ; idx < NbTh ; ++idx){
                                    int Index=idxLoop*NbidxN+idx;
                                    runtime1.task(
                                    SpRead(Index),
                                    [&](const int & k) -> bool {
                                        //if (QViewInfoSpecx) { std::cout <<"ST:"<<k<<"\n"; }
                                        //std::cout <<" "<<marker_thread_lists[k]<<"\n";
                                        multithreading_over_markers(marker_thread_lists[k],k,start_index_list[k]);
                                        //if (QViewInfoSpecx) { std::cout <<"EE\n"; }
                                    return true;
                                    }
                                    ).setTaskName("Op("+std::to_string(Index)+")");
                                    usleep(10);
                                    std::atomic_int counter(0);
                                    
                            }//End for
                            runtime1.waitAllTasks();
                            usleep(100);
                        }//End For idxLoop



                
                    runtime1.waitAllTasks();
                    runtime1.stopAllThreads();
                    timerTask1.stop();


                    auto timeSpecxTask1=timerTask1.getElapsed();
                    M_metadataJson["shadingMask"]["Timer Specx"]["MaskSaving Task 1"] = timeSpecxTask1;

                    if (QViewInfoSpecx) {
                        std::cout<<"\n";
                        std::cout<<"\n";
                        std::cout<<"[SPECX INFO] : STOP ALL SPECX THREAD\n";
                        std::cout<<"[SPECX INFO] : Task Time = "<<timeSpecxTask1<< " s" << std::endl; 
                    }

                    //std::cout<<"=END MARKERS============================================\n";
                    //getchar();

    /*
                //KILL this part
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
    */
                        auto timeComputation = toc("Shading masks computed using raytracing");
                        M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
                        M_metadataJson["shadingMask"]["Nthreads"] = M_Nthreads;
                        M_metadataJson["shadingMask"]["NraysPerElement"] = M_Nrays;

                        if (QSaveSpecxGenerateReports) {
                            tic();
                            if (QViewInfoSpecx) { std::cout<<"[SPECX INFO] : ******Save Generate.\n"; }
                            runtime1.generateDot("Runtime1.dot",true);
                            runtime1.generateTrace("Runtime1.svg");
                            auto timeSpecxGeneration = toc("Generate 1 SVG and DOT saved");
                            M_metadataJson["shadingMask"]["Timer Specx"]["MaskComputation Generate"] = timeSpecxGeneration;
                        }
                }    
            }

            // Divide the shading mask by the corresponding value of the angle table
            // If an angle combination has not been selected, suppose there is no shadow
            std::transform(SM_tables.begin(),SM_tables.end(),Angle_tables.begin(),SM_tables.begin(),std::divides<double>());            

            // Shading mask value 0 means that the surface is not shadowed, value 1 it is fully shadowed
            // Save the shading mask table to a csv file


            


            std::cout<<"[SPECX INFO] : ******Save Shading Mask.\n";
            std::cout<<"[SPECX INFO] : SAVE RESULTS IN PARALLEL MODE\n";

            if(M_saveMasks)
            {
                //QViewInfoSpecx=true;

                tic();
                if(j_["/Buildings"_json_pointer].contains("fileFaces") )
                {
                    if (QViewInfoSpecx) { 
                            std::cout<<"[SPECX INFO] : ******Save file Faces.\n"; 
                            std::cout<<"[SPECX INFO] : "<<M_listFaceMarkers.size()<<"\n"; 
                    }

                    if (QSaveWithSpecx)
                    {
                        NbObjects=M_listFaceMarkers.size();
                        int NbLoop=0; int NbTh=0; int NbCoor=0;
                        GetSpecxPreprocessingParameters(NbObjects,SpecxSaveNbThreadDesired,NbLoop,NbTh,NbCoor);
                        if (QViewInfoSpecx) {
                            std::cout<<"[SPECX INFO] : SpecxSaveNbThreadDesired="<<SpecxSaveNbThreadDesired<<"\n";
                            std::cout<<"[SPECX INFO] : Nb Objects="<<NbObjects<<"\n";
                            std::cout<<"[SPECX INFO] : Nb Loops="<<NbLoop<<" Th="<<NbTh<<" Coor="<<NbCoor<<"\n";
                            std::cout<<"[SPECX INFO] : NbThreads="<<NbTh<<" used\n";
                        }
                        SpRuntime runtime5(NbTh);
                        NbTh=runtime5.getNbThreads(); //CTRL if OK
                        if (QViewInfoSpecx) { std::cout<<"[SPECX INFO] : Num Thread in markers="<<NbTh<<"\n"; }
                        usleep(10);
                        SpTimer timerTask5;

                        int NbidxN=NbTh;
                        for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                            if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }
                            for(int idx = 0 ; idx < NbTh ; ++idx){
                                    int Index=idxLoop*NbidxN+idx;
                                    runtime5.task(
                                    SpRead(Index),
                                    [&](const int & k) -> bool {
                                        //if (QViewInfoSpecx) { std::cout <<"ST:"<<k<<"\n"; }
                                        std::string building_name = std::to_string(k);
                                        std::string marker = std::to_string(k);
                                        auto initial_index_SM = SM_tables.begin() +  k * matrixSize;
                                        auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(&(*initial_index_SM),M_azimuthSize, M_altitudeSize);                
                                        saveShadingMask(building_name,marker,shadingMatrix.matrix());
                                        //if (QViewInfoSpecx) { std::cout <<"EE\n"; }
                                    return true;
                                    }
                                    ).setTaskName("Op("+std::to_string(Index)+")");
                                    
                                    if (QSpecxLockConfigWaitSave) { usleep(10); std::atomic_int counter(0); }    
                                   
                                    
                            }//End for
                            if (QSpecxLockConfigWaitSave) { runtime5.waitAllTasks(); usleep(100); }
                        }//End For idxLoop


                        runtime5.waitAllTasks();
                        runtime5.stopAllThreads();
                        timerTask5.stop();

                        auto timeSpecxTask5 = timerTask5.getElapsed();
                        M_metadataJson["shadingMask"]["Timer Specx"]["Time save parallel mode"] = timeSpecxTask5;

                        if (QViewInfoSpecx) {
                            std::cout<<"\n";
                            std::cout<<"\n";
                            std::cout<<"[SPECX INFO] : STOP ALL SPECX THREAD\n";
                            std::cout<<"[SPECX INFO] : Task Time = "<<timeSpecxTask5<< " s" << std::endl;
                        }

                        if (QSaveSpecxGenerateReports) {
                            tic();
                            if (QViewInfoSpecx) { std::cout<<"[SPECX INFO] : ******Save Generate.\n"; }
                            runtime5.generateDot("Runtime5.dot",true);
                            runtime5.generateTrace("Runtime5.svg");
                            auto timeSpecxGeneration = toc("Generate 5 SVG and DOT saved");
                            M_metadataJson["shadingMask"]["Timer Specx"]["MaskSaving Generate"] = timeSpecxGeneration;
                        }
                    }
                    else
                    {
                        for(int i=0; i< M_listFaceMarkers.size(); i++)
                        {
                            std::string building_name = std::to_string(i);
                            std::string marker = std::to_string(i);
                            auto initial_index_SM = SM_tables.begin() +  i * matrixSize;
                            auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(&(*initial_index_SM),M_azimuthSize, M_altitudeSize);                
                            saveShadingMask(building_name,marker,shadingMatrix.matrix());
                        }
                    }
                }


                else if(j_["/Buildings"_json_pointer].contains("aggregatedMarkers") )
                {

                    if (QViewInfoSpecx) { 
                        std::cout<<"[SPECX INFO] : ******Save Aggregated Markers.\n"; 
                        std::cout<<"[SPECX INFO] : "<<M_listFaceMarkers.size()<<"\n"; 
                    }

                    if (QSaveWithSpecx)
                    {
                        
                        NbObjects=M_listFaceMarkers.size();
                        int NbLoop=0; int NbTh=0; int NbCoor=0;
                        GetSpecxPreprocessingParameters(NbObjects,SpecxSaveNbThreadDesired,NbLoop,NbTh,NbCoor);
                        if (QViewInfoSpecx) {
                            std::cout<<"[SPECX INFO] : SpecxSaveNbThreadDesired="<<SpecxSaveNbThreadDesired<<"\n";
                            std::cout<<"[SPECX INFO] : Nb Objects="<<NbObjects<<"\n";
                            std::cout<<"[SPECX INFO] : Nb Loops="<<NbLoop<<" Th="<<NbTh<<" Coor="<<NbCoor<<"\n";
                            std::cout<<"[SPECX INFO] : NbThreads="<<NbTh<<" used\n";
                        }
                        SpRuntime runtime5(NbTh);
                        NbTh=runtime5.getNbThreads(); //CTRL if OK
                        if (QViewInfoSpecx) { std::cout<<"[SPECX INFO] : Num Thread in markers="<<NbTh<<"\n"; }
                        usleep(10);
                        SpTimer timerTask5;

                        int NbidxN=NbTh;
                        for(int idxLoop = 0 ; idxLoop < NbLoop ; ++idxLoop){
                            if (idxLoop==NbLoop-1) { NbTh=NbTh+NbCoor; }
                            for(int idx = 0 ; idx < NbTh ; ++idx){
                                    int Index=idxLoop*NbidxN+idx;
                                    runtime5.task(
                                    SpRead(Index),
                                    [&](const int & k) -> bool {
                                        //if (QViewInfoSpecx) { std::cout <<"ST:"<<k<<"\n"; }
                                        std::string building_name = M_listFaceMarkers[k];
                                        std::string marker = "";
                                        auto initial_index_SM = SM_tables.begin() +  k * matrixSize;
                                        auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(&(*initial_index_SM),M_azimuthSize, M_altitudeSize);                
                                        saveShadingMask(building_name,marker,shadingMatrix.matrix());
                                        //if (QViewInfoSpecx) { std::cout <<"EE\n"; }
                                    return true;
                                    }
                                    ).setTaskName("Op("+std::to_string(Index)+")");

                                    if (QSpecxLockConfigWaitSave) { usleep(10); std::atomic_int counter(0); }    
                                    //usleep(1);                                
                            }//End for

                            if (QSpecxLockConfigWaitSave) { runtime5.waitAllTasks(); usleep(10); }
                            
                        }//End For idxLoop


                        runtime5.waitAllTasks();
                        runtime5.stopAllThreads();
                        timerTask5.stop();


                        auto timeSpecxTask5 = timerTask5.getElapsed();
                        M_metadataJson["shadingMask"]["Timer Specx"]["Time save parallel mode"] = timeSpecxTask5;

                        if (QViewInfoSpecx) {
                            std::cout<<"\n";
                            std::cout<<"\n";
                            std::cout<<"[SPECX INFO] : STOP ALL SPECX THREAD\n";
                            std::cout<<"[SPECX INFO] : Task Time = "<<timeSpecxTask5<< " s" << std::endl;
                        }
                        
                        if (QSaveSpecxGenerateReports) {
                            tic();
                            if (QViewInfoSpecx) { std::cout<<"[SPECX INFO] : ******Save Generate.\n"; }
                            runtime5.generateDot("Runtime5.dot",true);
                            runtime5.generateTrace("Runtime5.svg");
                            auto timeSpecxGeneration = toc("Generate 5 SVG and DOT saved");
                            M_metadataJson["shadingMask"]["Timer Specx"]["MaskSaving Generate"] = timeSpecxGeneration;
                        }
                        
                    }
                    else
                    {
                        //Part to be deleted in the future...
                        for(int i=0; i< M_listFaceMarkers.size(); i++)
                        {
                            std::string building_name = M_listFaceMarkers[i];
                            std::string marker = "";
                            auto initial_index_SM = SM_tables.begin() +  i * matrixSize;
                            auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(&(*initial_index_SM),M_azimuthSize, M_altitudeSize);                
                            saveShadingMask(building_name,marker,shadingMatrix.matrix());
                        }

                    }


                    //CTRL mode normal
                    if((M_saveMasks) && (QCTRL_SAVE_NORMAL))
                    {
                        tic();
                        cout<<"\n";
                        cout<<"[SPECX INFO] : SAVE RESULTS IN SEQUENTIAL MODE ===> ";
                        CONSOLE_SaveCursorPosition();
                        CONSOLE_CursorHidden(true); 

                        if(j_["/Buildings"_json_pointer].contains("fileFaces") )
                        {
                            for(int i=0; i< M_listFaceMarkers.size(); i++)
                            {
                                CONSOLE_PRINT_PERCENTAGE(i+1,M_listFaceMarkers.size());
                                std::string building_name = std::to_string(i);
                                std::string marker = std::to_string(i);
                                auto initial_index_SM = SM_tables.begin() +  i * matrixSize;
                                auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(&(*initial_index_SM),M_azimuthSize, M_altitudeSize);                
                                saveShadingMask("CTRL_"+building_name,marker,shadingMatrix.matrix());
                            }
                        }
                        else if(j_["/Buildings"_json_pointer].contains("aggregatedMarkers") )
                        {
                            for(int i=0; i< M_listFaceMarkers.size(); i++)
                            {
                                CONSOLE_PRINT_PERCENTAGE(i+1,M_listFaceMarkers.size());
                                std::string building_name = M_listFaceMarkers[i];
                                std::string marker = "";
                                auto initial_index_SM = SM_tables.begin() +  i * matrixSize;
                                auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(&(*initial_index_SM),M_azimuthSize, M_altitudeSize);                
                                saveShadingMask("CTRL_"+building_name,marker,shadingMatrix.matrix());
                            }
                        }
                        CONSOLE_CursorHidden(false); 
                        cout<<"\n"; 
                        M_metadataJson["shadingMask"]["Timer Specx"]["Time save sequential mode"] = toc("Time save sequential mode");
                    }

                    if((M_saveMasks) && (QCTRL_DATA))
                    {
                        tic();
                        cout<<"[SPECX INFO] : TASK CTRL DATA Level ===> ";
                        //CTRL 
                        if ((j_["/Buildings"_json_pointer].contains("aggregatedMarkers") ) && (true)) 
                        {
                            bool QTestOKLevel1=true;
                            bool QTestOKLevel2=true;

                            CONSOLE_SaveCursorPosition();
                            CONSOLE_CursorHidden(true); 

                            for(int i=0; i< M_listFaceMarkers.size(); i++)
                                {
                                    std::string building_name = M_listFaceMarkers[i];
                                    std::string marker = "";
                                    std::string shadingMaskFolder = (boost::filesystem::path(Environment::appRepository())/("shadingMasks")).string();
                                    QTestOKLevel1=QTestOKLevel1 && (testShadingMaskComparisonLevel1(shadingMaskFolder,building_name,marker));

                                    //QTestOKLevel2=QTestOKLevel2 && (testShadingMaskComparisonLevel2(shadingMaskFolder,building_name,marker));

                                    CONSOLE_PRINT_PERCENTAGE(i+1,M_listFaceMarkers.size());

                                }
                            //cout<<"\n";
                            if (QTestOKLevel1) { cout<<" >TRUE\n"; } else { cout<<" >ERROR\n"; }
                            M_metadataJson["shadingMask"]["Timer Specx"]["Time CTRL Masks"] = toc("Time CTRL Masks");
                        }
                        CONSOLE_CursorHidden(false); 

                    }

                    if (QViewInfoSpecx) { std::cout<<"[SPECX INFO] : ******Saved Aggregated Markers.\n"; }
                }


                //auto timeCSVsaving = toc("Mask CSV saved");
                M_metadataJson["shadingMask"]["Timer"]["MaskCSVsaving"] = toc("Masks CSV saved");


            }
        }
        auto end_computation = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end_computation);
        M_metadataJson["shadingMask"]["Timestamp"]["End"] = strtok(std::ctime(&end_time),"\n");

        //std::cout<<"[SPECX INFO] : ******Section Save Metadata\n";
        saveMetadata();   
    }

}


 