#include <future>
#include <feel/feelmesh/bvh.hpp>
#include <unordered_map>

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

    using mesh_entity_type = typename std::conditional<MeshType::nDim==MeshType::nRealDim,
                                            entity_range_t<elements_reference_wrapper_t<MeshType>>,
                                            entity_range_t<faces_reference_wrapper_t<MeshType>> >::type;

    using wrapper_type = elements_reference_wrapper_t<MeshType>;
                                        
    typedef typename matrix_node<double>::type matrix_node_type;

public:
    using value_type = double;

    ShadingMask(mesh_ptrtype mesh, nl::json const& specs, int intervalsAzimuth=72, int intervalsAltitude=10 )
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
        for(int i=0; i<M_Nrays; i++)
        {
            getRandomDirectionSM(random_direction,M_gen,M_gen2,index_azimuth,index_altitude);                    
            M_raysdirections[i] = std::make_tuple(random_direction,index_azimuth,index_altitude);
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
        }        
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

    // Get a random point from the surface of the triangle
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
        if( j_["/Buildings"_json_pointer].contains("list") ) // the list of volume markers is provided
        {
            auto markersVolume = j_["Buildings"]["list"].get<std::vector<std::string>>();
            for(std::string building_name : markersVolume)
            {
                computeMasksOneBuilding(building_name);//,M_bvh_tree_vector[building_name]);
            }
        }
        else if( j_["/Buildings"_json_pointer].contains("fileVolumes")) // a csv containing the volume markers is provided
        {
            std::ifstream fileVolumes(Environment::expand(j_["Buildings"]["fileVolumes"].get<std::string>()));
        }
            // read, line by line, the building marker
        else if( j_["/Buildings"_json_pointer].contains("fileSurfaces") ) // a csv containing the surface markers is provided
        {
            std::string building_name;
            std::ifstream fileSurfaces(Environment::expand(j_["Buildings"]["fileSurfaces"].get<std::string>()));

            // read, line by line, the building marker
            while ( getline(fileSurfaces,building_name) )
            {
                computeMasksOneBuilding(building_name);
            }
        }
        else if( j_["/Buildings"_json_pointer].contains("fileFaces") ) // a csv containing the face markers is provided
        {            
            std::vector<double> random_direction(3);            

            std::map<std::string, int> markerLineMap;

            int matrixSize = M_azimuthSize * M_altitudeSize;

            // Store large vectors containing all the shading mask matrices whose columns are stacked onto each other
            std::vector<double> SM_tables(M_listFaceMarkers.size() * matrixSize,0);
            std::vector<double> Angle_tables(M_listFaceMarkers.size() * matrixSize,0);

            int markerNumber = 0;   

            // Multithread over rays
            if (M_mthreadtype == "ray")
            {
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
                            }
                            return true;
                        };

                    // Store the index of M_listFaceMarkers where each thread will stop its computations
                    std::vector<int> marker_threads_list_length;

                    marker_threads_list_length.push_back(M_listFaceMarkers.size() - (M_Nthreads-1) * (int)(M_listFaceMarkers.size() / M_Nthreads));

                    for(int t= 1; t < M_Nthreads; ++t){
                    marker_threads_list_length.push_back( marker_threads_list_length[t-1] + M_listFaceMarkers.size() / M_Nthreads);
                    }

                    int n0 = 0;
                    int t = 0;

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
            }
            
            // Divide the shading mask by the corresponding value of the angle table
            // If an angle combination has not been selected, suppose there is no shadow
            std::transform(SM_tables.begin(),SM_tables.end(),Angle_tables.begin(),SM_tables.begin(),std::divides<double>());            

            // Shading mask value 0 means that the surface is not shadowed, value 1 it is fully shadowed
            // Save the shading mask table to a csv file
            if(M_saveMasks)
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
    }

    // Compute shading masks for one building only
    void computeMasksOneBuilding(std::string building_name)
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

    
    std::map<std::string,std::unique_ptr<BVH<typename tr_mesh_type::element_type>>> M_bvh_tree_vector;
    std::map<std::string,tr_mesh_ptrtype> M_submeshes;
    std::map<int,node_type> M_faces_to_normals;

    std::unique_ptr<BVH<typename tr_mesh_type::element_type>> M_bvh;

    std::unordered_map<int,std::string> M_mapEntityToBuildingFace;
    wrapper_type M_rangeFaces;
    std::vector<std::string> M_listFaceMarkers;

    Eigen::VectorXd M_azimuthAngles, M_altitudeAngles;

    std::vector< std::tuple< std::vector<double> /* random direction */,int /* index azimuth */, int /* index altitude */ > > M_raysdirections;

    std::map<std::string /* marker name */,std::vector<typename tr_mesh_type::element_type> /* list of references to face entities */> M_listMarkerFaceEntity;

    int M_azimuthSize;
    int M_altitudeSize;
    int M_Nrays;
    int M_Nthreads;
    std::string M_mthreadtype;
    bool M_saveMasks;

    nl::json j_;

    std::random_device M_rd;
    std::random_device M_rd2;
    std::mt19937 M_gen;
    std::mt19937 M_gen2;
};
} // namespace Feel
