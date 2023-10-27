#include<shading_mask.hpp>

namespace Feel {

template <typename MeshType>
class ShadingMaskAIdataset: public ShadingMask<MeshType>
{
    using mesh_type = MeshType;
    typedef typename MeshType::ptrtype mesh_ptrtype;

public:

ShadingMaskAIdataset(mesh_ptrtype mesh, nl::json const& specs, int intervalsAzimuth=72, int intervalsAltitude=10):
ShadingMask<MeshType>(mesh, specs, intervalsAzimuth, intervalsAltitude)
{
    // Create the folder structure for the dataset of modified shading masks in absence of one building
    // modShadingMasks/missingBuilding/modifiedShadingMask.csv
    std::string shadingMaskFolder = (boost::filesystem::path(Environment::appRepository())/("modShadingMasks")).string();
    if (!boost::filesystem::exists(shadingMaskFolder))
        boost::filesystem::create_directory(shadingMaskFolder);
    for(auto face : this->M_rangeFaces)
    {
        auto f = boost::unwrap_ref( face );
        std::vector<std::string> composite_marker; // collects all necessary marker substrings to compose the face marker using buildingId and faceId
        for( auto m : f.marker() )
        {
            auto markerName = mesh->markerName(m);
            if((markerName.find("buildingId_") != std::string::npos))
            {
                auto pos = markerName.find_last_of('_');
                // insert the building name at the beginning of the vector
                auto buildingId = markerName.substr(pos+1);
                std::string missingBuildingFolder = shadingMaskFolder+"/"+buildingId;
                if (!boost::filesystem::exists(missingBuildingFolder))
                    boost::filesystem::create_directory(missingBuildingFolder);
            }

        }
    }

}

void computeMasks();

void saveShadingMaskLocalModifications(const Eigen::Ref<const Eigen::MatrixXd>& M, std::string const& building_id, std::string const& mask_name);

std::map< std::string, std::map< std::string,std::vector<double> > > M_building_modified_masks; // buildingN:{ name_modified_mask_buildingM: modifiedMaskMatrix}
};


template <typename MeshType>
void
ShadingMaskAIdataset<MeshType>::computeMasks()
{
    if( this->j_["/Buildings"_json_pointer].contains("fileFaces") || this->j_["/Buildings"_json_pointer].contains("aggregatedMarkers") ) // a csv containing the face markers is provided, or they are computed using aggregated markers
    {
        std::vector<double> random_direction(3);

        std::map<std::string, int> markerLineMap;

        int matrixSize = this->M_azimuthSize * this->M_altitudeSize;

        // Store large vectors containing all the shading mask matrices whose columns are stacked onto each other
        std::vector<double> SM_tables(this->M_listFaceMarkers.size() * matrixSize,0);
        std::vector<double> Angle_tables(this->M_listFaceMarkers.size() * matrixSize,0);

        std::cout << "Allocated SM_tables and Angle_tables of size " << this->M_listFaceMarkers.size() * matrixSize << std::endl;

        int markerNumber = 0;

        // Multithread over rays
        if (this->M_mthreadtype == "ray")
        {
            tic();
            for(auto const &eltWrap : this->M_rangeFaces ) // from each element of the submesh, launch M_Nrays randomly oriented
            {
                auto const& el = unwrap_ref( eltWrap );

                auto elMarker = this->M_mapEntityToBuildingFace.at(el.id());

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
                            auto random_origin = this->get_random_point(el.vertices());

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

                            random_direction = std::get<0>(this->M_raysdirections[initial_index_rays + j]);
                            index_azimuth = std::get<1>(this->M_raysdirections[initial_index_rays + j]);
                            index_altitude = std::get<2>(this->M_raysdirections[initial_index_rays + j]);
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
                                auto rayIntersectionResult =  this->M_bvh->intersect(ray) ;
                                if ( !rayIntersectionResult.empty() )
                                    closer_intersection_element = 1;
                            }
                            // Compute the index associated to the entry to modify
                            // The vector is constituted of M_altitudeSize blocks of M_azimuthSize stacked onto each other
                            int vector_entry = index_azimuth + this->M_azimuthSize*index_altitude;

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
                n_rays_thread.push_back(this->M_Nrays - (this->M_Nthreads-1) * (int)(this->M_Nrays / this->M_Nthreads));

                for(int t= 1; t < this->M_Nthreads; ++t){
                n_rays_thread.push_back( this->M_Nrays / this->M_Nthreads);
                }

                // Used to store the future results
                std::vector< std::future< std::pair<Eigen::VectorXd, Eigen::VectorXd > > > futures;

                for(int t = 0; t < this->M_Nthreads; ++t){

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
            this->M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
            this->M_metadataJson["shadingMask"]["Nthreads"] = this->M_Nthreads;
            this->M_metadataJson["shadingMask"]["NraysPerElement"] = this->M_Nrays;
        }
// Multithread over markers
        if (this->M_mthreadtype == "markers")
        {

                auto multithreading_over_markers = [&](std::vector<std::string> marker_list_thread, int id_thread, int start_index){

                        int index_altitude;
                        int index_azimuth;

                        int initial_index_marker;
                        int i_marker = 0;

                        int len_marker_list_thread = marker_list_thread.size();

                        std::map< std::string, std::map< std::string,std::vector<double> > > local_building_modified_masks;

                        int vector_entry;

                        for( auto const& marker : marker_list_thread)
                        {
                            auto faces_with_marker = this->M_listMarkerFaceEntity[marker];

                            initial_index_marker = start_index + i_marker;

                            auto initial_index_SM = SM_tables.begin() +  initial_index_marker * matrixSize;
                            auto initial_index_Angles = Angle_tables.begin() +  initial_index_marker * matrixSize;

                            // Extract a view from the vectors SM_tables and Angle_tables
                            auto SM_vector = Eigen::Map<Eigen::VectorXd>( &(*initial_index_SM), matrixSize);
                            auto Angle_vector = Eigen::Map<Eigen::VectorXd>( &(*initial_index_Angles), matrixSize);


                            for(auto const& face : faces_with_marker)
                            {
                                for(int j=0;j<this->M_Nrays;j++)
                                {
                                    // Construct the ray emitting from a random point of the element
                                    auto random_origin = this->get_random_point(face.vertices());

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

                                    random_direction = std::get<0>(this->M_raysdirections[j]);
                                    index_azimuth = std::get<1>(this->M_raysdirections[j]);
                                    index_altitude = std::get<2>(this->M_raysdirections[j]);
                                    for(int i=0;i<3;i++)
                                    {
                                        rand_dir(i) = random_direction[i];
                                    }
                                    if(rand_dir.dot(element_normal)>=0)
                                    {
                                        inward_ray=true;
                                    }

                                    BVHRay<mesh_type::nRealDim> ray( origin, rand_dir, 1e-8 );

                                    std::string building_number_first_intersection;
                                    bool more_than_one_intersection;

                                    int closer_intersection_element = -1;
                                    if(inward_ray)
                                    {
                                        closer_intersection_element = 1;
                                    }
                                    else
                                    {
                                        auto rayIntersectionResult =  this->M_bvh->intersect(ray) ;
                                        if ( !rayIntersectionResult.empty() )
                                        {
                                            closer_intersection_element = 1;

                                            auto intersected_building_face = this->M_mapEntityToBuildingFace[rayIntersectionResult.front().primitive().meshEntity().id()];
                                            unsigned first = intersected_building_face.find("building_");
                                            unsigned last = intersected_building_face.find("_face_");
                                            unsigned beginning_building_id = first + std::string("building_").length();
                                            building_number_first_intersection = intersected_building_face.substr(beginning_building_id,last-beginning_building_id);
                                            if(rayIntersectionResult.size()>1)
                                            {
                                                more_than_one_intersection=true;
                                            }
                                            else
                                            {
                                                more_than_one_intersection=false;
                                            }
                                        }

                                    }
                                    // Compute the index associated to the entry to modify
                                    // The vector is constituted of M_altitudeSize blocks of M_azimuthSize stacked onto each other
                                    vector_entry = index_azimuth + this->M_azimuthSize*index_altitude;

                                    // If there is an intersection, increase the shading mask table entry by 1 and augment the angle table by 1 as well
                                    if ( closer_intersection_element >=0 )
                                    {
                                        SM_vector(vector_entry)++;
                                        Angle_vector(vector_entry)++;


                                        if(local_building_modified_masks[building_number_first_intersection][marker].size()!=matrixSize)
                                        {
                                            local_building_modified_masks[building_number_first_intersection][marker].resize(matrixSize);
                                        }

                                        // create a negative shading mask table for the building that shades the current surface
                                        // this way, a new shading mask for the surface is created in the case the aforementioned building were not there

                                        if( !more_than_one_intersection )
                                        // the ray is not shaded in absence of building "building_number_first_intersection"
                                        // otherwise the ray is still shaded by another building which is behind "building_number_first_intersection"
                                        {
                                            local_building_modified_masks[building_number_first_intersection][marker][vector_entry]--;
                                        }

                                    }
                                    else
                                    {
                                        Angle_vector(vector_entry)++;
                                    }
                                }
                            }
                            // Save the partial shading masks only if they are modified
                            for(auto [building_id,mask_structure] : local_building_modified_masks)
                            {
                                for(auto [mask_name,mask_matrix] : mask_structure)
                                {
                                    if (mask_name.find(building_id) != std::string::npos)
                                        continue;

                                    auto initial_index_marker_it = std::find(this->M_listFaceMarkers.begin(),this->M_listFaceMarkers.end(),mask_name);

                                    int initial_index_marker = initial_index_marker_it-this->M_listFaceMarkers.begin();

                                    bool matrix_is_zeros = std::all_of(mask_matrix.begin(), mask_matrix.end(), [](int i) { return i==0; });

                                    if(!matrix_is_zeros)
                                    {
                                        // Subtract from the complete mask the influence of each building
                                        std::transform(mask_matrix.begin(),mask_matrix.end(),SM_vector.begin(),mask_matrix.begin(),std::plus<double>());
                                        // Divide by the total amount of rays launched in that direction
                                        std::transform(mask_matrix.begin(),mask_matrix.end(),Angle_vector.begin(),mask_matrix.begin(),std::divides<double>());

                                        auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(mask_matrix.data(),this->M_azimuthSize, this->M_altitudeSize);
                                        saveShadingMaskLocalModifications(shadingMatrix, building_id, mask_name);
                                    }
                                }
                            }
                            local_building_modified_masks.clear();
                            i_marker += 1;
                            // std::cout << "I_marker " << i_marker << " thread number " << id_thread << " marker " << marker << std::endl;
                        }

                        return true;
                    };
                tic();
                // Store the index of M_listFaceMarkers where each thread will stop its computations
                std::vector<int> marker_threads_list_length;

                marker_threads_list_length.push_back(this->M_listFaceMarkers.size() - (this->M_Nthreads-1) * (int)(this->M_listFaceMarkers.size() / this->M_Nthreads));

                for(int t= 1; t < this->M_Nthreads; ++t){
                marker_threads_list_length.push_back( marker_threads_list_length[t-1] + this->M_listFaceMarkers.size() / this->M_Nthreads);
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
                        marker_thread_lists[t].push_back(this->M_listFaceMarkers[i]);
                    }
                    n0 = n;
                    t += 1;
                    start_index_list[t]=n;
                }

                // Used to store the future results
                // need to
                std::vector< std::future< bool > > futures;

                for(int t= 0; t < this->M_Nthreads; ++t){
                // Launch in parallel asynchronously
                futures.emplace_back(std::async(std::launch::async, multithreading_over_markers, marker_thread_lists[t], t, start_index_list[t]));
                }

                // Collect futures just to allow asynchronous executions
                for( auto& f : futures){
                // Wait for the result to be ready
                auto a =  f.get();
                }
                auto timeComputation = toc("Shading masks computed using raytracing");
                this->M_metadataJson["shadingMask"]["Timer"]["MaskComputation"] = timeComputation;
                this->M_metadataJson["shadingMask"]["Nthreads"] = this->M_Nthreads;
                this->M_metadataJson["shadingMask"]["NraysPerElement"] = this->M_Nrays;

        }

        // Divide the shading mask by the corresponding value of the angle table
        // If an angle combination has not been selected, suppose there is no shadow
        std::transform(SM_tables.begin(),SM_tables.end(),Angle_tables.begin(),SM_tables.begin(),std::divides<double>());

        // Shading mask value 0 means that the surface is not shadowed, value 1 it is fully shadowed
        // Save the shading mask table to a csv file
        if(this->M_saveMasks)
        {
            tic();
            if(this->j_["/Buildings"_json_pointer].contains("fileFaces") )
            for(int i=0; i< this->M_listFaceMarkers.size(); i++)
            {
                std::string building_name = std::to_string(i);
                std::string marker = std::to_string(i);
                auto initial_index_SM = SM_tables.begin() +  i * matrixSize;
                auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(&(*initial_index_SM),this->M_azimuthSize, this->M_altitudeSize);
                this->saveShadingMask(building_name,marker,shadingMatrix.matrix());
            }
            else if(this->j_["/Buildings"_json_pointer].contains("aggregatedMarkers") )
            {
                for(int i=0; i< this->M_listFaceMarkers.size(); i++)
                {
                    std::string building_name = this->M_listFaceMarkers[i];
                    std::string marker = "";
                    auto initial_index_SM = SM_tables.begin() +  i * matrixSize;
                    auto shadingMatrix = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(&(*initial_index_SM),this->M_azimuthSize, this->M_altitudeSize);
                    this->saveShadingMask(building_name,marker,shadingMatrix.matrix());
                }
            }
            auto timeCSVsaving = toc("Mask CSV saved");
            this->M_metadataJson["shadingMask"]["Timer"]["MaskCSVsaving"] = timeCSVsaving;
        }
    }
}

template <typename MeshType>
void
ShadingMaskAIdataset<MeshType>::saveShadingMaskLocalModifications(const Eigen::Ref<const Eigen::MatrixXd>& M, std::string const& building_id, std::string const& mask_name)
{
    int i_long=0;
    int j_alt=0;

    std::string shadingMaskFolder = (boost::filesystem::path(Environment::appRepository())/("modShadingMasks")).string();

    std::ofstream matrix_file;
    std::string matrix_filename = shadingMaskFolder+"/"+building_id+"/SM_Matrix_"+mask_name+".csv";
    matrix_file.open(matrix_filename,std::ios_base::out);
    for(int i=0; i<this->M_azimuthSize; i++)
    {
        for(int j=0; j<this->M_altitudeSize-1; j++)
        {
            matrix_file << M(i,j) << ",";
        }
        matrix_file << M(i,this->M_altitudeSize-1) << "\n";
    }
    matrix_file.close();

}

}