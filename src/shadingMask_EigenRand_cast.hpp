#include <feel/feelmesh/bvh.hpp>
#include <feel/feelcore/environment.hpp>
#include "../benchmark/extlibs/eigenrand/EigenRand/EigenRand"
// #include <eigenrand/EigenRand/EigenRand>

namespace Feel { 

template <typename MeshType>
class ShadingMaskERC
{
    typedef typename MeshType::ptrtype mesh_ptrtype;
    typedef typename MeshType::trace_mesh_ptrtype tr_mesh_ptrtype;
    typedef typename matrix_node<double>::type matrix_node_type;

public:    
    using value_type = double;
    
    ShadingMaskERC(mesh_ptrtype mesh, nl::json const& specs, int intervalsAzimuth=72, int intervalsAltitude=10 )
    {    
        // Read the number of rays per triangle and the number of threads
        j_ = specs;
        M_Nrays = specs["Nrays"];
        M_Nthreads = specs["Nthreads"].get<int>() ;

        // For each building, save the surface mesh and build the corresponding BVH tree for ray search
        for(std::string buildingName : specs["Buildings"])
        {
            std::cout << fmt::format("{}\n",buildingName);
            auto volumeSubmesh = createSubmesh(_mesh=mesh,_range=markedelements(mesh,buildingName),_update=0);    
            auto surfaceSubmesh = createSubmesh(_mesh=volumeSubmesh,_range=boundaryfaces(volumeSubmesh),_update=0);
            BVHTree<MeshType::nDim> bvhBuilding;
            bvhBuilding.buildPrimitivesInfo(surfaceSubmesh);
            bvhBuilding.buildRootTree();

            M_bvh_tree_vector.insert(std::make_pair( buildingName , bvhBuilding ));
            M_submeshes.insert(std::make_pair( buildingName , surfaceSubmesh ));
        }
        
        // Define the discretization of the azimuth and altitude vectors
        fixAzimuthAltitudeDiscretization(intervalsAzimuth, intervalsAltitude);
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
    void getRandomDirectionSM(std::vector<double> &random_direction, int& index_azimuth, int& index_altitude, Eigen::Rand::UniformRealGen<double>& unif_azi, Eigen::Rand::UniformRealGen<double>& unif_alti, Eigen::Rand::P8_mt19937_64 &urng1, Eigen::Rand::P8_mt19937_64 &urng2)
    {
        int size = random_direction.size();

        if(random_direction.size()==3)
        {
            index_azimuth = static_cast<int>(unif_azi(urng1));
            index_altitude = static_cast<int>(unif_alti(urng2));
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

    Eigen::VectorXd get_random_point(matrix_node_type const& element_points, Eigen::Rand::UniformRealGen<double>& unif_real1, Eigen::Rand::UniformRealGen<double>& unif_real2, Eigen::Rand::P8_mt19937_64 &urng1, Eigen::Rand::P8_mt19937_64 &urng2)
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
                double s = unif_real1(urng1);
                double t = unif_real2(urng2);
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
        if (math::abs(area-elem_area)<1e-5)
            return true;
        else
            return false;
    }
    
    // Compute shading masks for the buildings in the json file
    void computeMasks()
    {
        for(std::string building_name : j_["Buildings"])
        {
            computeMasksOneBuilding(building_name);//,M_bvh_tree_vector[building_name]);
        }
    }

    // Compute shading masks for one building only
    void computeMasksOneBuilding(std::string building_name)//, BVHTree<MeshType::nDim> bvh_tree)
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
                    auto rays_from_element = [&,marker=marker](int n_rays_thread){

                        Eigen::Rand::UniformRealGen<double> unif_gen_azi(0.0,static_cast<double>(M_azimuthSize));
                        Eigen::Rand::UniformRealGen<double> unif_gen_alti(0.0,static_cast<double>(M_altitudeSize));
                        Eigen::Rand::P8_mt19937_64 urng_azi{ std::random_device{}() };
                        Eigen::Rand::P8_mt19937_64 urng_alti{ std::random_device{}() };

                        Eigen::Rand::UniformRealGen<double> unif_gen_real1(0.,1.);
                        Eigen::Rand::UniformRealGen<double> unif_gen_real2(0.,1.);
                        Eigen::Rand::P8_mt19937_64 urng_real1{ std::random_device{}() };
                        Eigen::Rand::P8_mt19937_64 urng_real2{ std::random_device{}() };

                        Eigen::MatrixXd SM_table(M_azimuthSize,M_altitudeSize);
                        SM_table.setZero();

                        Eigen::MatrixXd Angle_table(M_azimuthSize,M_altitudeSize);
                        Angle_table.setZero();
                        
                        int index_altitude;
                        int index_azimuth;                                                 
                        for(int i=0;i<n_rays_thread;i++)
                        {              

                            // Construct the ray emitting from a random point of the element
                            auto random_origin = get_random_point(el.second.vertices(), unif_gen_real1, unif_gen_real2, urng_real1, urng_real2);
                                        
                            Eigen::VectorXd rand_dir(dim); 
                            Eigen::VectorXd p1(dim),p2(dim),p3(dim),origin(3);
                            bool inward_ray=false;
                            if(dim==3)
                            {
                                getRandomDirectionSM(random_direction,index_azimuth,index_altitude, unif_gen_azi, unif_gen_alti, urng_azi, urng_alti);
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
                                if(rand_dir.dot(element_normal)>=0)
                                {
                                    inward_ray=true;
                                }
                                                        
                            }

                            BVHRay ray(origin,rand_dir);
                            
                            int closer_intersection_element = -1;
                            if(inward_ray)
                            {
                                closer_intersection_element = 1;
                            }
                            else
                            {
                                for(auto& [building_name,bvh_building_tree] : M_bvh_tree_vector)
                                {
                                    closer_intersection_element = bvh_building_tree.raySearch(ray,"") ; 
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
        
        std::string matrix_filename = shadingMaskFolder+"/SM_Matrix_"+building_name+"_"+marker_name+"EIGENRAND_CAST.csv";
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

    std::map<std::string,BVHTree<MeshType::nDim>> M_bvh_tree_vector;
    std::map<std::string,tr_mesh_ptrtype> M_submeshes;
    std::map<int,node_type> M_faces_to_normals;

    Eigen::VectorXd M_azimuthAngles, M_altitudeAngles;

    int M_azimuthSize;
    int M_altitudeSize;
    int M_Nrays;
    int M_Nthreads;

    nl::json j_;
};
} // namespace Feel