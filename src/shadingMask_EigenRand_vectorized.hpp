#include <feel/feelmesh/bvh.hpp>
#include <feel/feelcore/environment.hpp>
#include "../benchmark/extlibs/eigenrand/EigenRand/EigenRand"
// #include <eigenrand/EigenRand/EigenRand>

namespace Feel { 

template <typename MeshType>
class ShadingMaskERV
{
    typedef typename MeshType::ptrtype mesh_ptrtype;
    typedef typename MeshType::trace_mesh_ptrtype tr_mesh_ptrtype;
    typedef typename matrix_node<double>::type matrix_node_type;

public:    
    using value_type = double;
    
    ShadingMaskERV(mesh_ptrtype mesh, nl::json const& specs, int intervalsAzimuth=72, int intervalsAltitude=10 )
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

    
    Eigen::MatrixXd get_random_directions(Eigen::VectorXi & index_azimuth, Eigen::VectorXi & index_altitude, int n_rays_thread, Eigen::Rand::UniformIntGen<int>& unif_gen_azi, Eigen::Rand::UniformIntGen<int>& unif_gen_alti, Eigen::Rand::P8_mt19937_64& urng_azi, Eigen::Rand::P8_mt19937_64& urng_alti )
    {
        Eigen::MatrixXd random_directions(n_rays_thread,3);

        index_azimuth = unif_gen_azi.generate<Eigen::VectorXi>(n_rays_thread,1,urng_azi);
        index_altitude = unif_gen_alti.generate<Eigen::VectorXi>(n_rays_thread,1,urng_alti);
        
        // Some critical errors occur when generating numbers with the pcg generator
        // Some generated numbers are way out of the range of the azimuth and altitude vectors
        // This is a temporary fix, even if not ideal, we verify all generated values and replace them if not suitable

        for(int i = 0; i < n_rays_thread ; i++)
        {
            double phi = -( M_azimuthAngles[index_azimuth(i)] ) + M_PI*0.5 ; // recover spherical coordinate from azimuth angle
            double theta = M_PI*0.5 - M_altitudeAngles[index_altitude(i)]; // recover spherical coordinate from altitude
            random_directions(i,0)=math::sin(theta)*math::cos(phi);
            random_directions(i,1)=math::sin(theta)*math::sin(phi);
            random_directions(i,2)=math::cos(theta);
        }
        return random_directions;
    }

    Eigen::MatrixXd get_random_points(matrix_node_type const& element_points, int n_rays_thread, Eigen::Rand::UniformRealGen<double>& unif_gen_real1, Eigen::Rand::UniformRealGen<double>& unif_gen_real2, Eigen::Rand::P8_mt19937_64& urng_real1, Eigen::Rand::P8_mt19937_64& urng_real2)
    {
        int dimension;
        dimension = column(element_points, 0).size();

        if(dimension==3)
        {
            Eigen::MatrixXd random_points(n_rays_thread, 3);
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
                Eigen::VectorXd S = unif_gen_real1.generate<Eigen::VectorXd>(n_rays_thread,1,urng_real1);
                Eigen::VectorXd T = unif_gen_real2.generate<Eigen::VectorXd>(n_rays_thread,1,urng_real2);
                // If the point is on the left of the diagonal, keep it, else take the symmetric one
                for (int j=0 ; j<n_rays_thread ; j++)
                {
                    bool in_triangle = (S(j) + T(j) <= 1);
                    if(in_triangle)
                        p = p1 + S(j) * u + T(j) * v;  
                    else 
                        p= p1 + (1 - S(j)) * u + (1 - T(j)) * v;

                    if (isOnSurface(p,p1,p2,p3))
                        random_points.row(j) = p;
                    else
                    {
                        throw std::logic_error("Point not on triangle, but it must be");
                        random_points.row(j) = p1;
                    }  
                }
                return random_points;
            }
        }
        else
        {
            Eigen::MatrixXd p1(n_rays_thread, 3);
            p1.setZero();
            for(int i=0; i < n_rays_thread; i++)
            {
                for(int j=0;j<dimension;j++)
                {
                    p1(i,j)=column(element_points, 0)[j];
                }
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
        Eigen::MatrixXd random_directions(M_Nrays/M_Nthreads,dim);

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

                        Eigen::Rand::UniformIntGen<int> unif_gen_azi(0,M_azimuthSize-1);
                        Eigen::Rand::UniformIntGen<int> unif_gen_alti(0,M_altitudeSize-1);
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

                        Eigen::VectorXi index_altitude(n_rays_thread);
                        Eigen::VectorXi index_azimuth(n_rays_thread);

                        Eigen::MatrixXd random_directions(n_rays_thread,3);
                        Eigen::MatrixXd random_origins(n_rays_thread,3);

                        random_origins = get_random_points(el.second.vertices(),n_rays_thread, unif_gen_real1, unif_gen_real2, urng_real1, urng_real2);
                        random_directions = get_random_directions(index_azimuth, index_altitude, n_rays_thread, unif_gen_azi, unif_gen_alti, urng_azi, urng_alti);
                        for(int j=0;j<n_rays_thread;j++)
                        {
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
                                    rand_dir(i) = random_directions(j,i);
                                    origin(i) = random_origins(j,i);
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
                                // the index values were verified, they are correct           
                                SM_table(index_azimuth(j),index_altitude(j))++;
                                Angle_table(index_azimuth(j),index_altitude(j))++;
                            }    
                            else
                            {
                                Angle_table(index_azimuth(j),index_altitude(j))++;
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
        
        std::string matrix_filename = shadingMaskFolder+"/SM_Matrix_"+building_name+"_"+marker_name+"EIGENRANDVECTORIZED.csv";
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