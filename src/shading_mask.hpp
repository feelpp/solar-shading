#include <future>
#include <feel/feelmesh/bvh.hpp>
#include <unordered_map>



auto GetListNameObjects(std::string ChName)
{
    std::ifstream FICH(ChName);
    std::string lineA;
    std::vector<std::string> ListObjects;
    while (getline(FICH, lineA)) { ListObjects.push_back(lineA); } 
    FICH.close();
    return ListObjects;
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

    using mesh_entity_type = typename std::conditional<MeshType::nDim==MeshType::nRealDim,
                                            entity_range_t<elements_reference_wrapper_t<MeshType>>,
                                            entity_range_t<faces_reference_wrapper_t<MeshType>> >::type;

    using wrapper_type = elements_reference_wrapper_t<MeshType>;
                                        
    typedef typename matrix_node<double>::type matrix_node_type;

public:
    using value_type = double;

    ShadingMask(mesh_ptrtype mesh, nl::json const& specs, int intervalsAzimuth=72, int intervalsAltitude=10 );
    
    // Create the random number generators
    void makeRandomNumberGeneratorsSeed();

     // Create and store the directions of the M_Nrays
    void makeCreateM_NraysMatrix(int intervalsAzimuth, int intervalsAltitude);

    // For each building, save the surface mesh and build the corresponding BVH tree for ray search
    void loadMeshDataSubPart1(mesh_ptrtype mesh,int numOp);

    void loadMeshDataSubPart2(mesh_ptrtype mesh,int numOp);

    void loadMeshDataSubPart3(mesh_ptrtype mesh);

    void loadMeshData(mesh_ptrtype mesh, nl::json const& specs);
    
    // Subdivide the azimuth angles [0,360]° and altitude angles [0,90]° in subsets for easier computation of the shading masks
    void fixAzimuthAltitudeDiscretization(int intervalsAzimuth=72, int intervalsAltitude=10);

    // Choose a random pair of indices in the discretized azimuth and altitude vectors
    void getRandomDirectionSM(std::vector<double> &random_direction, std::mt19937 & M_gen, std::mt19937 & M_gen2, int& index_azimuth, int& index_altitude);

    // Get a random point from the surface of the triangle
    Eigen::VectorXd get_random_point(matrix_node_type const& element_points);

    // 3D case
    // Compute the sum of the areas of three subtriangles
    double elementArea(Eigen::VectorXd const& point,Eigen::VectorXd const& el_p1,Eigen::VectorXd const& el_p2,Eigen::VectorXd const& el_p3);

    // 3D case
    // Compare the area of the 2d simplex V1V2V3 (as sum of 3 subtriangles V_iV_jB) and the 2d triangle
    // created by the intersection P of the ray with the plane of  V1V2V3 (as sum of V_iV_jP)
    bool isOnSurface(Eigen::VectorXd const &point,Eigen::VectorXd const &el_p1,Eigen::VectorXd const &el_p2,Eigen::VectorXd const &el_p3);

    void computeMasksSubPart0(int numOp);

    void computeMasksSubPart1(int numOp);

    void computeMasksSubPart2(int numOp);

    // Compute shading masks for the buildings in the json file
    void computeMasksMaster();

    // Compute shading masks for the buildings in the json file
    void computeMasks();

    // Save Compute shading masks
    void computeSaveMasks(std::vector<double> SM_tables);

    // Compute shading masks for one building only
    void computeMasksOneBuilding(std::string building_name);

    // Save the shading mask table to a CSV file
    void saveShadingMask(std::string building_name, std::string marker_name, const Eigen::Ref<const Eigen::MatrixXd>& M);

    // Save the shading mask table metadata to json
    void saveMetadata();

    

private:
    
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

    int matrixSize;

    nl::json j_;

    nl::json M_metadataJson;

    std::random_device M_rd;
    std::random_device M_rd2;
    std::mt19937 M_gen;
    std::mt19937 M_gen2;
};
} // namespace Feel

#include<shading_mask_init.hpp>
#include<shading_mask_geometry.hpp>
#include<shading_mask_compute.hpp>
#include<shading_mask_save.hpp>