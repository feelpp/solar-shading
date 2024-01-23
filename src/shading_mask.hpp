#include <future>
#include <feel/feelmesh/bvh.hpp>
#include <unordered_map>

#include <functional>



auto GetListNameObjects(std::string ChName)
{
    std::ifstream FICH(ChName);
    std::string lineA;
    std::vector<std::string> ListObjects;
    while (getline(FICH, lineA)) { ListObjects.push_back(lineA); } 
    FICH.close();
    return ListObjects;
}

//BEGIN::CONSOLE TOOLS FOR LINUX
void CONSOLE_ClearScreen()                     { printf("\033[2J"); }
void CONSOLE_SaveCursorPosition()              { printf("\033[s"); }
void CONSOLE_RestoreCursorPosition()           { printf("\033[u");}
void CONSOLE_SetCursorPosition(int x, int y)   { printf("\033[%d;%dH",y+1,x+1); }
void CONSOLE_GetCursorPosition(int* x, int* y) { printf("\033[6n");  int err=scanf("\033[%d;%dR", x, y); }
void CONSOLE_CursorBlinkingOnOff()             { printf("\033[?12l"); }
void CONSOLE_CursorHidden(bool q)              { q ? printf("\e[?25l"):printf("\e[?25h"); }

void CONSOLE_PRINT_PERCENTAGE(int x, int y, int k, int nb)
{
	float v = ((float(k)) / float(nb)) * 100.0;
	CONSOLE_SetCursorPosition(x,y); printf(" %3.1f %%", v);
}

void CONSOLE_PRINT_PERCENTAGE(int k, int nb)
{
	float v = ((float(k)) / float(nb)) * 100.0;
	CONSOLE_RestoreCursorPosition(); printf(" %3.1f %%", v);
}

void CONSOLE_Color(int n)
{
    switch(n) {
        case 0:
            printf("\e]P0000000"); ///black
        break;
        case 1:
            printf("\e]P1D75F5F"); //darkred
        break;
        case 2:
            printf("\e]P287AF5F"); //darkgreen
        break;
        case 3:
            printf("\e]P3D7AF87"); //brown
        break;
        case 4:
            printf("\e]P48787AF"); //darkblue
        break;
        case 5:
            printf("\e]P5BD53A5"); //darkmagenta
        break;
        case 6:
            printf("\e]P65FAFAF"); //darkcyan
        break;
        case 7:
            printf("\e]P7E5E5E5"); //lightgrey
        break;
        case 8:
            printf("\e]P82B2B2B"); //darkgrey
        break;
        case 9:
            printf("\e]P9E33636"); //red
        break;
        case 10:
            printf("\e]PA98E34D"); //green
        break;
        case 11:
            printf("\e]PBFFD75F"); //yellow
        break;
        case 12:
            printf("\e]PC7373C9"); //blue
        break;
        case 13:
            printf("\e]PDD633B2"); //magenta
        break;
        case 14:
            printf("\e]PE44C9C9"); //cyan
        break;
        case 15:
            printf("\e]PFFFFFFF"); //white
        break;
    }
}
//1-CONSOLE_SaveCursorPosition();
//2-CONSOLE_CursorHidden(true); 
//3-LOOP   CONSOLE_PRINT_PERCENTAGE(i+1,M_listFaceMarkers.size());
//4-CONSOLE_CursorHidden(false); 
//END::CONSOLE TOOLS



//BEGIN::FISRT PART OF TASKDISPACH MULTITHREAD
class MyTaskDispach
{
    public:
        int numTypeTh;
        int nbTh;
        bool QSave;
        std::string FileName;
        void init(int numType,int nbThread,bool QsaveInfo);
        void setFileName(std::string s);
        template<class Function>
            Function run(Function myFunc);
        MyTaskDispach(void);
};


MyTaskDispach::MyTaskDispach() { 
    numTypeTh=2; 
    nbTh=6;
    QSave=false;
    FileName="TestDispach";
}

void MyTaskDispach::init(int numType,int nbThread,bool QsaveInfo)
{
    numTypeTh=numType; nbTh=nbThread; QSave=QsaveInfo;
}

void MyTaskDispach::setFileName(std::string s)
{
    FileName=s;
}

template<class Function>
Function MyTaskDispach::run(Function myFunc)
{
    bool QInfo=false; 
    if (numTypeTh==1) //with std::async
    {
        std::vector< std::future< bool > > futures;
        for(int k= 0; k < nbTh; ++k){ 
            auto const& idk = k;
            if (QInfo) { std::cout<<"Call num Thread futures="<<k<<"\n"; }
            futures.emplace_back(std::async(std::launch::async,myFunc,idk));
        }
        for( auto& r : futures){ auto a =  r.get(); }
        if (QInfo) { std::cout<<"\n"; }
    }

    if (numTypeTh==2) //With Specx
    {
        SpRuntime runtime(nbTh);  
        nbTh= runtime.getNbThreads();
        int iValue=0;
        for(int k= 0; k < nbTh; ++k)
        { 
            auto const& idk = k;
            runtime.task(SpRead(idk),myFunc).setTaskName("Op("+std::to_string(k)+")");
            usleep(1);
            std::atomic_int counter(0);
        }
        runtime.waitAllTasks();
        runtime.stopAllThreads();
        if (QSave)
        {
            runtime.generateDot(FileName+".dot", true);
            runtime.generateTrace(FileName+".svg");   
        }
        if (QInfo) { std::cout<<"\n"; }
    }
    return myFunc;
}

//END::FISRT PART OF TASKDISPACH MULTITHREAD

//SECOND PART OF TASKDISPACH MULTITHREAD
//...


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
    bool QSaveTypeThreadDotON;
    bool QSaveControlFiles;
    int  numTypeThread;


    ShadingMask(int num,mesh_ptrtype mesh, nl::json const& specs, int intervalsAzimuth=72, int intervalsAltitude=10 );

    auto commonComputePartCTRL(int NumOption,matrix_node_type const& element_points,int n_rays_thread, int id_thread);

    bool computePartMarker(std::vector<std::string> marker_list_thread, int id_thread, int start_index);



    // Create the random number generators
    void makeRandomNumberGeneratorsSeed(bool QCTRL_SAVE_SEED,bool QCTRL_LOAD_SEED);

     // Create and store the directions of the M_Nrays
    void makeCreateM_NraysMatrix(int intervalsAzimuth, int intervalsAltitude);

    // For each building, save the surface mesh and build the corresponding BVH tree for ray search
    void loadMeshData(mesh_ptrtype mesh, nl::json const& specs);
        void loadMeshDataSubPartList          (mesh_ptrtype mesh);
        void loadMeshDataSubPartVolumes       (mesh_ptrtype mesh,int numOp);
        void loadMeshDataSubPartSurfacesFaces (mesh_ptrtype mesh,int numOp);
        void loadMeshDataSubPartMarkers       (mesh_ptrtype mesh);

    // Subdivide the azimuth angles [0,360]° and altitude angles [0,90]° in subsets for easier computation of the shading masks
    void fixAzimuthAltitudeDiscretization(int intervalsAzimuth=72, int intervalsAltitude=10);

    // Choose a random pair of indices in the discretized azimuth and altitude vectors
    void getRandomDirectionSM(std::vector<double> &random_direction, std::mt19937 & M_gen, std::mt19937 & M_gen2, int& index_azimuth, int& index_altitude);

    // Get a random point from the surface of the triangle
    Eigen::VectorXd get_random_point(matrix_node_type const& element_points);

    Eigen::VectorXd get_element_normal (Eigen::VectorXd p1,Eigen::VectorXd p2,Eigen::VectorXd p3);
   

    // 3D case
    // Compute the sum of the areas of three subtriangles
    double elementArea(Eigen::VectorXd const& point,Eigen::VectorXd const& el_p1,Eigen::VectorXd const& el_p2,Eigen::VectorXd const& el_p3);

    // 3D case
    // Compare the area of the 2d simplex V1V2V3 (as sum of 3 subtriangles V_iV_jB) and the 2d triangle
    // created by the intersection P of the ray with the plane of  V1V2V3 (as sum of V_iV_jP)
    bool isOnSurface(Eigen::VectorXd const &point,Eigen::VectorXd const &el_p1,Eigen::VectorXd const &el_p2,Eigen::VectorXd const &el_p3);


    // Compute shading masks for the buildings in the json file
    void computeMasksMaster();
        void computeMasksSubPartList();
        void computeMasksSubPartSurfaceVolumes(int numOp);
        void computeMasksSubPartMarkersCTRL();

    void computeThread(Eigen::MatrixXd SM_table_marker,Eigen::MatrixXd Angle_table_marker,matrix_node_type const& element_points);

    // Save Compute shading masks
    void computeSaveMasks(std::vector<double> SM_tables);

    // Compute shading masks for one building only
    void computeMasksOneBuildingCTRL(std::string building_name);

    //void computeMasksOneBuildingOld(std::string building_name);

    // Save the shading mask table to a CSV file
    void saveShadingMask(std::string prefix_name,std::string building_name, std::string marker_name, const Eigen::Ref<const Eigen::MatrixXd>& M);

    // Save the shading mask table metadata to json
    void saveMetadata(std::string name);

    void saveMetadataInfoPart();

    void testComparisonAllMasksValidation();
    bool testComparisonMaskValidationLevel1(std::string matrix_filename_NEW,std::string matrix_filename_CTRL);
    bool testComparisonMaskValidationLevel2(std::string matrix_filename_NEW,std::string matrix_filename_CTRL);


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
    int dim=3;

    nl::json j_,M_metadataJson;

    std::random_device M_rd,M_rd2;
    std::mt19937 M_gen,M_gen2;

    std::vector<double> SM_tables_Alpha;
    std::vector<double> Angle_tables_Alpha;

    std::time_t beginning_time;

    
};
} // namespace Feel

#include<shading_mask_init.hpp>
#include<shading_mask_geometry.hpp>
#include<shading_mask_compute.hpp>
#include<shading_mask_save.hpp>
#include<shading_mask_control.hpp>






