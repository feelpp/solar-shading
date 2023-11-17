#include <future>
#include <feel/feelmesh/bvh.hpp>
#include <unordered_map>


#include "SpDataAccessMode.hpp"
#include "Utils/SpUtils.hpp"

#include "Task/SpTask.hpp"
#include "Legacy/SpRuntime.hpp"
#include "Utils/SpTimer.hpp"
#include "Utils/small_vector.hpp"
#include "Utils/SpConsumerThread.hpp"


#include <ncurses.h>


void GetSpecxPreprocessingParameters(int NbObjects,int NbTasks,int & NbLoop, int & NbTh, int & NbCoor,bool QViewInfo)
{
    if (QViewInfo) { std::cout<<"[SPECX INFO] : SpecxSaveNbThreadDesired="<<NbTasks; }
        NbTasks=std::min(NbTasks,NbObjects);
        int NbTasksToSubmit = std::min(NbTasks,std::min(NbObjects,SpUtils::DefaultNumThreads()));
        NbLoop=std::max(1,NbObjects/NbTasks);
        NbCoor=round((float(NbObjects)/float(NbTasks)-float(NbObjects/NbTasks))*float(NbTasks));
        NbTh=NbTasks; 
    if (QViewInfo) { std::cout<<" Nb Objects="<<NbObjects<<" Nb Loops="<<NbLoop<<" Coor="<<NbCoor<<" NbThreads="<<NbTh<<" used\n"; }
}


auto GetListNameObjects(std::string ChName)
{
    std::ifstream FICH(ChName);
    std::string lineA;
    std::vector<std::string> ListObjects;
    while (getline(FICH, lineA)) { ListObjects.push_back(lineA); } 
    FICH.close();
    return ListObjects;
}


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

    //ShadingMask(mesh_ptrtype mesh, nl::json const& specs, int intervalsAzimuth=72, int intervalsAltitude=10 );

    ShadingMask(int TestNbThread,int TestNbRayon,mesh_ptrtype mesh, nl::json const& specs,int intervalsAzimuth=72, int intervalsAltitude=10);
    
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

    // Compute shading masks for the buildings in the json file
    void computeMasks();

    // Compute shading masks for one building only
    void computeMasksOneBuilding(std::string building_name);

    // Save the shading mask table to a CSV file
    void saveShadingMask(std::string building_name, std::string marker_name, const Eigen::Ref<const Eigen::MatrixXd>& M);

    // Save the shading mask table metadata to json
    void saveMetadata();

    // Test if it the same
    bool testShadingMaskComparisonLevel1(std::string shadingMaskFolder,std::string building_name,std::string marker_name);

    bool testShadingMaskComparisonLevel2(std::string shadingMaskFolder,std::string building_name,std::string marker_name);


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

    nl::json j_;

    nl::json M_metadataJson;

    std::random_device M_rd;
    std::random_device M_rd2;
    std::mt19937 M_gen;
    std::mt19937 M_gen2;

    int  NbObjects; 
    int  SpecxNbThreadDesired;
    int  SpecxSaveNbThreadDesired;
    bool QSaveSpecxGenerateReports;
    bool QSpecxLockConfigWaitSave;
    bool QSaveWithSpecx;
    bool QViewInfoSpecx;
    bool QCTRL_SAVE_NORMAL;
    bool QCTRL_DATA;
    bool QCTRL_SAVE_SEED;
    bool QCTRL_LOAD_SEED;
    bool QMAKE_WITH_SPECX;

    int  SpecxActivateSubNumberPart;

    std::time_t beginning_time;
};
} // namespace Feel

#include<Specx_Shading_mask_init.hpp>
#include<Specx_Shading_mask_geometry.hpp>
#include<Specx_Shading_mask_compute.hpp>
#include<Specx_Shading_mask_save.hpp>
