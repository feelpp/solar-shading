#include <feel/feelpython/pybind11/pybind11.h>
#include <feel/feelpython/pybind11/stl.h>
#include <feel/feelpython/pybind11/json.h>

#include <shading_mask.hpp>


namespace py = pybind11;

using namespace Feel;

template<typename MeshType>
void defSM(py::module &m)
{
    using namespace Feel;
    using sm_t = ShadingMask<MeshType> ;

    using mesh_type = MeshType;
    typedef typename MeshType::ptrtype mesh_ptrtype;

    std::string pyclass_name = std::string("ShadingMask_") + std::to_string(mesh_type::nDim) + std::string("D") + std::to_string(mesh_type::nRealDim) + std::string("D") ;
    py::class_<sm_t,std::shared_ptr<sm_t>>(m,pyclass_name.c_str())
        .def(py::init<mesh_ptrtype, nl::json const&, int , int  >(),
             py::arg("mesh"),
             py::arg("json"),
             py::arg("intervalsAzimuth")=72,
             py::arg("intervalsAltitude")=10,
             "Initialize the shading mask class"
             )
        .def("computeMasks",&sm_t::computeMasks, "Compute the shading masks for all the specified buildings");
}

PYBIND11_MODULE(_shadingmask, m) {

    using namespace Feel;
    
    defSM<Mesh<Simplex<2,1,3>>>(m);
    defSM<Mesh<Simplex<3,1,3>>>(m);

}