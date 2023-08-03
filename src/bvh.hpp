#ifndef FEELPP_BVH_HPP
#define FEELPP_BVH_HPP

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feel.hpp>

namespace Feel
{
    class BVHRay
    {
    public:
        using vec_t = Eigen::VectorXd;
        BVHRay(const vec_t &orig, const vec_t &dir);
        BVHRay(const BVHRay &r);
        BVHRay();
        vec_t origin, dir;
    };

    inline constexpr float robust_const(int n);

    template<int nDim>
    class BVHTree
    {
        using self_type = BVHTree<nDim>;
        typedef Simplex<nDim,1> convex_type;
        typedef Mesh<convex_type> mesh_type;
        typedef typename mesh_type::trace_mesh_type trace_mesh_type;
        typedef typename mesh_type::trace_mesh_ptrtype trace_mesh_ptrtype;
        typedef typename matrix_node<double>::type matrix_node_type;

        public:
        struct BVHPrimitiveInfo
        {
            BVHPrimitiveInfo(int primitiveNumber, Eigen::VectorXd bounds_min, Eigen::VectorXd bounds_max);
            int M_primitiveNumber;
            Eigen::VectorXd M_bound_min;
            Eigen::VectorXd M_bound_max;
            Eigen::VectorXd M_centroid;
        };

        void buildPrimitivesInfo(trace_mesh_ptrtype const& mesh);
        std::vector<BVHPrimitiveInfo> M_primitiveInfo;
        trace_mesh_ptrtype M_mesh;
        thread_local static inline std::vector<int> M_intersected_leaf;
        thread_local static inline std::vector<double> M_lengths;
        std::vector<int> orderedPrims;

        class BVHNode
        {
            public:
            void buildLeaf(BVHNode * current_parent, int first, int n, Eigen::VectorXd bounds_min, Eigen::VectorXd bounds_max);
            void buildInternalNode(BVHNode * current_parent, int splitaxisIn, BVHNode *child0, BVHNode *child1);
            Eigen::VectorXd newBoundsMin(Eigen::VectorXd bounds1,Eigen::VectorXd bounds2);
            Eigen::VectorXd newBoundsMax(Eigen::VectorXd bounds1,Eigen::VectorXd bounds2);
            bool isLeaf();
            BVHNode * nearChild(BVHRay const& ray);
            BVHTree::BVHNode * otherChild(BVHTree::BVHNode * parent);
            bool checkIntersection(BVHRay const& rayon);
            std::pair<bool,double> checkIntersectionWithSegment(matrix_node_type const& nodes, BVHRay const& ray);
            std::pair<bool,double> checkIntersectionWithTriangle( BVHRay const& ray, trace_mesh_ptrtype mesh, std::vector<BVHPrimitiveInfo>& primitiveInfo);
            std::pair<bool,double> checkLeafIntersection(BVHRay const& rayon, trace_mesh_ptrtype mesh, std::vector<BVHPrimitiveInfo>& primitiveInfo);
            BVHNode *children[2];
            BVHNode *parent;
            int splitaxis, nPrimitives,firstPrimOffset;
            Eigen::VectorXd M_bounds_min,M_bounds_max,M_centroid;
        };
        BVHNode *  M_root_tree;

        BVHNode * getRootNode();
        BVHNode * buildRootTree();
        BVHNode * recursiveBuild(BVHNode * current_parent, int cut_dimension, int start_index_primitive, int end_index_primitive, std::vector<int> &orderedPrims);
        int raySearch( Feel::BVHRay const& rayon,std::string s);
        void traverse_stackless(BVHNode * tree, BVHRay const& rayon);
    };
}

#endif // FEELPP_BVH_HPP
