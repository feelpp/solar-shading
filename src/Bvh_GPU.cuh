#ifndef MY_KERNELS_CUH
#define MY_KERNELS_CUH

#include "bvh.hpp"
#include <vector>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <device_atomic_functions.h>
#include <cuda.h>

// Error checking macro
#define CHECK_CUDA_ERRORS(call) checkCudaErrors((call), __FILE__, __LINE__)

// Functions and Class declarations
__host__ void checkCudaErrors(cudaError_t error, const char* file, int line);

struct GPURay
{
    float origin[3];
    float dir[3];

    GPURay();
    GPURay(const float o[3], const float d[3]);
    GPURay(const Feel::BVHRay R);
};

struct GPUNode
{
    GPUNode* parent;
    GPUNode* leftchild;
    GPUNode* rightchild;
    int nPrimitives;
    int firstPrimOffset;
    int splitAxis;
    float bounds_min[3];
    float bounds_max[3];
    float centroid[3];

    GPUNode(const Feel::BVHTree<3>::BVHNode& bvhNode);

    __host__ void buildLeaf(int first, int n, const float bmin[3], const float bmax[3], const float cent[3]);
    __host__ void buildInternalNode(int splitaxisIn, GPUNode* child0, GPUNode* child1);
    __host__ GPUNode* DeepCopyToGPU();
    __device__ GPUNode* nearChild(GPURay const& ray);
    __device__ GPUNode* otherChild(GPUNode* parent);
    __device__ bool checkIntersection(GPURay const& rayon);
    __device__ bool isLeaf();
    __device__ int getfirstPrimOffset();
};

struct GPUBVH
{
    GPUNode* M_root_gpu_tree;

    __host__ GPUNode* buildRootTree(Feel::BVHTree<3>* tree);
    __device__ void GPU_traverse_stackless(GPUNode* tree, GPURay const& ray, int* results, int& result_count);
    __global__ void GPU_traverse_kernel(GPUNode* tree, GPURay const* rays, int* results, int numRays);
    __host__ std::vector<int> GPUraySearch(std::vector<Feel::BVHRay> const& rays, const Feel::BVHTree<3>* tree);
};

std::vector<int> GPUraySearchWrapper(std::vector<Feel::BVHRay> const& rays, const Feel::BVHTree<3> * tree);

#endif
