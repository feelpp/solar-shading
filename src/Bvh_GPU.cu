#include "bvh.hpp"
#include "Bvh_GPU.cuh"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <device_atomic_functions.h>
#include <cuda.h>
#include <vector>
#include <algorithm>
#include <cfloat>


__host__ void checkCudaErrors(cudaError_t error, const char* file, int line)
{
    if (error != cudaSuccess)
    {
        printf("\nCUDA Error at %s:%d: %s\n", file, line, cudaGetErrorString(error));
        exit(-1);
    }
}
#define CHECK_CUDA_ERRORS(call) checkCudaErrors((call), __FILE__, __LINE__)

struct GPURay
{
    float origin[3];
    float dir[3];

    GPURay() : origin{0.0f, 0.0f, 0.0f}, dir{0.0f, 0.0f, 0.0f} {}

    GPURay(const float o[3], const float d[3]) : origin{o[0], o[1], o[2]}, dir{d[0], d[1], d[2]} {}

    GPURay( const Feel::BVHRay R)
    {
        for (int i=0;i<3;i++)
        {
            origin[i]=R.origin[i];
            dir[i]=R.dir[i];
        }
    }
};

struct GPUNode
{
    GPUNode* parent;      // using pointers may cause enormous problems, since they are instanciated on the CPU, 
    GPUNode* leftchild;   // the pointers will still refer to the CPU memory, and not the GPU memory, so we could 
    GPUNode* rightchild;  // either use an array of indices ( for example giving the nodes indices during constructions),
    int nPrimitives;      // But since everythin is already implemented, we will use pointers, and implement a deep copy
    int firstPrimOffset;  // function that will copy the whole tree from the CPU to the GPU, and adjust the pointers
    int splitAxis;
    float bounds_min[3]; // one float per dimension
    float bounds_max[3];
    float centroid[3];

    GPUNode(const Feel::BVHTree<3>::BVHNode& bvhNode) 
    {
        parent = nullptr; // Parent needs to be set by whoever creates this node

        if (bvhNode.isLeaf()) 
        {
            buildLeaf(bvhNode.firstPrimOffset, 
                    bvhNode.nPrimitives, 
                    bvhNode.M_bounds_min.data(), 
                    bvhNode.M_bounds_max.data(), 
                    bvhNode.M_centroid.data());
        } 
        else 
        {
            this->splitaxis = bvhNode.splitaxis;
            if (bvhNode.children[0] != nullptr)
            {
                this->leftchild = new GPUNode(*(bvhNode.children[0]));
                this->leftchild->parent = this;
            }
            if (bvhNode.children[1] != nullptr) 
            {
                this->rightchild = new GPUNode(*(bvhNode.children[1]));
                this->rightchild->parent = this;
            }
            buildInternalNode(this->splitaxis,
                            this->leftchild,
                            this->rightchild);
        }
    }

    __host__ void buildLeaf(int first, int n, const float bmin[3], const float bmax[3], const float cent[3]) {
        firstPrimOffset = first;
        nPrimitives = n;
        leftchild = nullptr;
        rightchild = nullptr;
        for (int i=0; i<3; i++) {
            bounds_min[i] = bmin[i];
            bounds_max[i] = bmax[i];
            centroid[i] = cent[i];
        }
    }

    __host__ void buildInternalNode(int splitaxisIn, GPUNode* child0, GPUNode* child1) {
        leftchild = child0;
        rightchild = child1;
        splitAxis = splitaxisIn; 
        nPrimitives = 0;
        firstPrimOffset = -1;
        for (int i=0; i<3; i++) {
            bounds_min[i] = std::min(child0->bounds_min[i], child1->bounds_min[i]);
            bounds_max[i] = std::max(child0->bounds_max[i], child1->bounds_max[i]);
            centroid[i] = (bounds_min[i] + bounds_max[i]) / 2.0f;
        }
    }

        
    // Deep Copy recursive function
    __host__ GPUNode* DeepCopyToGPU() 
    {
        if (this == NULL)
        {
            return nullptr;
        }

        // Allocate memory for this node on the GPU
        GPUNode *node_copy;
        cudaMalloc(&node_copy, sizeof(GPUNode));

        // Copy the node data to the GPU
        cudaMemcpy(node_copy, this, sizeof(GPUNode), cudaMemcpyHostToDevice);

        // Make recursive calls to copy children
        GPUNode *d_leftChild = nullptr;
        GPUNode *d_rightChild = nullptr;

        if (this->leftchild)
        {
            d_leftChild = this->leftchild->DeepCopyToGPU();
            cudaMemcpy(&(node_copy->leftchild), &d_leftChild, sizeof(GPUNode*), cudaMemcpyHostToDevice);
        }
        if (this->rightchild)
        {
            d_rightChild = this->rightchild->DeepCopyToGPU();
            cudaMemcpy(&(node_copy->rightchild), &d_rightChild, sizeof(GPUNode*), cudaMemcpyHostToDevice);
        }

        return node_copy;
    }



    // some methods for the traversal
    __device__ GPUNode * nearChild(GPURay const& ray)
    {
                
        if(ray.dir[this->splitaxis]>0)
            return this->child0;
        else
            return this->child1;

    }
            
    __device__ GPUNode * otherChild(GPUNode * parent)
    {
        if (this==parent->child0)
            return parent->child1;
         else
            return parent->child0;
    }

    __device__ bool checkIntersection(GPURay const& rayon)
    {
        float tmin = 0.0;
        float tmax = FLT_MAX;

        for(int i=0; i<nDim; i++)
        {                    
            float ratio = 1.0/(rayon.dir[i]+2*FLT_MIN);
            float t1 = (this->bounds_min[i]-rayon.origin[i]) * ratio;
            float t2 = (this->bounds_max[i]-rayon.origin[i]) * ratio;
            if (t1 > t2) 
            {
                float tTemp = t1;
                t1 = t2;
                t2 = tTemp;
            }
            if ( t1 > tmin)
                tmin = t1;
            if (t2 > tmax)
                tmax = t2;
            if (tmin > tmax)
                return false;
        }
        
        return true;         
    }

    __device__ bool isLeaf()
    {
        return (this->nPrimitives>0);
    }

    __device__ int getfirstPrimOffset()
    {
        return this->firstPrimOffset;
    }
};
        
struct GPUBVH
{
    GPUNode * M_root_gpu_tree;

    __host__ GPUNode * buildRootTree(Feel::BVHTree<3> * tree)
    {
        // Copy the necessary data from the BVHTree into the GPUBVH
        Feel::BHVTree::BVHNode * root_cpu_tree = tree->getRootNode();
        GPUNode * root_gpu_node = new GPUNode(root_cpu_tree);
        return root_gpu_node;
    }
    
    // this function only returns the first prim offset, we can link the first prim offset (stored in results)
    // to the ray which intersected it by using its index in the results array. If no intersection is found, the
    // array will display -2.
    __device__ void GPU_traverse_stackless(GPUNode * tree, GPURay const& ray, int * results, int & result_count)
    {
        auto current_node = tree -> nearChild(ray);
        char state = 'P'; 

        result_count = 0;       

        while (true)
        {
            switch (state)
            {
                case 'C':
                    if (current_node == M_root_gpu_tree) return;

                    if (current_node == current_node->parent->nearChild(ray))
                    {
                        current_node = current_node->otherChild(current_node->parent);
                        state = 'S'; // from Sibling
                    }
                    else 
                    {
                        current_node = current_node->parent;
                        state = 'C'; // the current node has been accessed from its sibling
                    }
                    break;

                case 'S': // the node is being traversed from its sibling
                    if (current_node->checkintersection(ray)==false) // go back to parent
                    {
                            current_node = current_node->parent;
                        state = 'C';
                    }
                    else if (current_node->isLeaf())
                    {
                        // either perform the Ray/Primitive intersection test here, or return 
                        // the primitive indices and let the CPU do the intersection tests
                        // if we pushback the firstPrimOffset, then we can use it to retrieve the primitive
                        // in the M_primitiveinfo array
                        int index = atomicAdd(&result_count, 1);            // Increment the result_count atomically and get the previous value as the index. This way, even if 
                        results[index] = current_node->getfirstPrimOffset();// multiple threads try to write to the same index, we will not lose any results since they'll be queued
                        current_node = current_node->parent;
                        state = 'C';
                    }
                    else
                    {
                        current_node = current_node->parent;
                        state = 'P';
                    }
                    break;
                
                case 'P':
                    if (current_node->checkIntersection(ray)==false)
                    {
                        current_node = current_node->otherChild(current_node->parent);
                        state = 'S';
                    }
                    else if (current_node->isLeaf())
                    {
                        int index = atomicAdd(&result_count, 1); // Increment the result_count atomically and get the previous value as the index
                        results[index] = current_node->getfirstPrimOffset();    
                        current_node = current_node->otherChild(current_node->parent);
                        state = 'S';
                    }
                    else
                    {
                        current_node = current_node->nearChild(ray);
                        state = 'P';
                    }
                    break;
                
                default:

                    break;
            }
        }
    }

    __global__ void GPU_traverse_kernel(GPUNode* tree, GPURay const* rays, int* results, int numRays)
    {
        int index = threadIdx.x + blockIdx.x * blockDim.x;
        int N = 10 ; // this is the max number of intersections we can have per ray
        if (index < numRays)
        {
            int thread_results[N]; // Local array specific to each thread
            int result_count = 0; // Local variable specific to each thread
            GPU_traverse_stackless(tree, rays[index], thread_results, result_count);

            for (int i = 0; i < result_count; i++)
            {
                results[index * 10 + i] = thread_results[i]; // Store results in the global results array
            }
        }
    }

    __host__ std::vector<int> GPUraySearch(std::vector<Feel::BVHRay> const& rays, const Feel::BVHTree<3> * tree)
    {
        int totalRays = rays.size();
        int numDevices;
        CHECK_CUDA_ERRORS(cudaGetDeviceCount(&numDevices));
        int raysPerDevice = totalRays / numDevices; // Assuming totalRays is divisible by numDevices here.
        std::vector<double> lengths; // no distances are computed on the GPU
        std::vector<GPURay> rayons;
        std::vector<int> results(totalRays, -2);

        // convert the BVHRays to GPURays here
        for (int i = 0; i < totalRays; i++)
        {
            rayons.push_back(GPURay(rays[i]));
        }

        // Get all informations on the devices needed to perform the ray search
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);// Assumes all devices are identical
        int maxThreadsDim = prop.maxThreadsDim[0];  
        int maxGridSize = prop.maxGridSize[0];
        int maxThreadsPerBlock = prop.maxThreadsPerBlock;
        size_t totalGlobalMem = prop.totalGlobalMem;// Total global memory (in bytes) ===> can be used to check whether the tree fits in the GPU memory
        int threadsPerBlock = std::min(totalRays, maxThreadsPerBlock);
        int blocks = (totalRays + threadsPerBlock - 1) / threadsPerBlock; 
        int blocksPerGrid = (raysPerDevice + threadsPerBlock - 1) / threadsPerBlock; // Round up division

        M_root_gpu_tree = buildRootTree(tree);


        cudaStream_t stream[numDevices];
        GPUNode *d_tree[numDevices];
        GPURay *d_rays[numDevices];
        int *d_results[numDevices];

        for (int i = 0; i < numDevices; i++) {
            CHECK_CUDA_ERRORS(cudaSetDevice(i));
            CHECK_CUDA_ERRORS(cudaStreamCreate(&stream[i]));

            // Allocate device memory for rays and copy from host to device
            CHECK_CUDA_ERRORS(cudaMalloc(&d_rays[i], sizeof(GPURay) * raysPerDevice));
            CHECK_CUDA_ERRORS(cudaMemcpyAsync(d_rays[i], rayons.data() + i * raysPerDevice, sizeof(GPURay) * raysPerDevice, cudaMemcpyHostToDevice, stream[i]));

            CHECK_CUDA_ERRORS(cudaMalloc(&d_results[i], sizeof(int) * raysPerDevice));

            d_tree[i] = M_root_gpu_tree->DeepCopy(M_root_gpu_tree);

            // Launch the kernel with one block per ray
            GPU_traverse_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_tree[i], d_rays[i], d_results[i], raysPerDevice);

            // Copy back results, the first fustrum wil be place from results[0] to results[raysPerDevice - 1] and so on
            CHECK_CUDA_ERRORS(cudaMemcpyAsync(results.data() + i * raysPerDevice, d_results[i], sizeof(int) * raysPerDevice, cudaMemcpyDeviceToHost, stream[i]));
        }

        for (int i = 0; i < numDevices; i++) 
        {
            CHECK_CUDA_ERRORS(cudaSetDevice(i));
            CHECK_CUDA_ERRORS(cudaStreamSynchronize(stream[i]));
            CHECK_CUDA_ERRORS(cudaStreamDestroy(stream[i]));
            CHECK_CUDA_ERRORS(cudaFree(d_tree[i]));
            CHECK_CUDA_ERRORS(cudaFree(d_rays[i]));
            CHECK_CUDA_ERRORS(cudaFree(d_results[i]));
        }

        return results;
    }
};

// Wrapper to be able to call the GPU function from the CPU
std::vector<int> GPUraySearchWrapper(std::vector<Feel::BVHRay> const& rays, const Feel::BVHTree<3> * tree)
{
    return GPUBVH::GPUraySearch(rays, tree);
}
