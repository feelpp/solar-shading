
namespace Feel
{

    template <typename MeshType>
    void 
    ShadingMask<MeshType>::saveShadingMask(std::string building_name, std::string marker_name, const Eigen::Ref<const Eigen::MatrixXd>& M)
    {
        int i_long=0;
        int j_alt=0;

        Eigen::MatrixXd matrix_sm(M_azimuthSize,M_altitudeSize);

        // Save the matrix into a CSV file, inside the shadingMasks subfolder of the results folder
        std::ofstream matrix_file;
        std::string shadingMaskFolder = (boost::filesystem::path(Environment::appRepository())/("shadingMasks")).string();
        if (!boost::filesystem::exists(shadingMaskFolder))
            boost::filesystem::create_directory(shadingMaskFolder);

        std::string matrix_filename = shadingMaskFolder+"/SM_Matrix_"+building_name+"_"+marker_name+".csv";
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

}