
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
        //matrix_file.open(matrix_filename,std::ios_base::out);
        matrix_file.open(matrix_filename,std::ios_base::trunc);
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


    template <typename MeshType>
    bool 
    ShadingMask<MeshType>::testShadingMaskComparisonLevel1(
        std::string shadingMaskFolder,std::string building_name,std::string marker_name)
    {
        bool QCTRL=false; 
        std::string matrix_filename_NEW  = shadingMaskFolder+"/SM_Matrix_"+building_name+"_"+marker_name+".csv";
        std::string matrix_filename_CTRL = shadingMaskFolder+"/SM_Matrix_CTRL_"+building_name+"_"+marker_name+".csv";

        //cout<<"Test mat:"<<building_name+"_"+marker_name+".csv \n"; 

        std::ifstream FICH_NEW,FICH_CTRL;
        FICH_NEW.open(matrix_filename_NEW);
        FICH_CTRL.open(matrix_filename_CTRL);
        std::string MatLine_NEW,MatLine_CTRL;
        bool QRep=true; int i=0;
        while ((!FICH_CTRL.eof())) 
        {
            std::getline(FICH_NEW,MatLine_NEW);
            std::getline(FICH_CTRL,MatLine_CTRL);
            QRep=QRep && (MatLine_NEW.compare(MatLine_CTRL)==0); i++;
            if (!QRep) { QCTRL=true; }
            /*if (QCTRL) {
                cout<<"*** ERROR !!! ***\n";
                cout<<"  Name NEW  = ["<<matrix_filename_NEW<<"]\n";
                cout<<"  Name CTRL = ["<<matrix_filename_CTRL<<"]\n";    
                cout<<"  In matrix line "<<i<<"\n";  
                cout<<"  Lgn Mat CRTL ="<<MatLine_CTRL<<"\n";
                cout<<"  Lgn Mat NEW  ="<<MatLine_NEW<<"\n";
                cout<<(MatLine_NEW.compare(MatLine_CTRL)==0);
            }
            */
        }
        FICH_NEW.close(); FICH_CTRL.close();
        if (QCTRL) { 
                    //cout<<"Test mat:"<<building_name+"_"+marker_name+".csv ===>"<<QRep<<"\n"; 
                    //getchar();
                     }
        return QRep;
    }


    template <typename MeshType>
    bool 
    ShadingMask<MeshType>::testShadingMaskComparisonLevel2(
        std::string shadingMaskFolder,std::string building_name,std::string marker_name)
    {
        bool QCTRL=false;
        std::string matrix_filename_NEW  = shadingMaskFolder+"/SM_Matrix_"+building_name+"_"+marker_name+".csv";
        std::string matrix_filename_CTRL = shadingMaskFolder+"/SM_Matrix_CTRL_"+building_name+"_"+marker_name+".csv";

        //cout<<"Test mat:"<<building_name+"_"+marker_name+".csv \n"; 

        std::ifstream FICH_NEW,FICH_CTRL;
        FICH_NEW.open(matrix_filename_NEW);
        FICH_CTRL.open(matrix_filename_CTRL);
        std::string MatLine_NEW,MatLine_CTRL;
        bool QRep=true; int l=0;
        std::vector<std::vector<std::string> > parsedCsv_CTRL;
        std::vector<std::vector<std::string> > parsedCsv_NEW;

        int LiCTRL,LjCTRL;
        int LiNEW,LjNEW;

        double averageCTRL=0.0;
        double averageNEW=0.0;

        while ((!FICH_CTRL.eof())) 
        {   l++;
            std::getline(FICH_NEW,MatLine_NEW);
            std::getline(FICH_CTRL,MatLine_CTRL);
            std::stringstream lineStream_CTRL(MatLine_CTRL);
            std::stringstream lineStream_NEW(MatLine_NEW);
            std::string cell_CTRL;
            std::string cell_NEW;
            std::vector<std::string> parsedRow_CTRL;
            std::vector<std::string> parsedRow_NEW;
            while(std::getline(lineStream_CTRL,cell_CTRL,',')) { parsedRow_CTRL.push_back(cell_CTRL);}
            parsedCsv_CTRL.push_back(parsedRow_CTRL);

            while(std::getline(lineStream_NEW,cell_NEW,',')) { parsedRow_NEW.push_back(cell_NEW); }
            parsedCsv_NEW.push_back(parsedRow_NEW);
        }

        LjCTRL=parsedCsv_CTRL.size()-1; LiCTRL=parsedCsv_CTRL[0].size();
        LjNEW=parsedCsv_NEW.size()-1;   LiNEW=parsedCsv_NEW[0].size();
        for (int i=0; i<LiCTRL;i++) 
        {
                for (int j=0; j<LjCTRL;j++) {  averageCTRL=averageCTRL+strtod(parsedCsv_CTRL[j][i].c_str(), NULL);
                    //cout<<parsedCsv_CTRL[j][i]<<" "; 
                }
                //cout<<"\n";
                for (int j=0; j<LjCTRL;j++) {  averageNEW=averageNEW+strtod(parsedCsv_NEW[j][i].c_str(), NULL);
                    //cout<<parsedCsv_NEW[j][i]<<" "; 
                }
                //cout<<"\n";
        }

        averageCTRL=averageCTRL/(double(LiCTRL)*double(LjCTRL));
        averageNEW=averageNEW/(double(LiNEW)*double(LjNEW));

         
        if (QCTRL) {
            //cout<<"*** ERROR !!! ***\n";
            cout<<"  Name NEW  = ["<<matrix_filename_NEW<<"]\n";
            cout<<"  Name CTRL = ["<<matrix_filename_CTRL<<"]\n";    
            cout<<"  In matrix line "<<l<<"\n";  
            cout<<"  Lgn Size CTRL ="<<LiCTRL<<":"<<LjCTRL<<"\n";
            cout<<"  Lgn Size NEW  ="<<LiNEW<<":"<<LjNEW<<"\n"; 
            cout<<"  averageCTRL   ="<<std::fixed<<std::setprecision(9)<<averageCTRL<<"\n"; 
            cout<<"  averageNEW    ="<<std::fixed<<std::setprecision(9)<<averageNEW<<"\n"; 
        }
        FICH_NEW.close(); FICH_CTRL.close();
        //QCTRL=true;
        if (QCTRL) { cout<<"Test mat:"<<building_name+"_"+marker_name+".csv ===>"<<QRep<<"\n"; 
                    //getchar();
        }
        return QRep;
    }

    template <typename MeshType>
    void 
    ShadingMask<MeshType>::saveMetadata()
    {
        std::string folder = (boost::filesystem::path(Environment::appRepository())).string();
        if (!boost::filesystem::exists(folder))
            boost::filesystem::create_directory(folder);
        
        //std::string json_filename = folder+"/shadingmask_metadata.json";
        std::string json_filename = folder+"/shadingmask_metadata_"
            +std::to_string(M_Nthreads)+"_"
            +std::to_string(M_Nrays)+".json";

        std::ofstream o(json_filename);
        o << std::setw(6) << M_metadataJson << std::endl;
                
    }
}