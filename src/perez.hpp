#include <iostream>
#include <cpr/cpr.h>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <curl/curl.h>
#include <filesystem>
#include <fmt/core.h>
#include <feel/feelmesh/sphere.hpp>
#include <feel/feelcore/environment.hpp>

/* library developped by NREL in order to compute the sun's position in the sky at a specific time and location */
extern "C" {
#include "spa.h"
}

/* define the normal incident extraterrestrial irradiance (W/m^2)*/
const double SOLAR_CONSTANT = 1366.0; 

/* define the Perez Coefficients table */
const std::vector< std::vector<double> > PerezCoeff = {
    // epsilon in [ 1.000 ; 1.065 ]
    {1.3525, -0.2576, -0.2690, -1.4366, -0.7670,
     0.0007, 1.2734, -0.1233, 2.8000, 0.6004,
     1.2375, 1.0000, 1.8734, 0.6297, 0.9738,
     0.2809, 0.0356, -0.1246, -0.5718, 0.9938},
    // epsilon in [1.065 ; 1.230]
    {-1.2219, -0.7730, 1.4148, 1.1016, -0.2054,
     0.0367, -3.9128, 0.9156, 6.9750, 0.1774,
     6.4477, -0.1239, -1.5798, -0.5081, -1.7812,
     0.1080, 0.2624, 0.0672, -0.2190, -0.4285},
    // epsilon in [1.230 ; 1.500]
    {-1.1000, -0.2515, 0.8952, 0.0156, 0.2782,
     -0.1812, -4.5000, 1.1766, 24.7219, -13.0812,
     -37.7000, 34.8438, -5.0000, 1.5218, 3.9229,
     -2.6204, -0.0156, 0.1597, 0.4199, -0.5562},
    // epsilon in [1.500 ; 1.950]
    {-0.5484, -0.6654, -0.2672, 0.7117, 0.7234,
     -0.6219, -5.6812, 2.6297, 33.3389, -18.3000,
     -62.2500, 52.0781, -3.5000, 0.0016, 1.1477,
     0.1062, 0.4659, -0.3296, -0.0876, -0.0329},
    // epsilon in [1.950 ; 2.800]
    {-0.6000, -0.3566, -2.5000, 2.3250, 0.2937,
     0.0496, -5.6812, 1.8415, 21.0000, -4.7656,
     -21.5906, 7.2492, -3.5000, -0.1554, 1.4062,
     0.3988, 0.0032, 0.0766, -0.0656, -0.1294},
    // epsilon in [2.800 ; 4.500]
    {-1.0156, -0.3670, 1.0078, 1.4051, 0.2875,
     -0.5328, -3.8500, 3.3750, 14.0000, -0.9999,
     -7.1406, 7.5469, -3.4000, -0.1078, -1.0750,
     1.5702, -0.0672, 0.4016, 0.3017, -0.4844},
    // epsilon in [4.500 ; 6.200]
    {-1.0000, 0.0211, 0.5025, -0.5119, -0.3000,
     0.1922, 0.7023, -1.6317, 19.0000, -5.0000,
     1.2438, -1.9094, -4.0000, 0.0250, 0.3844,
     0.2656, 1.0468, -0.3788, -2.4517, 1.4656},
    // epsilon in [6.200 ; ... ]
    {-1.0500, 0.0289, 0.4260, 0.3590, -0.3250,
     0.1156, 0.7781, 0.0025, 31.0625, -14.5000,
     -46.1148, 55.3750, -7.2312, 0.4050, 13.3500,
     0.6234, 1.5000, -0.6426, 1.8564, 0.5636}};

namespace Feel { 
namespace FS = std::filesystem;
/* Class enabling computation of the Perez all-weather Sky model, presented in 1993 */
class PerezSkyModel
{
public:
    /* Class members */
    std::string M_BuildingName;
    Eigen::VectorXd M_PerezParameters; // vector containing the 5 Perez parameters for each sun position (~for each hour)
    Eigen::MatrixXd M_SkyModel; // matrix containing the values of the sky model for a unique sun position
    std::vector<int> M_Time; // vector containing the hour as an int for each sun position where the sun is above the horizon

    /* hourly values to compute one model per hour, will all have the same size = round(sunset-sunrise) - 1 */
    std::vector<double> M_DirectRadiation, M_DiffuseRadiation, M_DirectNormalIrradiance, M_DiffuseIlluminance, M_Brightness, M_Clearness, M_AirMass, M_ZenithAngles, M_SolarAzimuth;
    int M_azimuthSize, M_altitudeSize;

    /* Constructor */
    PerezSkyModel()
    {
        M_azimuthSize = 72;
        M_altitudeSize = 10;
        // Initialize all matrices and vectors to default without giving any size
        M_PerezParameters = Eigen::VectorXd::Zero(5);
        M_SkyModel = Eigen::MatrixXd::Zero(0, 0);
        M_Time = std::vector<int>();
        M_DirectRadiation = std::vector<double>();
        M_DiffuseRadiation = std::vector<double>();
        M_DirectNormalIrradiance = std::vector<double>();
        M_DiffuseIlluminance = std::vector<double>();
    }

    /* Retrieve data using curl */
    void RetrieveData(double longitude, double latitude, std::string start_date)
    {
        /* setting the path for the resulting file */
        std::ofstream file;
        FS::path resultPath = FS::current_path().parent_path().parent_path().parent_path();
        FS::path filePath = resultPath / "solar-shading-perez" / "solar-shading" / "results" / "meteo" / "openmeteo.json";
        ::CURL *curl;
        FILE *fp;
        CURLcode res;
        std::string url = "https://archive-api.open-meteo.com/v1/archive?latitude=" + std::to_string(latitude) + "&longitude=" + std::to_string(longitude) + "&start_date=" + start_date.c_str() + "&end_date=" + start_date.c_str() + "&hourly=is_day,direct_radiation,diffuse_radiation,direct_normal_irradiance";
        curl = curl_easy_init();
        if (curl) {
            fp = fopen(filePath.string().c_str(),"wb");
            curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, NULL);
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
            res = curl_easy_perform(curl);
            /* always cleanup */
            curl_easy_cleanup(curl);
            fclose(fp);
        }
        else {
            throw std::logic_error( "Unable to open file" );
        }
    }

    /* Retrieve data using cpr */
    // void RetrieveData(double longitude, double latitude, std::string start_date)
    // {
    //     /* setting the path for the resulting file */
    //     std::ofstream file;
    //     FS::path resultPath = FS::current_path().parent_path().parent_path();
    //     FS::path filePath = resultPath / "solar-shading-perez" / "solar-shading" / "results" / "meteo" / "openmeteo.json";
    //     std::cout << "The path is : " << filePath << std::endl;

    //     std::string url = "https://archive-api.open-meteo.com/v1/archive?latitude=" + std::to_string(latitude) + "&longitude=" + std::to_string(longitude) + "&start_date=" + start_date.c_str() + "&end_date=" + start_date.c_str() + "&hourly=is_day,direct_radiation,diffuse_radiation,direct_normal_irradiance";

    //     auto response = cpr::Get(cpr::Url{url});

    //     if (response.status_code == 200) { // HTTP status 200 means OK
    //         std::ofstream outfile(filePath);
    //         if(outfile.is_open()){
    //             outfile << response.text;
    //             outfile.close();
    //         } else {
    //             std::cout << "Unable to open file at " << filePath << std::endl;
    //         }
    //     } else {
    //         std::cout << "HTTP request failed with status " << response.status_code << std::endl;
    //     }
    // }

    /* Verify correctness of the data, if the direct / diffuse radiation or direct normal irradiance is equal to zero, replace with 1 */
    void checkData()
    {
        for (size_t i = 0; i < M_DirectRadiation.size(); i++)
        {
            if (M_DirectRadiation[i] == 0)
            {
                M_DirectRadiation[i] = 1.;
            }
            if (M_DiffuseRadiation[i] == 0)
            {
                M_DiffuseRadiation[i] = 1.;
            }
            if (M_DirectNormalIrradiance[i] == 0)
            {
                M_DirectNormalIrradiance[i] = 1.;
            }
        }
    }

    /* Correct the solar Zenith Angle, if it is greater than π/2 , set it to π/2 */
    void checkSolarAngles()
    {
        for (size_t i = 0; i < M_ZenithAngles.size(); i++)
        {
            if (M_ZenithAngles[i] > M_PI / 2.0)
            {
                M_ZenithAngles[i] = M_PI / 2.0 + 1e-6;
            }
        }
    }
    
    /* convert radians to degrees */
    double deg(double rad)
    {
        return rad * 180.0 / M_PI;
    }

    /* convert degrees to radians */
    double rad(double deg)
    {
        return deg * M_PI / 180.0;
    }

    /* Compute the sun's Position thanks to the NREL method */
    std::vector<std::pair<double, double>> computeSunPosition(int year, int month, int day, float latitude, float longitude, int timezone, const std::vector<int>& Time) 
    {
        std::vector<std::pair<double, double>> sun_positions;
        for (int hour : Time) {
            spa_data spa;  // declare the SPA structure
            
            // Set the SPA structure values
            spa.year = year;
            spa.month = month;
            spa.day = day;
            spa.hour = hour;
            spa.minute = 0;
            spa.second = 0;
            spa.timezone = timezone;
            spa.longitude = longitude;
            spa.latitude = latitude;
            spa.elevation = 0;
            spa.pressure = 820; // average value in millibar
            spa.temperature = 11; // average value in Celsius
            spa.slope = 0;
            spa.azm_rotation = 0;
            spa.atmos_refract = 0.5667; // average value
            spa.function = SPA_ALL;

            // Call the SPA calculate function and print the results
            int result = spa_calculate(&spa);
            if (result == 0) {  // check for SPA errors
                sun_positions.push_back(std::make_pair(spa.azimuth, spa.zenith));
            } else {
                throw std::logic_error("SPA Error Code: " + std::to_string(result));
            }
        }
        return sun_positions;
    }

    /* Rewrite the meteorological data to std::vectors, the startAppending term is used in order */
    /* to start appending the data only when the sun is in sight */
    std::vector<int> getTime(nl::json j) 
    {
        std::vector<std::string> allTime = j["hourly"]["time"].get<std::vector<std::string>>();
        std::vector<int> is_day = j["hourly"]["is_day"].get<std::vector<int>>();
        std::vector<int> Time;

        bool startAppending = false;
        
        for (size_t i = 0; i < is_day.size(); i++) {
            if (is_day[i] == 1) {
                if (startAppending) {
                    std::string strTime = allTime[i].substr(allTime[i].find('T') + 1, 2);
                    Time.push_back(std::stoi(strTime));
                } else {
                    startAppending = true;
                }
            }
        }
        return Time;
    }

    /* input : nl::json file || output : std::vector<double> containing the Direct Normal Irradiance */
    std::vector<double> getDirectNormalIrradiance(nl::json j) 
    {
        std::vector<int> is_day = j["hourly"]["is_day"].get<std::vector<int>>();
        std::vector<double> allDirectNormalIrradiance = j["hourly"]["direct_normal_irradiance"].get<std::vector<double>>();
        std::vector<double> DirectNormalIrrandiance;

        bool startAppending = false;
        for (size_t i = 0; i < is_day.size(); i++) {
            if (is_day[i] == 1) {
                if (startAppending) {
                    DirectNormalIrrandiance.push_back(allDirectNormalIrradiance[i]);
                } else {
                    startAppending = true;
                }
            }
        }
        return DirectNormalIrrandiance;
    }
    /* input : nl::json file || output : std::vector<double> containing the Direct Radiation */
    std::vector<double> getDirectRadiation(nl::json j) 
    {
        std::vector<int> is_day = j["hourly"]["is_day"].get<std::vector<int>>();
        std::vector<double> allDirectRadiation = j["hourly"]["direct_radiation"].get<std::vector<double>>();
        std::vector<double> DirectRadiation;

        bool startAppending = false;
        for (size_t i = 0; i < is_day.size(); i++) {
            if (is_day[i] == 1) {
                if (startAppending) {
                    DirectRadiation.push_back(allDirectRadiation[i]);
                } else {
                    startAppending = true;
                }
            }
        }
        return DirectRadiation;
    }

    /* input : nl::json file || output : std::vector<double> containing the Diffuse Radiation */
    std::vector<double> getDiffuseRadiation(nl::json j) 
    {
        std::vector<int> is_day = j["hourly"]["is_day"].get<std::vector<int>>();
        std::vector<double> allDiffuseRadiation = j["hourly"]["diffuse_radiation"].get<std::vector<double>>();
        std::vector<double> DiffuseRadiation;

        bool startAppending = false;
        for (size_t i = 0; i < is_day.size(); i++) {
            if (is_day[i] == 1) {
                if (startAppending) {
                    DiffuseRadiation.push_back(allDiffuseRadiation[i]);
                } else {
                    startAppending = true;
                }
            }
        }
        return DiffuseRadiation;
    }

    /* input : clearness (ε) || output : int giving the index in order to set the Perez parameters */
    int getRange(double epsilon)
    {
        if (epsilon < 1.0 )
        {
            throw std::logic_error("epsilon must be greater than 1.0");
        }
        else if (epsilon >= 1.000 && epsilon < 1.065 )
        {
            return 0;
        }
        else if (epsilon >= 1.065 && epsilon < 1.230 )
        {
            return 1;
        }
        else if (epsilon >= 1.230 && epsilon < 1.500 )
        {
            return 2;
        }
        else if (epsilon >= 1.500 && epsilon < 1.950 )
        {
            return 3;
        }
        else if (epsilon >= 1.950 && epsilon < 2.800 )
        {
            return 4;
        }
        else if (epsilon >= 2.800 && epsilon < 4.500 )
        {
            return 5;
        }
        else if (epsilon >= 4.500 && epsilon < 6.200 )
        {
            return 6;
        }
        else if (epsilon >= 6.200 )
        {
            return 7;
        }
        else
        {
            throw std::logic_error("Problem occured during epsilon range determination, epsilon is equal to : " + std::to_string(epsilon) + "\n");
        }
    }

    /* Computation of the five Perez parameters a, b, c, d and e given a specific */
    /* brightness (Δ), clearness (ε) and ZenithAngle (Z) are stored as class members */
    void setPerezParameters(double epsilon, double delta, double ZenithAngle) 
    {
        int range = getRange(epsilon);
        Eigen::MatrixXd temp(5, 4); /* i = 0 : a || i = 1 : b || i = 2 : c || i = 3 : d || i = 4 : e */
        for (int i = 0 ; i < 5 ; i++)
        {
            for (int j = 0 ; j < 4 ; j++)
            {
                temp(i, j) = PerezCoeff[range][i * 4 + j];
            }
        }
        // Compute Perez Parameters
        for (int i = 0 ; i < 5 ; i++)
        {
            M_PerezParameters(i) = temp(i, 0) + temp(i, 1) * ZenithAngle + delta * (temp(i, 2) + temp(i, 3) * ZenithAngle);
        }
        // Compute special cases for epsilon in [1.000 ; 1.065] for c and d
        if (range == 0)
        {
            M_PerezParameters(2) = std::exp(std::pow(delta * (temp(2, 0) + temp(2, 1) * ZenithAngle), temp(2, 2))) - temp(2, 3);
            M_PerezParameters(3) = -std::exp(delta*(temp(3, 0) + temp(3, 1) * ZenithAngle)) + temp(3, 2) + delta * temp(3, 3);
        }
    }

    /* Computation of the optical air mass thanks to a simple model. Takes into account the */
    /* curvature of the earth so that the optical air mass at sunrise wouldn't be infinite */
    /* Kasten and Young's formula for air mass */
    double getAirMass( double ZenithAngle )
    {
        if (ZenithAngle > M_PI / 2.0)
        {
            ZenithAngle = M_PI / 2.0;
        }
        else if (ZenithAngle < 0.0)
        {
            throw std::logic_error("Zenith angle must be greater than 0");
        }
        /* The zenith angle is assumed to be given in radians */
        return (1.0 / (std::cos(ZenithAngle) + 0.50572 * std::pow(96.07995 - deg(ZenithAngle), -1.6364)));
    }

    /* Computation of the optical air mass thanks to a simple model. Takes into account the */
    /* curvature of the earth so that the optical air mass at sunrise wouldn't be infinite */
    /* Kasten and Young's formula for air mass , presented in : */
    /* F. Kasten and A. T. Young. 1989. "Revised optical air mass tables and approximation formula." Applied Optics 28, 4735-4738 */
    std::vector<double> computeAirMass(std::vector<double> ZenithAngles)
    {
        std::vector<double> airMass;
        for (size_t i = 0; i < ZenithAngles.size(); i++)
        {
            airMass.push_back(getAirMass(ZenithAngles[i]));
        }
        return airMass;
    }

    /* Computation of the sky's clearness (ε) with respect to the horizontal diffuse */
    /* irradiance (Eed), the normal incident direct irradiance(Ees) and the solar zenith angle (Z) */
    double getClearness(double Eed, double Ees, double Z)
    {
        return (( Eed + Ees ) / Eed + 1.041 * std::pow( Z , 3.0)) / ( 1.0 + 1.041 * std::pow( Z , 3.0));
    }
    
    /* Computation of the sky's clearness (ε) with respect to the horizontal diffuse */
    /* irradiance (Eed), the normal incident direct irradiance(Ees) and the solar zenith angle (Z) */
    std::vector<double> computeClearness(std::vector<double> Eed, std::vector<double> Ees, std::vector<double> Z)
    {
        std::vector<double> clearness;
        for (size_t i = 0; i < Eed.size(); i++)
        {
            clearness.push_back(getClearness(Eed[i], Ees[i], Z[i]));
        }
        return clearness;
    }

    /* Computation of the sky's brightness (Δ) with respect to the horizontal diffuse */
    /* irradiance (Eed), the optical air mass (m) and the normal incident extraterrestrial */
    /* irradiance (Eeso, here SOLAR_CONSTANT) */
    double getBrightness(double Eed, double m)
    {
        return m * Eed / SOLAR_CONSTANT;
    }

    /* Computation of the sky's brightness (Δ) with respect to the horizontal diffuse */
    /* irradiance (Eed), the optical air mass (m) and the normal incident extraterrestrial */
    /* irradiance (Eeso, here SOLAR_CONSTANT) */
    std::vector<double> computeBrightness(std::vector<double> Eed, std::vector<double> m)
    {
        std::vector<double> brightness;
        for (size_t i = 0; i < Eed.size(); i++)
        {
            brightness.push_back(getBrightness(Eed[i], m[i]));
        }
        return brightness;
    }

    /* Compute the relative sky luminance (lv) */
    double lv(double zeta, double gamma)
    {
        // check if zeta is near pi/2
        if (std::abs(zeta) > M_PI/2 - 1e-2)
        {
            zeta = M_PI/2;
            // std::cout << "Warning, zeta is near pi/2 , zeta = " << zeta << " and cos(zeta) = " << std::cos(zeta) << " and b = " << M_PerezParameters(1) << std::endl;
        }
        double lv = (1.0 + M_PerezParameters(0) * std::exp(M_PerezParameters(1) / (std::cos(zeta) + 1e-6))) * (1.0 + M_PerezParameters(2) * std::exp(M_PerezParameters(3) * gamma) + M_PerezParameters(4) * std::pow(std::cos(gamma), 2.0));
        return lv;
    }

    /* Compute the integral of the relative sky luminance over the whole sky hemisphere */
    /* in order to normalize the resulting sky model */
    double compute_integral(double solarAzimuth, double solarZenith, int num_points = 200)
    {
        double integral_lv = 0.0;
        double delta_theta = M_PI / (2.0 * num_points); 
        double delta_phi = 2.0 * M_PI / num_points; 

        for(int i = 0; i < num_points; ++i)
        {
            double theta = (i + 0.5) * delta_theta;
            double zeta = M_PI / 2.0 - theta;

            for(int j = 0; j < num_points; ++j)
            {
                double phi = (j + 0.5) * delta_phi;
                double cosGamma = std::sin(solarZenith) * std::sin(zeta) * std::cos(std::abs(solarAzimuth - phi)) + std::cos(solarZenith) * std::cos(zeta);
                double gamma = std::acos(cosGamma);
                double lv_value = lv(zeta, gamma);
                
                // Check for NaN
                // if (std::isnan(lv_value))
                // {
                //     throw std::logic_error("lv() returned NaN for these values: Zeta = " + std::to_string(zeta) + "; gamma = " + std::to_string(gamma) + "; solarAzimuth = " + std::to_string(solarAzimuth) + "; solarZenith = " + std::to_string(solarZenith));
                // }
                // // Check for inf
                // if (std::isinf(lv_value))
                // {
                //     throw std::logic_error("lv() returned inf for these values: Zeta = " + std::to_string(zeta) + "; gamma = " + std::to_string(gamma) + "; solarAzimuth = " + std::to_string(solarAzimuth) + "; solarZenith = " + std::to_string(solarZenith));
                // }

                integral_lv += lv_value * std::cos(zeta);
            }
        }
        integral_lv *= delta_theta * delta_phi;

        // Check for NaN
        // if (std::isnan(integral_lv))
        // {
        //     throw std::logic_error("Integral is NaN for these values: solarAzimuth = " + std::to_string(solarAzimuth) + "; solarZenith = " + std::to_string(solarZenith));
        // }
        // // Check for inf
        // if (std::isinf(integral_lv))
        // {
        //     throw std::logic_error("Integral is inf for these values: solarAzimuth = " + std::to_string(solarAzimuth) + "; solarZenith = " + std::to_string(solarZenith));
        // }

        return integral_lv;
    }

    double getNormalizedLuminance(double zeta, double gamma, double Evd, double solarAzimuth, double solarZenith, double integral_value)
    {
        double lv_value = lv(zeta, gamma);
        double result = lv_value * Evd / integral_value;

        return result;
    }


    /* method to compute the diffuse illuminance (Evd), which can be estimated from the diffuse */
    /* solar radiation (Eed) and depends on sky conditions, but this can't be handled using the open-meteo */
    /* api since the available data isn't sufficient for the computation of the conversion factor */
    /* namely the measurements of the solar spectral irradiance (i.e. the irradiance at different wavelengths) */
    /* instead, we will use a commonly used average for the conversion of sunlight irradiance to illuminance*/
    std::vector<double> getDiffuseIlluminance(std::vector<double> Eed)
    {
        std::vector<double> Evd;
        for (size_t i = 0; i < Eed.size(); i++)
        {
            Evd.push_back(Eed[i] * 120.);
        }
        return Evd;
    }

    /* Saves the current Sky Model M_SkyModel (the hour is written in the name of the file)*/
    void saveSkyModel(int time)
    {
        int i_long=0;
        int j_alt=0;
        
        Eigen::MatrixXd matrix_sm(M_azimuthSize,M_altitudeSize);
        
        // Save the matrix into a CSV file, inside the shadingMasks subfolder of the results folder
        std::ofstream matrix_file;        
        std::string SkyModelFolder = (boost::filesystem::path(Environment::appRepository())/("skyModels")).string();
        if (!boost::filesystem::exists(SkyModelFolder))
            boost::filesystem::create_directory(SkyModelFolder);
        
        std::string matrix_filename = SkyModelFolder+"/SkyModel_Matrix_"+M_BuildingName+"_"+std::to_string(time)+"H.csv";
        matrix_file.open(matrix_filename,std::ios_base::out);
        for(int i=0; i<M_azimuthSize; i++)
        {
            for(int j=0; j<M_altitudeSize-1; j++)
            {                
                matrix_file << M_SkyModel(i,j) << ",";
            }              
            matrix_file << M_SkyModel(i,M_altitudeSize-1) << "\n";     
        }
        matrix_file.close();
    }

    /* Initialize and compute the model */
    void computePerezSkyModel(double longitude, double latitude, std::string start_date, nl::json const& specs, int AzimuthSize = 72, int AltitudeSize = 10)
    {
        RetrieveData(longitude, latitude, start_date);
        M_BuildingName = specs["Buildings"][0];
        std::cout << fmt::format("Start Computation of the Sky Model for the building {}\n", M_BuildingName);

        FS::path resultPath = FS::current_path().parent_path().parent_path().parent_path();
        FS::path filePath = resultPath / "solar-shading-perez" / "solar-shading" / "results" / "meteo" / "openmeteo.json";
        std::ifstream jsonFile(filePath);

        // Check if the file was opened successfully
        if (!jsonFile) {
            throw std::logic_error( "Unable to open JSON file" );
        }

        // Parse the JSON file into a json object
        nlohmann::json j;
        jsonFile >> j;

        // Close the file stream
        jsonFile.close();

        int year = std::stoi(start_date.substr(0, 4)); // Year: 2023
        int month = std::stoi(start_date.substr(5, 2)); // Month: 07
        int day = std::stoi(start_date.substr(8, 2)); // Day: 21  

        // get the time zone, supposing the longitude and latitude are in France, taking into account the summer time
        int timezone = 1;
        if (month >= 3 && month <= 10) {
            timezone = 2;
        }

        // Compute the sun's position
        std::vector<std::pair<double, double>> sun_positions = computeSunPosition(year, month, day, latitude, longitude, timezone, getTime(j));

        // set the Angles, and convert them to radians
        for (size_t i = 0; i < sun_positions.size(); i++) {
            M_SolarAzimuth.push_back(rad(sun_positions[i].first));
            M_ZenithAngles.push_back(rad(sun_positions[i].second));
        }
        checkSolarAngles();
        
        // Instanctiate 
        M_SkyModel = Eigen::MatrixXd::Zero(AzimuthSize, AltitudeSize);

        // Access the time, direct radiation, diffuse radiation, direct normal irradiance
        M_azimuthSize = AzimuthSize;
        M_altitudeSize = AltitudeSize;
        M_Time = getTime(j);
        M_DirectRadiation = getDirectRadiation(j);
        M_DiffuseRadiation = getDiffuseRadiation(j);
        M_DirectNormalIrradiance = getDirectNormalIrradiance(j);

        // Verify the correctness of the downloaded data
        checkData();

        // Compute the diffuse illuminance, air mass, clearness and brightness
        M_DiffuseIlluminance = getDiffuseIlluminance(M_DiffuseRadiation);
        M_AirMass = computeAirMass(M_ZenithAngles);
        M_Clearness = computeClearness(M_DiffuseRadiation, M_DirectNormalIrradiance, M_ZenithAngles);
        M_Brightness = computeBrightness(M_DiffuseRadiation, M_AirMass);

        // Check the data's values by doing a std::cout
        std::cout << "The time is : ";
        for (size_t i = 0; i < M_Time.size(); i++) {
            std::cout << M_Time[i] << " ";
        }
        std::cout << "\nThe direct radiation is : ";
        for (size_t i = 0; i < M_DirectRadiation.size(); i++) {
            std::cout << M_DirectRadiation[i] << " ";
        }
        std::cout << "\nThe diffuse radiation is : ";
        for (size_t i = 0; i < M_DiffuseRadiation.size(); i++) {
            std::cout << M_DiffuseRadiation[i] << " ";
        }
        std::cout << "\nThe direct normal irradiance is : ";
        for (size_t i = 0; i < M_DirectNormalIrradiance.size(); i++) {
            std::cout << M_DirectNormalIrradiance[i] << " ";
        }
        std::cout << "\nThe diffuse illuminance is : ";
        for (size_t i = 0; i < M_DiffuseIlluminance.size(); i++) {
            std::cout << M_DiffuseIlluminance[i] << " ";
        }
        std::cout << "\nThe air mass is : ";
        for (size_t i = 0; i < M_AirMass.size(); i++) {
            std::cout << M_AirMass[i] << " ";
        }
        std::cout << "\nThe clearness is : ";
        for (size_t i = 0; i < M_Clearness.size(); i++) {
            std::cout << M_Clearness[i] << " ";
        }
        std::cout << "\nThe brightness is : ";
        for (size_t i = 0; i < M_Brightness.size(); i++) {
            std::cout << M_Brightness[i] << " ";
        }
        std::cout << "\nThe solar azimuth is : ";
        for (size_t i = 0; i < M_SolarAzimuth.size(); i++) {
            std::cout << M_SolarAzimuth[i] << " ";
        }
        std::cout << "\nThe zenith angle is : ";
        for (size_t i = 0; i < M_ZenithAngles.size(); i++) {
            std::cout << M_ZenithAngles[i] << " ";
        }
        std::cout << "\n";
        
        std::cout << "==============================================================================================================\n";
        // Loop over the hours of the day in order to compute all models 
        for (int i = 0; i < M_Time.size(); i++)
        {
            int hour = M_Time[i];
            double delta = M_Brightness[i];  
            double epsilon = M_Clearness[i];  
            double A = M_SolarAzimuth[i];  
            double Z = M_ZenithAngles[i];
            setPerezParameters(epsilon, delta, Z);
            double integral_value = compute_integral(A, Z);
            std::cout << "The Perez Parameters at hour " << hour << " are : " << M_PerezParameters(0) << "   " << M_PerezParameters(1) << "   " << M_PerezParameters(2) << "   " << M_PerezParameters(3) << "   " << M_PerezParameters(4) << "\n";

            for (int j = 0; j < AltitudeSize; j++) 
            {
                for (int k = 0; k < AzimuthSize; k++) 
                {
                    double z = M_PI / 2.0 - (double)j / (double)AltitudeSize * M_PI / 2.0;  // zenith of the sky element
                    double a = (double)k / (double)AzimuthSize * 2.0 * M_PI;  // azimuth of the sky element
                    if ( j == 0 )
                        z = M_PI / 2.0 - 1e-6; // avoid division by zero
                    
                    // Calculate gamma with error checking
                    double gamma;
                    double cos_gamma = std::cos(Z) * std::cos(z) + std::sin(Z) * std::sin(z) * std::cos(std::abs(A - a));
                    if (cos_gamma > 1 && cos_gamma < 1.1)
                        gamma = 0;
                    else if (cos_gamma > 1.1)
                    {
                        throw std::logic_error ("Error in calculation of gamma (angle between point and sun)");
                    }
                    else
                    {
                        gamma = std::acos(cos_gamma);
                    }
                    M_SkyModel(k, j) += getNormalizedLuminance(z, gamma, M_DiffuseIlluminance[i], A, Z, integral_value);
                }
            }
            saveSkyModel(hour);
        }
    }
};
}