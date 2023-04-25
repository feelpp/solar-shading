/**
 * @file laplacian.cpp
 * @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 * @brief Laplacian
 * @version 0.1
 * @date 2022-10-22
 *
 * @copyright Copyright (c) 2022 Université de Strasbourg
 *
 */
#include <iostream>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/minmax.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelts/bdf.hpp>

namespace Feel
{
inline const int FEELPP_DIM=2;
inline const int FEELPP_ORDER=1;

static inline const bool do_print = true;
static inline const bool dont_print = false;

/**
 * @brief compute the summary of a container
 * 
 * @tparam Container type of the container
 * @param c container
 * @param print boolean, true print the summary, false otherwise
 * @return nl::json json object containing the summary
 */
template<typename Container>
nl::json summary( Container const& c, bool print = do_print )
{
    using namespace Feel;
    using namespace Feel::vf;
    nl::json j;
    j["size"] = c.size();
    auto r = minmaxelt(_range = elements(support(c.functionSpace())), _element = c);
    j["min"] = r[0];
    j["max"] = r[1];
    j["mean"] = mean( _range = elements( c.mesh() ), _expr = idv( c ) );

    if (print)
    {
        if (Environment::isMasterRank())        
            std::cout << j.dump(2) << std::endl;
    }
    return j;
}
inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "rht options" );
    options.add_options()

        // mesh parameters
        ( "specs", Feel::po::value<std::string>(),
          "json spec file for rht" )

        ( "steady", Feel::po::value<bool>()->default_value( 1 ),
          "if 1: steady else unsteady" );

    return options.add( Feel::feel_options() );
}


template <int Dim = FEELPP_DIM, int Order = FEELPP_ORDER>
int runLaplacian( nl::json const& specs )
{
    using mesh_t = Mesh<Simplex<Dim>>;
    auto mesh = loadMesh( _mesh = new mesh_t, _filename = specs["/Meshes/laplacian/Import/filename"_json_pointer].get<std::string>() );
    Pch_ptrtype<mesh_t, Order> Xh;
    // define Xh on a marked region 
    if ( specs["/Spaces/laplacian/Domain"_json_pointer].contains("marker") )
        Xh = Pch<Order>(mesh, markedelements(mesh, specs["/Spaces/laplacian/Domain/marker"_json_pointer].get<std::vector<std::string>>()));
    // define Xh via a levelset phi where phi < 0 defines the Domain and phi = 0 the boundary
    else if (specs["/Spaces/laplacian/Domain"_json_pointer].contains("levelset"))
        Xh = Pch<Order>(mesh, elements(mesh, expr(specs["/Spaces/laplacian/Domain/levelset"_json_pointer].get<std::string>())));
    // define Xh on the whole mesh
    else
        Xh = Pch<Order>(mesh);

    auto u = Xh->element();
    auto v = Xh->element();

    auto a = form2( _test = Xh, _trial = Xh );
    auto at = form2( _test = Xh, _trial = Xh );
    auto l = form1( _test = Xh );
    auto lt = form1( _test = Xh );

    auto M_bdf = bdf( _space = Xh );

    M_bdf->start();

    // from now if the option "steady" is set to True then M_ bdf-setSteady will set time-step=time-final
    if ( boption("steady") )
        M_bdf->setSteady();

    for ( auto [key, material] : specs["/Models/laplacian/Materials"_json_pointer].items() )
    {
        LOG( INFO ) << fmt::format( "Material {} found", material );
        std::string mat = fmt::format( "/Materials/{}/k", material.get<std::string>() );
        auto k = specs[nl::json::json_pointer( mat )].get<std::string>();

        a += integrate( _range = markedelements( support( Xh ), material.get<std::string>() ), 
                _expr = M_bdf->polyDerivCoefficient( 0 ) * expr( k ) * gradt( u ) * trans( grad( v ) ) );
    }

    // BC Neumann
    if ( specs["/BoundaryConditions/laplacian"_json_pointer].contains( "flux" ) )
    {
        for ( auto& [bc, value] : specs["/BoundaryConditions/laplacian/flux"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "flux {}: {}", bc, value.dump() );
            auto flux = value["expr"].get<std::string>();

            l += integrate( _range = markedfaces( support( Xh ), bc ),
                    _expr = M_bdf->polyDerivCoefficient( 0 ) * expr( flux ) * id( v ) );
        }
    }

    // BC Robin
    if ( specs["/BoundaryConditions/laplacian"_json_pointer].contains( "convective_laplacian_flux" ) )
    {
        for ( auto& [bc, value] : specs["/BoundaryConditions/laplacian/convective_laplacian_flux"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "convective_laplacian_flux {}: {}", bc, value.dump() );
            auto h = value["h"].get<std::string>();
            auto Text = value["Text"].get<std::string>();

            a += integrate( _range = markedfaces( support( Xh ), bc ),
                    _expr = M_bdf->polyDerivCoefficient( 0 ) * expr( h ) * id( v ) * idt( u ) );
            l += integrate( _range = markedfaces( support( Xh ), bc ),
                    _expr = M_bdf->polyDerivCoefficient( 0 ) * expr( h ) * expr( Text ) * id( v ) );
        }
    }

    M_bdf->initialize( u );

    if ( boption("steady") )
        std::cout << "\n***** Steady state *****" << std::endl;
    else
        std::cout << "The step is  " << M_bdf->timeStep() << "\n"
                  << "The initial time is " << M_bdf->timeInitial() << "\n"
                  << "The final time is " << M_bdf->timeFinal() << "\n"
                  << "BDF order :  " << M_bdf->timeOrder() << "\n" << std::endl;

    // exporter mesh and fields for visualisation in paraview                  
    auto e = exporter(_mesh = mesh);

    // time loop
    for ( M_bdf->start(); M_bdf->isFinished()==false; M_bdf->next(u) )
    {
        at = a;
        lt = l;

        for ( auto [key, material] : specs["/Models/laplacian/Materials"_json_pointer].items() )
        {
            std::string matRho = fmt::format( "/Materials/{}/rho", material.get<std::string>() );
            std::string matCp = fmt::format( "/Materials/{}/Cp", material.get<std::string>() );
            auto Rho = specs[nl::json::json_pointer( matRho )].get<std::string>();
            auto Cp = specs[nl::json::json_pointer( matCp )].get<std::string>();

            lt += integrate( _range = markedelements( support( Xh ), material.get<std::string>() ), 
                    _expr = expr( Rho ) * expr( Cp ) * idv( M_bdf->polyDeriv() ) * id( v ) );
        }

        at.solve( _rhs = lt, _solution = u );
        // compute summary 
        summary(u);

        e->step(M_bdf->time())->addRegions();
        e->step(M_bdf->time())->add("u", u);
        e->save();
    }
    return 0;
}
} // namespace Feel

int main( int argc, char** argv )
{
    using namespace Feel;
    int status;
    try
    {
        Environment env( _argc = argc, _argv = argv,
                         _desc = makeOptions(),
                         _about = about( _name = fmt::format( "laplacian-{}dp{}", FEELPP_DIM, FEELPP_ORDER ),
                                         _author = "Feel++ Consortium",
                                         _email = "feelpp@cemosis.fr" ) );
        auto jsonfile = removeComments( readFromFile( Environment::expand( soption( "specs" ) ) ) );
        std::istringstream istr( jsonfile );
        json specs = json::parse( istr );
        return runLaplacian<FEELPP_DIM, FEELPP_ORDER>( specs );
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 0;
}