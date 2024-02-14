
#include <thread>
#include <vector>
#include <array>
#include <typeinfo>
#include <iostream>
#include <mutex>
#include <sched.h>
#include <pthread.h>

#include<algorithm> 
#include <string>
#include <utility>
#include <functional>
#include <future>
#include <cassert>
#include <chrono>
#include <type_traits>
#include <list>
#include <ranges>

#include "SpDataAccessMode.hpp"
#include "Utils/SpUtils.hpp"

#include "Task/SpTask.hpp"
#include "Legacy/SpRuntime.hpp"
#include "Utils/SpTimer.hpp"
#include "Utils/small_vector.hpp"
#include "Utils/SpConsumerThread.hpp"

#include "napp.hpp"

//=======================================================================================================================
//...
//=======================================================================================================================

constexpr auto& _parameters = NA::identifier<struct parameters_tag>;
constexpr auto& _task = NA::identifier<struct task_tag>;

namespace Backend{

    template<typename ...T, size_t... I>
    auto extractParametersAsTuple( std::tuple<T...> && t, std::index_sequence<I...>)
    {
        return std::forward_as_tuple( std::get<I>(t).getValue()...);
    }

    struct Runtime{
        template <typename ... Ts>
        void task(Ts && ... ts ) {
            auto t = std::make_tuple( std::forward<Ts>(ts)... );
            auto callback = std::get<sizeof...(Ts) - 1>(t);
            auto parameters = extractParametersAsTuple( std::move(t), std::make_index_sequence<sizeof...(Ts)-1>{} );
            std::apply( callback, std::move(parameters) );
        }
    };

    template <typename T,bool b>
    class SpData
    {
        static_assert(std::is_reference<T>::value,
                    "The given type must be a reference");
    public:
        using value_type = T;
        static constexpr bool isWrite = b;

        template <typename U, typename = std::enable_if_t<std::is_convertible_v<U,T>> >
        constexpr explicit SpData( U && u ) : M_val( std::forward<U>(u) ) {}

        constexpr value_type getValue() { return M_val; }
    private:
        value_type M_val;
    };

    template <typename T>
    auto spRead( T && t )
    {
        return SpData<T,false>{ std::forward<T>( t ) };
    }
    template <typename T>
    auto spWrite( T && t )
    {
        return SpData<T,true>{ std::forward<T>( t ) };
    }

    template<typename T>
    auto toSpData( T && t )
    {
        if constexpr ( std::is_const_v<std::remove_reference_t<T>> )
            return spRead( std::forward<T>( t ) );
        else
            return spWrite( std::forward<T>( t ) );
    }

    template<typename ...T, size_t... I>
    auto makeSpDataHelper( std::tuple<T...>& t, std::index_sequence<I...>)
    {
        return std::make_tuple( toSpData(std::get<I>(t))...);
    }
    template<typename ...T>
    auto makeSpData( std::tuple<T...>& t ){
        return makeSpDataHelper<T...>(t, std::make_index_sequence<sizeof...(T)>{});
    }

    template<typename T>
    auto toSpDataSpecx( T && t )
    {
        if constexpr ( std::is_const_v<std::remove_reference_t<T>> )
            return SpRead(std::forward<T>( t ));
        else
            return SpWrite(std::forward<T>( t ));
    }

    template<typename ...T, size_t... I>
    auto makeSpDataHelperSpecx( std::tuple<T...>& t, std::index_sequence<I...>)
    {
        return std::make_tuple( toSpDataSpecx(std::get<I>(t))...);
    }
    template<typename ...T>
    auto makeSpDataSpecx( std::tuple<T...>& t ){
        return makeSpDataHelperSpecx<T...>(t, std::make_index_sequence<sizeof...(T)>{});
    }
}


namespace Frontend
{
    template <typename ... Ts>
    void
    runTask( Ts && ... ts )
    {
        auto args = NA::make_arguments( std::forward<Ts>(ts)... );
        auto && task = args.get(_task);
        auto && parameters = args.get_else(_parameters,std::make_tuple());
        Backend::Runtime runtime;

        std::apply( [&runtime](auto... args){ runtime.task(args...); }, std::tuple_cat( Backend::makeSpData( parameters ), std::make_tuple( task ) ) );
    }

    template <typename ... Ts>
    auto parameters(Ts && ... ts)
    {
        //Construit un tuple de références aux arguments dans args pouvant être transmis en tant qu'argument à une fonction
        return std::forward_as_tuple( std::forward<Ts>(ts)... );
    }
}


//=======================================================================================================================
//...
//=======================================================================================================================


void *WorkerInNumCPU(void *arg) {
    std::function<void()> *func = (std::function<void()>*)arg;
    (*func)();
    pthread_exit(NULL);
}

class TasksDispach
{
    private:
        int nbThTotal;   
        std::vector<int> indice;
        void initIndice();

    public:
        int  numTypeTh;
        int  nbTh;
        bool qSave;
        bool qSlipt;
        bool qDeferred;
        bool qViewChrono;
        bool qInfo;
        std::string FileName;

        auto begin(); 
        auto end(); 
        auto size(); 

        int  getNbMaxThread();
        void init(int numType,int nbThread,bool qsaveInfo);
        void setFileName(std::string s);
        void setNbThread(int v);

        //BEGIN::multithread specx part
        template<class Function>
            Function run(Function myFunc);
                template<class Function>
                    Function sub_run_multithread(Function myFunc);
                template<class Function>
                    Function sub_run_async(Function myFunc);
                template<class Function>
                    Function sub_run_specx_R(Function myFunc);

                template<class Function>
                    std::vector<double> sub_run_multithread_beta(Function myFunc);
                template<class Function> 
                    std::vector<double> sub_run_specx(Function myFunc);
                template<class Function>
                    std::vector<double> sub_run_async_beta(Function myFunc);
        //END::multithread specx part

        //BEGIN::Detach part
        template<typename FctDetach>
                    auto sub_detach_future_alpha(FctDetach&& func) -> std::future<decltype(func())>;
        template<typename FctDetach>
                    auto sub_detach_future_beta(FctDetach&& func) -> std::future<decltype(func())>;
        template <typename result_type, typename FctDetach>
                    std::future<result_type> sub_detach_future_gamma(FctDetach func);

        template<class FunctionLambda,class FunctionLambdaDetach>
                    void sub_detach_specx_beta(FunctionLambda myFunc,int nbThreadsA,FunctionLambdaDetach myFuncDetach,int nbThreadsD);
        //END::Detach part

        //BEGIN::Thread affinity part
        template<class Function>
            void RunTaskInNumCPU(int idCPU,Function myFunc);
        template<class Function>
            void RunTaskInNumCPUs(const std::vector<int> & numCPU ,Function myFunc);
        //END::Thread affinity part

        template<class InputIterator, class Function>
            Function for_each(InputIterator first, InputIterator last,Function myFunc);

        TasksDispach(void);
        ~TasksDispach(void);
};


TasksDispach::TasksDispach() { 
    nbThTotal=std::thread::hardware_concurrency();
    nbTh=nbThTotal;
    numTypeTh=2; 
    qSave=false;
    FileName="TestDispach";
    qSlipt=false;
    qDeferred=false;
    qViewChrono=true;
    qViewChrono=false;
    qInfo=false;  
    initIndice();
}

TasksDispach::~TasksDispach(void) { 
}


auto TasksDispach::begin()
{
    return(indice.begin());
}

auto TasksDispach::end()
{
    return(indice.end());
}

auto TasksDispach::size()
{
    return(indice.size());
}

void TasksDispach::initIndice()
{
    indice.clear(); for (int i = 1; i <= nbTh; ++i)  { indice.push_back(i); }
}

void TasksDispach::setNbThread(int v)
{
    nbTh=std::min(v,nbThTotal);
    initIndice();
}

void TasksDispach::init(int numType,int nbThread,bool qsaveInfo)
{
    numTypeTh=numType; nbTh=nbThread; qSave=qsaveInfo; qInfo=false;  
    initIndice();
}


void TasksDispach::setFileName(std::string s)
{
    FileName=s;
}

int TasksDispach::getNbMaxThread()
{
    nbThTotal=std::thread::hardware_concurrency();
    return(nbThTotal);
}


template<class Function>
void TasksDispach::RunTaskInNumCPU(int idCPU,Function myFunc)
{
    const std::vector<int> v={idCPU};
    RunTaskInNumCPUs(v,myFunc);
}

template<class Function>
void TasksDispach::RunTaskInNumCPUs(const std::vector<int> & numCPU ,Function myFunc)
{
  int nbTh=numCPU.size();
  std::function<void()> func =myFunc;
  pthread_t thread_array[nbTh];
  pthread_attr_t pta_array[nbTh];

  for (int i = 0; i < nbTh; i++) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(numCPU[i], &cpuset);
    std::cout<<"Num CPU="<< numCPU[i] <<" activated"<<std::endl;
    pthread_attr_init(&pta_array[i]);
    pthread_attr_setaffinity_np(&pta_array[i], sizeof(cpuset), &cpuset);
    if (pthread_create(&thread_array[i],&pta_array[i],WorkerInNumCPU,&func)) { std::cerr << "Error in creating thread" << std::endl; }
  }

  for (int i = 0; i < nbTh; i++) {
        pthread_join(thread_array[i], NULL);
  }

  for (int i = 0; i < nbTh; i++) {
        pthread_attr_destroy(&pta_array[i]);
  }
}


template<typename FctDetach>
auto TasksDispach::sub_detach_future_alpha(FctDetach&& func) -> std::future<decltype(func())>
{
    using result_type = decltype(func());
    auto promise = std::promise<result_type>();
    auto future  = promise.get_future();
    std::thread(std::bind([=](std::promise<result_type>& promise)
    {
        try
        {
            promise.set_value(func()); 
        }
        catch(...)
        {
            promise.set_exception(std::current_exception());
        }
    }, std::move(promise))).detach();
    return future;
}

template<typename FctDetach>
auto TasksDispach::sub_detach_future_beta(FctDetach&& func) -> std::future<decltype(func())>
{
    auto task   = std::packaged_task<decltype(func())()>(std::forward<FctDetach>(func));
    auto future = task.get_future();
    std::thread(std::move(task)).detach();
    return std::move(future);
}


template <typename result_type, typename FctDetach>
    std::future<result_type> TasksDispach::sub_detach_future_gamma(FctDetach func) {
    std::promise<result_type> pro;
    std::future<result_type> fut = pro.get_future();
    std::thread([func](std::promise<result_type> p){p.set_value(func());},std::move(pro)).detach();
    return fut;
}

template<class FunctionLambda,class FunctionLambdaDetach>
void TasksDispach::sub_detach_specx_beta(FunctionLambda myFunc,int nbThreadsA,FunctionLambdaDetach myFuncDetach,int nbThreadsD)
{ 
    if (qInfo) { std::cout<<"Call Specx Detach="<<"\n"; }
    SpRuntime runtimeA(nbThreadsA);
    SpRuntime runtimeD(nbThreadsD);
    int idData=1;
    runtimeA.task(SpWrite(idData),
        [&,&runtimeD](int & depFakeData)
        {
            runtimeD.task([&,&depFakeData]()
            {
                myFuncDetach();
            });
            myFunc();
        });
    runtimeA.waitAllTasks();
    runtimeA.stopAllThreads();
    if (qSave)
    {
        runtimeA.generateDot("DetachA.dot", true);
        runtimeA.generateTrace("DetachA.svg");   
    }
    runtimeD.waitAllTasks();
    runtimeD.stopAllThreads();
    if (qSave)
    {
        runtimeD.generateDot("DetachD.dot", true);
        runtimeD.generateTrace("DetachD.svg");   
    }
    if (qInfo) { std::cout << std::endl; }
} 


template<class InputIterator, class Function>
Function TasksDispach::for_each(InputIterator first, InputIterator last,Function myFunc)
{
        if (numTypeTh==1) {
            std::vector< std::future< bool > > futures;
            for ( ; first!=last; ++first )
            { 
                auto const& idk = *first;
                if (qInfo) { std::cout<<"Call num Thread futures="<<idk<<"\n"; }
                if (qDeferred) { futures.emplace_back(std::async(std::launch::deferred,myFunc,idk)); }
                else { futures.emplace_back(std::async(std::launch::async,myFunc,idk)); }
            }
            for( auto& r : futures){ auto a =  r.get(); }
        }

        if (numTypeTh==2) {
            SpRuntime runtime(nbTh);  
            nbTh= runtime.getNbThreads();
            for ( ; first!=last; ++first )
            { 
                auto const& idk = *first;
                if (qInfo) { std::cout<<"Call num Thread Read Specx="<<idk<<"\n"; }
                runtime.task(SpRead(idk),myFunc).setTaskName("Op("+std::to_string(idk)+")");
                usleep(1);
                std::atomic_int counter(0);
            }
            runtime.waitAllTasks();
            runtime.stopAllThreads();
            if (qSave)
            {
                runtime.generateDot(FileName+".dot", true);
                runtime.generateTrace(FileName+".svg");   
            }
        }
        if (qInfo) { std::cout<<"\n"; }
    return myFunc;
}


template<class Function>
Function TasksDispach::sub_run_multithread(Function myFunc)
{
        auto begin = std::chrono::steady_clock::now();
        std::vector<std::thread> mythreads;
        for(int k= 0; k < nbTh; ++k){ 
            auto const& idk = k;
            if (qInfo) { std::cout<<"Call num Multithread ="<<k<<"\n"; }
            std::thread th(myFunc,idk);
            mythreads.push_back(move(th));
        }
        for (std::thread &t : mythreads) {
            t.join();
        }
        auto end = std::chrono::steady_clock::now();
        if (qInfo) { std::cout<<"\n"; }
        if (qViewChrono) {  std::cout << "===> Elapsed microseconds: "<< std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()<< " us\n"; std::cout<<"\n"; }
    return myFunc;
}

template<class Function>
std::vector<double> TasksDispach::sub_run_multithread_beta(Function myFunc)
{
        std::vector<double> valuesVec(nbTh,0);
        auto begin = std::chrono::steady_clock::now();
        auto LF=[&](const int& k) {  myFunc(k,valuesVec.at(k)); return true;};
        std::vector<std::thread> mythreads;
        for(int k= 0; k < nbTh; ++k){ 
            auto const& idk = k;
            if (qInfo) { std::cout<<"Call num Multithread ="<<k<<"\n"; }
            std::thread th(LF,idk);
            mythreads.push_back(move(th));
        }
        for (std::thread &t : mythreads) {
            t.join();
        }
        auto end = std::chrono::steady_clock::now();
        if (qInfo) { std::cout<<"\n"; }
        if (qViewChrono) {  std::cout << "===> Elapsed microseconds: "<< std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()<< " us\n"; std::cout<<"\n"; }
    return valuesVec;
}


template<class Function>
Function TasksDispach::sub_run_async(Function myFunc)
{
        auto begin = std::chrono::steady_clock::now();
        std::vector< std::future< bool > > futures;
        for(int k= 0; k < nbTh; ++k){ 
            auto const& idk = k;
            if (qInfo) { std::cout<<"Call num Thread futures="<<k<<"\n"; }
            if (qDeferred) { futures.emplace_back(std::async(std::launch::deferred,myFunc,idk)); }
            else { futures.emplace_back(std::async(std::launch::async,myFunc,idk)); }
        }
        for( auto& r : futures){ auto a =  r.get(); }
        auto end = std::chrono::steady_clock::now();
        if (qInfo) { std::cout<<"\n"; }
        if (qViewChrono) {  std::cout << "===> Elapsed microseconds: "<< std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()<< " us\n";std::cout<<"\n"; }
    return myFunc;
}

template<class Function>
std::vector<double> TasksDispach::sub_run_async_beta(Function myFunc)
{
        std::vector<double> valuesVec(nbTh,0);
        auto begin = std::chrono::steady_clock::now();
        auto LF=[&](const int& k) {  myFunc(k,valuesVec.at(k)); return true;};
        std::vector< std::future< bool > > futures;
        for(int k= 0; k < nbTh; ++k){ 
            auto const& idk = k;
            if (qInfo) { std::cout<<"Call num Thread futures="<<k<<"\n"; }
            futures.emplace_back(std::async(std::launch::async,LF,idk)); 
        }
        for( auto& r : futures){ auto a =  r.get(); }
        auto end = std::chrono::steady_clock::now();
        if (qInfo) { std::cout<<"\n"; }
        if (qViewChrono) {  std::cout << "===> Elapsed microseconds: "<< std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()<< " us\n";std::cout<<"\n"; }
    return valuesVec;
}




template<class Function>
std::vector<double> TasksDispach::sub_run_specx(Function myFunc)
{
        std::vector<double> valuesVec(nbTh,0); 
        SpRuntime runtime(nbTh);  
        auto begin = std::chrono::steady_clock::now();
        nbTh= runtime.getNbThreads();
        for(int k= 0; k < nbTh; ++k)
        { 
            auto const& idk = k;
            if (qInfo) { std::cout<<"Call num Thread Write Specx="<<idk<<"\n"; }
            runtime.task(SpRead(idk),SpWrite(valuesVec.at(idk)),myFunc).setTaskName("Op("+std::to_string(idk)+")");
            usleep(0); std::atomic_int counter(0);
        }
        runtime.waitAllTasks();
        runtime.stopAllThreads();
        auto end = std::chrono::steady_clock::now();
        if (qSave)
        {
            runtime.generateDot(FileName+".dot", true);
            runtime.generateTrace(FileName+".svg");   
        }
        if (qInfo) { std::cout<<"\n"; }
        if (qViewChrono) {  std::cout << "===> Elapsed microseconds: "<< std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()<< " us\n";std::cout<<"\n"; }

    return valuesVec;
}


template<class Function>
Function TasksDispach::sub_run_specx_R(Function myFunc)
{
        SpRuntime runtime(nbTh);  
        auto begin = std::chrono::steady_clock::now();
        nbTh= runtime.getNbThreads();
        int iValue=0;
        for(int k= 0; k < nbTh; ++k)
        { 
            auto const& idk = k;
            if (qInfo) { std::cout<<"Call num Thread Read Specx="<<idk<<"\n"; }
            runtime.task(SpRead(idk),myFunc).setTaskName("Op("+std::to_string(idk)+")");
            usleep(0); std::atomic_int counter(0);
        }
        runtime.waitAllTasks();
        runtime.stopAllThreads();
        auto end = std::chrono::steady_clock::now();
        if (qSave)
        {
            runtime.generateDot(FileName+".dot", true);
            runtime.generateTrace(FileName+".svg");   
        }
        
        if (qInfo) { std::cout<<"\n"; }
        if (qViewChrono) {  std::cout << "===> Elapsed microseconds: "<< std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()<< " us\n"; std::cout<<"\n"; }
    return myFunc;
}




template<class Function>
Function TasksDispach::run(Function myFunc)
{
  switch(numTypeTh) {
    case 1: return(sub_run_async(myFunc));
    break;
    case 2: return(sub_run_specx_R(myFunc));
    break;
    default:
        return(sub_run_multithread(myFunc));
  }
}



//=======================================================================================================================
//...
//=======================================================================================================================




class TasksDispachComplex 
{
    int nbThTotal;   
    std::string FileName;

    public:
        int nbTh;
        bool qViewChrono;
        bool qInfo;
        bool qSave;

        void setNbThread(int v);
        int  getNbMaxThread();

        template <typename ... Ts>
            auto parameters(Ts && ... ts);

        template <typename ... Ts>
            void runTask( Ts && ... ts );

        template <typename ... Ts>
            void runTaskLoopAsync( Ts && ... ts );

        template <typename ... Ts>
            void runTaskLoopSpecx( Ts && ... ts );
        
        TasksDispachComplex(void);

        void setFileName(std::string s);
};


TasksDispachComplex::TasksDispachComplex()
{
    nbThTotal=std::thread::hardware_concurrency();
    nbTh=nbThTotal;
    qViewChrono=false;
    qInfo=false;
    qSave=false;
    FileName="TestDispachComplex";
}

void TasksDispachComplex::setFileName(std::string s)
{
    FileName=s;
}

void TasksDispachComplex::setNbThread(int v)
{
    nbTh=std::min(v,nbThTotal);
}

int TasksDispachComplex::getNbMaxThread()
{
    nbThTotal=std::thread::hardware_concurrency();
    return(nbThTotal);
}

template <typename ... Ts>
auto TasksDispachComplex::parameters(Ts && ... ts)
{
    return std::forward_as_tuple( std::forward<Ts>(ts)... );
}

template <typename ... Ts>
void TasksDispachComplex::runTask( Ts && ... ts )
{
    auto begin = std::chrono::steady_clock::now();
    auto args = NA::make_arguments( std::forward<Ts>(ts)... );
    auto && task = args.get(_task);
    auto && parameters = args.get_else(_parameters,std::make_tuple());
    Backend::Runtime runtime;
    auto tp=std::tuple_cat( 
					Backend::makeSpData( parameters ), 
					std::make_tuple( task ) 
				);
    std::apply( [&runtime](auto... args){ runtime.task(args...); }, tp );
    auto end = std::chrono::steady_clock::now();
    if (qViewChrono) {  std::cout << "===> Elapsed microseconds: "<< std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()<< " us\n"; std::cout<<"\n"; }
}

template <typename ... Ts>
void TasksDispachComplex::runTaskLoopAsync( Ts && ... ts )
{
    auto begin = std::chrono::steady_clock::now();
    auto args = NA::make_arguments( std::forward<Ts>(ts)... );
    auto && task = args.get(_task);
    auto && parameters = args.get_else(_parameters,std::make_tuple());
    Backend::Runtime runtime;
    std::vector< std::future< bool > > futures;
    std::cout <<"nbTh="<<nbTh<< "\n";

    for (int k = 0; k < nbTh; k++) {
        if (qInfo) { std::cout<<"Call num Thread futures="<<k<<"\n"; }
		auto tp=std::tuple_cat( 
					Backend::makeSpData( parameters ), 
					std::make_tuple( task ) 
				);
		auto LamdaTransfert = [&]() {
			std::apply([&runtime](auto... args){ runtime.task(args...); }, tp);
            return true; 
		};

        /*
        futures.emplace_back(
            std::async(std::launch::async,LamdaTransfert));
        */

        futures.emplace_back(
            std::async(std::launch::deferred,LamdaTransfert));
    }
    for( auto& r : futures){ auto a =  r.get(); }
    auto end = std::chrono::steady_clock::now();
    if (qViewChrono) {  std::cout << "===> Elapsed microseconds: "<< std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()<< " us\n"; std::cout<<"\n"; }
}

/*
template <typename ... Ts>
void TasksDispachComplex::runTaskLoopMultithread( Ts && ... ts )
{
    //add mutex
    auto begin = std::chrono::steady_clock::now();
    auto args = NA::make_arguments( std::forward<Ts>(ts)... );
    auto && task = args.get(_task);
    auto && parameters = args.get_else(_parameters,std::make_tuple());
    Backend::Runtime runtime;
    std::vector<std::thread> mythreads;
    std::cout <<"nbTh="<<nbTh<< "\n";
    std::mutex mtx; 
    for (int k = 0; k < nbTh; k++) {
        std::cout<<"Call num Thread futures="<<k<<"\n";
		auto tp=std::tuple_cat( 
					Backend::makeSpData( parameters ), 
					std::make_tuple( task ) 
				);
		auto LamdaTransfert = [&]() {
			std::apply([&runtime](auto... args){ runtime.task(args...); }, tp);
            return true; 
		};
        std::thread th(LamdaTransfert);
        mythreads.push_back(move(th));
    }
    for (std::thread &t : mythreads) { t.join();}
    auto end = std::chrono::steady_clock::now();
    if (qViewChrono) {  std::cout << "===> Elapsed microseconds: "<< std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()<< " us\n"; std::cout<<"\n"; }
}
*/


template <typename ... Ts>
void TasksDispachComplex::runTaskLoopSpecx( Ts && ... ts )
{
    auto begin = std::chrono::steady_clock::now();
    auto args = NA::make_arguments( std::forward<Ts>(ts)... );
    auto && task = args.get(_task);
    auto && parameters = args.get_else(_parameters,std::make_tuple());
    Backend::Runtime runtime;
    std::cout <<"nbTh="<<nbTh<< "\n";

    SpRuntime runtime_Specx(nbTh); 


    auto tpSpecxFront=Backend::makeSpData( parameters );
        int NbtpSpecxFront=std::tuple_size<decltype(tpSpecxFront)>::value;
        std::cout <<"Size Tuple Front="<<NbtpSpecxFront<< std::endl;

    auto LambdaExpression=std::make_tuple( task );
        int NbLambdaExpression=std::tuple_size<decltype(LambdaExpression)>::value;
        std::cout <<"Size Tuple Parameters="<<NbLambdaExpression<< std::endl;
        
    //auto tpBackend=Backend::makeSpData( parameters );
    auto tpBackend=Backend::makeSpDataSpecx( parameters );
         int NbtpBackend=std::tuple_size<decltype(tpBackend)>::value;
        std::cout <<"Size tpBackend="<<NbtpBackend<< std::endl;

    auto tpSpecx=std::tuple_cat( 
					Backend::makeSpDataSpecx( parameters ), 
					std::make_tuple( task ) 
				);

    for (int k = 0; k < nbTh; k++) {
        if (qInfo) { std::cout<<"Call num Thread specx="<<k<<"\n"; }
        std::apply([&](auto &&... args) { runtime_Specx.task(args...); },tpSpecx);
    }

    runtime_Specx.waitAllTasks();
    runtime_Specx.stopAllThreads();

    auto end = std::chrono::steady_clock::now();

    if (qSave)
    {
        runtime_Specx.generateDot("Test.dot", true);
        runtime_Specx.generateTrace("Test.svg");   
    }

    if (qViewChrono) {  std::cout << "===> Elapsed microseconds: "<< std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()<< " us\n"; std::cout<<"\n"; }
}



//=======================================================================================================================
//...
//=======================================================================================================================


void testScanAllThreadMethods()
{
    int qPutLittleTroublemaker=true;
    int time_sleep= 100000;
    auto P001=[time_sleep,qPutLittleTroublemaker](const int i,double& s) {  
            double sum=0.0; 
            for(int j=0;j<100;j++) { sum+=double(j); }
            if (qPutLittleTroublemaker) {
                srand((unsigned)time(0)); int time_sleep2 = rand() % 5000 + 1; usleep(time_sleep2); 
            }
            usleep(time_sleep);
            s=sum+i;      
        return true;
    };

    bool qChrono=false;

    TasksDispach Fg1; 
    int nbThreads = Fg1.nbTh;
    //int nbThreads = 96;
    Color(7); std::cout<<"Test Scan [";
    Color(3); std::cout<<nbThreads;
    Color(7); std::cout<<"] Thread Methods >>> ";
    Fg1.setFileName("Results"); 
    Fg1.init(1,nbThreads,true); Fg1.qViewChrono=qChrono; 
    std::vector<double> valuesVec1=Fg1.sub_run_multithread_beta(P001);
    double Value1=std::reduce(valuesVec1.begin(),valuesVec1.end()); 

    TasksDispach Fg2; 
    Fg2.setFileName("Results"); 
    Fg2.init(1,nbThreads,true); Fg2.qViewChrono=qChrono; 
    std::vector<double> valuesVec2=Fg2.sub_run_async_beta(P001);
    double Value2=std::reduce(valuesVec2.begin(),valuesVec2.end()); 

    TasksDispach Fg3; 
    Fg3.setFileName("Results"); 
    Fg3.init(2,nbThreads,true); Fg3.qViewChrono=qChrono; 
    std::vector<double> valuesVec3=Fg3.sub_run_specx(P001);
    double Value3=std::reduce(valuesVec3.begin(),valuesVec3.end()); 
    if ((Value1==Value2) && (Value1==Value3)) {
        Color(2); std::cout <<"OK"<< "\n"; 
    } 
    else 
    {
        Color(1); std::cout <<"ERROR "<<"m1:"<<Value1<<" m2:"<<Value2<<" m3:"<<Value3<< "\n"; 
    }
    std::cout << "\n"; 
    Color(7);
    std::cout << "\n"; 
}



