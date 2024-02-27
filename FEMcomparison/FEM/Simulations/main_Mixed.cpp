#include "InputTreatment.h"
#include "MeshInit.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>
#include <TPZTimer.h>
using namespace std;
#ifdef FEMCOMPARISON_TIMER
bool ContributeVOL=true;
bool ContributeBC=true;
double solveTime =0.;
double assembleTime =0.;
extern double calcstiffTime;
extern double contributeTime; //  Total contribute time
double contributeTimeVol=0.;
double contributeTimeBoundary=0.;
double contributeTimeInterface=0.;
double interfaceTime=0.;
extern int64_t contributeCounter;
int64_t contributeMaterialCounter=0;
int64_t contributeBoundaryCounter=0;
double solveglobaltime;
#endif

int main(int argc, char *argv[]) {
#ifdef FEMCOMPARISON_TIMER
    TPZTimer timer;
    timer.start();
    calcstiffTime=0.;
    contributeTime =0;
    contributeCounter=0;
#endif
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    PreConfig pConfig;
    pConfig.k = 1;//
    pConfig.n = 2;
    pConfig.problem = "ESinSin";        //// {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfig.approx = "Mixed";                    //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";          //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel = 2;                        //// How many refinements
    pConfig.postProcess = false;                    //// Print geometric and computational mesh

    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);

    pConfig.exp *= pow(2,pConfig.refLevel-1);
    pConfig.h = 1./pConfig.exp;
    ProblemConfig config;
    Configure(config,pConfig.refLevel ,pConfig,argv);
    Solve(config,pConfig);
    pConfig.hLog = pConfig.h;
    
    if(pConfig.postProcess)
    {
        std::string command = "cp ErroMixed.txt " + pConfig.plotfile + "/Erro.txt";
        system(command.c_str());
        FlushTable(pConfig,argv);
    }
    std::cout<<"solveglobaltime= "<< solveglobaltime << std::endl;
    std::cout<<"solveTime= "<< solveTime << std::endl;
    std::cout<<"assembleTime= "<< assembleTime << std::endl;
    std::cout<<"calcstiffTime= "<< calcstiffTime << std::endl;
    std::cout<<"contributeTimeInterface= "<< contributeTimeInterface << std::endl;
    std::cout<<"contributeTimeVol= "<< contributeTimeVol << std::endl;
    std::cout<<"contributeTimeBoundary= "<< contributeTimeBoundary << std::endl;
    //#ifdef OPTIMIZE_RUN_TIME
   std::cout<<"contributeTime= "<< contributeTime << std::endl;
//#endif
    //cout<<"contadorTimeMaterial= "<<contadorTimeMaterial<<endl;
    cout<<"contributeCounter= "<<contributeCounter<<endl;
    cout<<"contributeMaterialCounter= "<<contributeMaterialCounter<<endl;
    cout<<"contributeBoundaryCounter= "<<contributeBoundaryCounter<<endl;

    return 0.;
}







