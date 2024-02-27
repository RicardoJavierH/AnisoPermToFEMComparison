//
// Created by victor on 12/09/2022.
//

#include "InputTreatment.h"
#include "MeshInit.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>
#include "pzvisualmatrix.h"
#include <TPZTimer.h>
#include "computStatist.h"
#include "TPZCreateHybridizedMixedSpaces.h"
using namespace std;

#ifdef FEMCOMPARISON_TIMER
//#include <chrono>
#include <math.h>//used only for timer statistics
#include <fstream>//used only for timer statistics
#include <vector>//used only for timer statistics

#endif
int main(int argc, char *argv[]) {
#ifdef FEMCOMPARISON_TIMER
    //TPZTimer timer;
    //timer.start();

    bool atypical1=false;

    bool MKL_contribute;
#ifdef FEMCOMPARISON_USING_MKL
    MKL_contribute=true;
#endif
#endif

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    PreConfig pConfig;
    pConfig.k = 1;//
    pConfig.n = 2;
    pConfig.problem = "ESinSin";              //// {"ESinSin","EArcTan",ESteklovNonConst", "ESteepWave"}
    pConfig.approx = "Mixed";                //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";       //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel =1;                     //// How many refinements
    pConfig.postProcess = false;                  //// Print geometric and computational mesh
    pConfig.shouldColor =false;
    pConfig.isTBB = false;
    pConfig.tData.nThreads = 0;


    if(argc == 2) {
        argc = 1;
        pConfig.tData.nThreads = atoi(argv[1]);
    }
    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);

    pConfig.exp *= pow(2,pConfig.refLevel-1);
    pConfig.h = 1./pConfig.exp;
    ProblemConfig config;
    Configure(config,pConfig.refLevel,pConfig,argv);

    std::set<int> matids = {1};
    std::set<int> bcmatids ={-1,-2};
    TPZCreateHybridizedMixedSpaces createHMixed(config.gmesh,matids,bcmatids);
    TPZMultiphysicsCompMesh* hybridizedMixed = createHMixed.GenerateMesh();


    return 0;
}


