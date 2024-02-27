
//  ProblemConfig.h
//  ErrorEstimation
//
//  Created by Philippe Devloo on 09/05/18.
//

#ifndef ProblemConfig_h
#define ProblemConfig_h

#include <set>
#include "TPZAnalyticSolution.h"

struct PreConfig;

/// class to guide the error estimator
struct ProblemConfig
{
    /// geometric mesh on which the computational meshes are based
    TPZGeoMesh *gmesh = 0;
    /// polynomial order of the original mesh
    int k = -1;
    int n = -1;
    
    /// option to compute the error based on continuous pressures or not
    bool makepressurecontinuous = true;
    
    /// number of uniform refinements applied to the mesh
    int ndivisions = -1;
    int adaptivityStep = -1;
    int dimension = 0;
    bool prefine = false;
    bool steklovexample = false;
    bool GalvisExample = false;
    bool TensorNonConst = false;
    bool MeshNonConvex = false;
    
    STATE alpha=1;
    STATE Km = 0.;
    STATE coefG = 0.;
    /// directory where the files will be stored
    std::string dir_name = ".";
    /// name identifying the problem
    std::string problemname;
    /// set of materialids in the mesh
    std::set<int> materialids;
    /// set of boundary condition material ids
    std::set<int> bcmaterialids;
    /// exact solution

    TPZAutoPointer<TLaplaceExample1> exact;

    TPZFMatrix<STATE> perm, invperm;
    
    ProblemConfig() = default;

    ProblemConfig(const ProblemConfig &cp) = default;

    ProblemConfig &operator=(const ProblemConfig &cp) = default;
    
    void InitializePermTensor(int dim, PreConfig &pConfig);
};

struct MultithreadData{

    int nThreads;
    int maxThreads, maxRef;
    unsigned long int assembleTime;
    unsigned long int solveTime;

};

struct Statistics{
    int nLoops;
    int iterNum;
    TPZFMatrix<double> timeVec;
    TPZVec<double> avg, spu, cvar, serialTime;
    std::vector<std::ofstream*> csv, txt;
};

struct RSimulation {
    int nThreads, maxThreads, iterNum, maxRef, nRef, nDof;
    double assembleTime, solveTime, totalTime, L2Error,energyError;
};

struct Target{
    bool automated = false;
    bool timeEfficiency = false;
    bool errorMeasurement = false;
};

struct PreConfig{
    std::ofstream Erro, timer;
    std::ofstream *speedUpOfstream, *speeedUpPath;

    TPZVec<REAL> *rate, *Log;
    std::vector<double> errorRate, errorLog;
    int refLevel = -1;

    int k = 1;
    int n = 1;
    int dim = 1;
    int topologyMode = -1;

    std::string problem;
    std::string approx;
    std::string topology;           //Topology' name typed as input
    std::string topologyFileName;   //Simplified name used for naming files/directories

    REAL perm_Q1 = 5;      /// Permeability coefficient of even quadrants (Steklov only)
    REAL perm_Q2 = 1;

    REAL hLog = -1, h = -1000, href0 = 1;
    int numErrors = 4;

    std::string plotfile;
    std::string speedUpFilePath;
    std::string automatedFileName;
    std::string automatedFilePath;
    
    int mode = -1;           // 0 = "H1"; 1 = "Hybrid"; 2 = "Mixed";
    int argc = 1;
    int type= -1;

    bool makeScript = false;
    bool isSimplify = false;
    bool postProcess = true;
    int exp = 2; // Initial exponent of mesh refinement (numElem = 2*2^exp)
    
    bool shouldColor = true;
    bool isTBB = true;
    
    MultithreadData tData;
    Statistics stat;
    RSimulation rSimulation;
    Target target;
    
    int ref2D = -1;
    int ref3D = -1;
    int ref0 = -1;
};

inline void ProblemConfig::InitializePermTensor(int dim, PreConfig &pConfig)
{
    perm = TPZFMatrix<STATE>(3, 3, 0);
    bool isIsotropic = (pConfig.type != 4) ? true : false;
    if (isIsotropic ){
        if(dim==2){
        perm = TPZFMatrix<STATE>({{1.,0.},{0.,1.}});
        }
        else if(dim == 3){
        perm = TPZFMatrix<STATE>({{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}});
        }
        else DebugStop();
        return;
    }
    else
    {
        if(dim == 2 ){
            TPZVec<TPZFMatrix<STATE>>permVec(11);
            {
                STATE theta = (25./180.)*M_PI;
                TPZFMatrix<STATE> rot({{std::cos(theta), -std::sin(theta)},
                    {std::sin(theta), std::cos(theta)}});
                TPZFMatrix<STATE> invRot(2,2,0.);
                {
                    const STATE det = rot(0,0)*rot(1,1)-rot(1,0)*rot(0,1);
                    invRot(0,0) = rot(1,1)/det;
                    invRot(1,1) = rot(0,0)/det;
                    invRot(1,0) = -rot(1,0)/det;
                    invRot(0,1) = -rot(0,1)/det;
                }
                
                
                const STATE bigEigen = 0.3;
                for(int ialpha=0; ialpha<11; ialpha++){
                    TPZFMatrix<STATE> aux(2,2,0.),sol(2,2,0), diag(2,2,0);
                    const STATE alpha = bigEigen/(std::pow(2,ialpha));
                    diag(0,0) = bigEigen;
                    diag(1,1) = alpha;
                    
                    rot.Multiply(diag,aux);
                    aux.Multiply(invRot,sol);
                    
                    permVec[ialpha] = TPZFMatrix<STATE>(sol);
                }
                
            }
            
            //perm = TPZFMatrix<STATE>({{0.24643600207249652,0.11486836424569073},
                                     // {0.11486836424569073,0.05366399792750343},
                //});
           
            perm = permVec[1];
            
            TPZFNMatrix<100, STATE> copy(perm);
            invperm.Redim(copy.Rows(),copy.Rows());
            int64_t r;
            for(r=0; r<copy.Rows(); r++){
                invperm(r,r) = 1.;
            }
            copy.Solve_LU(&invperm);
        }
        else if(dim==3){
            TPZFMatrix<STATE> P = {{-1./3.,2./3.,-2./3.},
                {-2./3.,1./3.,2./3.},
                {-2./3.,-2./3.,-1./3.}};
            
            TPZFMatrix<STATE> PInv(3,3);
            PInv.Identity();
            P.Solve_LU(&PInv);
            invperm = PInv;
            
            TPZFMatrix<STATE> D = {{0.3,0.,0.},{0.,0.001,0.},{0.,0.,0.00001}};
            TPZFMatrix<STATE> Ktemp(3,3,0),K(3,3,0);;

            P.SetIsDecomposed(ENoDecompose);//Just to fix a mistake
            P.Multiply(D, Ktemp);
            Ktemp.Multiply(PInv, K);
            //perm = K;
            perm = TPZFMatrix<STATE>({{0.04266666666666666, 0.07066666666666666,
                0.05799999999999998}, {0.07066666666666666, 0.13599999999999995,
                0.12866666666666662}, {0.05799999999999998, 0.12866666666666662,
                0.1423333333333333}});
        }
        else DebugStop();
    }
    
    this->exact->SetPermeabilyTensor(perm,invperm);
}

#endif /* ProblemConfig_h */
