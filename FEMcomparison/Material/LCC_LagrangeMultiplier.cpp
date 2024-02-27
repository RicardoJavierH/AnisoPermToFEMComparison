//
//  LCC_LagrangeMultiplier.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//
#include "fstream"
#include "LCC_LagrangeMultiplier.h"
#include "pzaxestools.h"
#ifdef FEMCOMPARISON_USING_MKL
#include "mkl.h"
#endif
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("LagrangeMultipliersData"));
static LoggerPtr logerror(Logger::getLogger("LagrangeMultipliersError"));
#endif
#include "TPZTimer.h"
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger loggerCTI("contributeTimeInterface");
#endif


/** @brief Unique identifier for serialization purposes */
int LCC_LagrangeMultiplier::ClassId() const{
    return Hash("LCC_LagrangeMultiplier") ^ TPZMaterial::ClassId() << 1;
}

/** @brief Saves the element data to a stream */
void LCC_LagrangeMultiplier::Write(TPZStream &buf, int withclassid) const
{
    TPZMaterial::Write(buf, withclassid);
    buf.Write(&fNStateVariables);
}

/** @brief Reads the element data from a stream */
void LCC_LagrangeMultiplier::Read(TPZStream &buf, void *context)
{
    TPZMaterial::Read(buf, context);
    buf.Read(&fNStateVariables);
    
}

/**
 * @brief Computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation
 * @param data [in]
 * @param dataleft [in]
 * @param dataright [in]
 * @param weight [in]
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @since June 5, 2012
 */



// print the data in human readable form
void LCC_LagrangeMultiplier::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZLagrangeMultiplierCS<STATE>::Print(out);
    out << "NStateVariables " << this->fNStateVariables << std::endl;
    out << "fDimension " << this->fDimension << std::endl;
    out << "fMultiplier " << this->fMultiplier << std::endl;
}

