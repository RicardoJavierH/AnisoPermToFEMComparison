//
//  LCC_LagrangeMultiplier.h
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//

#ifndef __PZ__LCC_LagrangeMultiplier__
#define __PZ__LCC_LagrangeMultiplier__

#include <iostream>
#include "TPZLagrangeMultiplierCS.h"
#include "pzfmatrix.h"
#include "TPZMaterialDataT.h"
/// Material which implements a Lagrange Multiplier
class LCC_LagrangeMultiplier : public TPZLagrangeMultiplierCS<STATE>
{

    /** @brief Pointer to a blocked matrix object*/
    TPZManVector<TPZFNMatrix<1000, STATE>,8> fStiffnessMatrix;
    
    public :
	/** @brief Simple constructor */
    LCC_LagrangeMultiplier() : TPZRegisterClassId(&LCC_LagrangeMultiplier::ClassId),
    TPZLagrangeMultiplierCS<STATE>()
    {
        fStiffnessMatrix.Resize(14);
        for(int iMatrix = 0; iMatrix < 14; iMatrix ++){
            fStiffnessMatrix[iMatrix].Resize(0,0);
        }
    }
	/** @brief Constructor with the index of the material object within the vector */
    LCC_LagrangeMultiplier(int nummat, int dimension, int nstate) : TPZRegisterClassId(&LCC_LagrangeMultiplier::ClassId),
    TPZLagrangeMultiplierCS<STATE>(nummat,dimension,nstate)
    {
        fStiffnessMatrix.Resize(14);
        for(int iMatrix = 0; iMatrix < 14; iMatrix ++){
            fStiffnessMatrix[iMatrix].Resize(0,0);
        }
    }
	
	/** @brief Copy constructor */
	LCC_LagrangeMultiplier(const LCC_LagrangeMultiplier &copy) : TPZRegisterClassId(&LCC_LagrangeMultiplier::ClassId),
    TPZLagrangeMultiplierCS<STATE>(copy)
    {
	    if(fStiffnessMatrix.size() != 14){
	        DebugStop();
	    }
        for(int iMatrix = 0; iMatrix < 14; iMatrix ++){
            fStiffnessMatrix[iMatrix] = copy.fStiffnessMatrix[iMatrix];
        }
    }
    
    LCC_LagrangeMultiplier &operator=(const LCC_LagrangeMultiplier &copy)
    {
        TPZLagrangeMultiplierCS<STATE>::operator=(copy);
        fNStateVariables = copy.fNStateVariables;
        fDimension = copy.fDimension;
        fMultiplier = copy.fMultiplier;
        if(fStiffnessMatrix.size() != 14){
            DebugStop();
        }
        for(int iMatrix = 0; iMatrix < 14; iMatrix ++){
            fStiffnessMatrix[iMatrix] = copy.fStiffnessMatrix[iMatrix];
        }
        return *this;
    }
    
    TPZMaterial *NewMaterial() const override
    {
        return new LCC_LagrangeMultiplier(*this);
    }
    
	/** @brief Destructor */
	virtual ~LCC_LagrangeMultiplier()
    {
        
    }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const override
    {
        return fDimension;
    }
    
    virtual void SetMultiplier(STATE mult) override
    {
        fMultiplier = mult;
    }
    
    STATE Multiplier()
    {
        return fMultiplier;
    }
    
    int NStateVariables()
    {
        return fNStateVariables;
    }
	
	virtual std::string Name() const override
    {
        return "LCC_LagrangeMultiplier";
    }

    void GetStiffnessMatrix(TPZFMatrix<STATE> &S, int iMatrix){
        if(iMatrix <0 || iMatrix >13){
            DebugStop();
        }
        S = fStiffnessMatrix[iMatrix];
    }

    void FillStiffnessMatrix(const TPZFMatrix<STATE> &S,int iMatrix){
        if(iMatrix <0 || iMatrix >13){
            DebugStop();
        }
        fStiffnessMatrix[iMatrix] = S;
    }
	
    // print the data in human readable form
    virtual void Print(std::ostream &out) const override;
	/**
	 * @brief Fill material data parameter with necessary requirements for the ContributeInterface method.
     * @since April 10, 2007
	 */
	/**
	 * Here, in base class, all requirements are considered as necessary. \n
	 * Each derived class may optimize performance by selecting only the necessary data.
	 */
//	virtual void FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data) const override
//    {
//        data.SetAllRequirements(false);
//    }
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of interface elements */
    virtual void FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavec_left, std::map<int, TPZMaterialDataT<STATE>> &datavec_right) override
    {
        data.SetAllRequirements(false);
//        data.fNeedsNormal = true;
        for(auto &it : datavec_left)
        {
            it.second.SetAllRequirements(false);
        }
        for(auto &it : datavec_right)
        {
            it.second.SetAllRequirements(false);
        }
        for(auto it = datavec_left.begin() ; it != datavec_left.end() ; it++){
            it->second.fNeedsSol = true;
        }
        for(auto it = datavec_right.begin() ; it != datavec_right.end() ; it++){
            it->second.fNeedsSol = true;
        }
    }

	virtual int NStateVariables() const override
    {
        return fNStateVariables;
    }
	
		
	
    /** @{
     * @name Save and Load methods
     */
    
	/** @brief Unique identifier for serialization purposes */
	public:
virtual int ClassId() const override;

	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context) override;
	
    /**
     * @}
     */
};

#endif /* defined(__PZ__LCC_LagrangeMultiplier__) */
