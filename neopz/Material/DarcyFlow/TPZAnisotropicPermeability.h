//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZANISOTROPICPERMEABILITY_H
#define TPZANISOTROPICPERMEABILITY_H

#include <functional>
#include "pzreal.h"
#include "pzvec.h"
#include "pzfmatrix.h"

// Alias to improve readability of the permeability function type
using PermTensorFunctionType = std::function<TPZFMatrix<STATE>(const TPZVec<REAL> &coord)>;

// Forward declaration of dummy BC interface class
class TPZAnisotropicPermeabilityBC;

/**
 * @brief  This class implements the interface with the methods required to
 * handle the permeability field of an isotropic material.
 */
class TPZAnisotropicPermeability : public virtual TPZSavable {

public:
    using TInterfaceBC = TPZAnisotropicPermeabilityBC;

    TPZAnisotropicPermeability() : fHeterogenPermeability({{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}}), fPermeabilityFunction(NULL){
    }
    
    TPZAnisotropicPermeability(const TPZAnisotropicPermeability &copy) = default;
    TPZAnisotropicPermeability &operator=(const TPZAnisotropicPermeability &copy) = default;
    /**
     * @brief Set a constant permeability to the material
     * @param [in] constant permeability value
     */
    void SetPermeability( TPZFMatrix<STATE> &permT);

    /**
     * @brief Set a varying permeability field to the material
     * @param [in] perm_function function that describes the permeability field
     */
    void SetPermeabilityFunction(PermTensorFunctionType &perm_function);

    /**
     * @brief Return the permeability value at a coordinate
     * @param [in] coord coordinate of interest
     */
    TPZFMatrix<STATE> GetPermeability(const TPZVec<REAL> &coord);
    
    void GetPermeability(const TPZVec<REAL> &coord, TPZFMatrix<STATE> &perm, TPZFMatrix<STATE> &invperm);


    [[nodiscard]] int ClassId() const override;

    void Read(TPZStream &buf, void *context) override {};

    void Write(TPZStream &buf, int withclassid) const override {};

private:

    // Member variable to describe a constant permeability field
    TPZFMatrix<STATE> fHeterogenPermeability, fInvPerm;

    // Member variable to describe a varying permeability field
    PermTensorFunctionType fPermeabilityFunction = NULL;
};

// Dummy BC interface class
class TPZMaterial;
class TPZAnisotropicPermeabilityBC : public TPZAnisotropicPermeability {
protected:
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    static void SetMaterialImpl(TPZMaterial *) {}
};

#endif //TPZISOTROPICPERMEABILITY_H
