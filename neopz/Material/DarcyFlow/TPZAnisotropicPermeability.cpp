//
// Created by Gustavo Batistela on 5/13/21.
//
#include "pzfmatrix.h"

#include "TPZAnisotropicPermeability.h"

void TPZAnisotropicPermeability::SetPermeability( TPZFMatrix<STATE> &permT) {
    fHeterogenPermeability.Resize(permT.Rows(),permT.Cols());
    
    for(int irow=0; irow<permT.Rows(); irow++){
        for(int icol=0; icol<permT.Cols(); icol++){
            fHeterogenPermeability(irow,icol) = permT(irow,icol);
        }
    }
    
    //fHeterogenPermeability = permT;
    TPZFNMatrix<100, STATE> copy(permT);
    fInvPerm.Redim(copy.Rows(),copy.Rows());
    int64_t r;
    for(r=0; r<copy.Rows(); r++){
        fInvPerm(r,r) = 1.;
    }
    
    copy.Solve_LU(&fInvPerm);
    const int spacedim = 3;
    fHeterogenPermeability.Resize(spacedim,spacedim);
    fInvPerm.Resize(spacedim,spacedim);
}

void TPZAnisotropicPermeability::SetPermeabilityFunction(PermTensorFunctionType &perm_function) {
    fPermeabilityFunction = perm_function;
}

TPZFMatrix<STATE> TPZAnisotropicPermeability::GetPermeability(const TPZVec<REAL> &coord)
{
    return fPermeabilityFunction ? fPermeabilityFunction(coord) : fHeterogenPermeability;
}

void TPZAnisotropicPermeability::GetPermeability(const TPZVec<REAL> &coord, TPZFMatrix<STATE> &perm, TPZFMatrix<STATE> &invperm){
    perm = fHeterogenPermeability;
    
    invperm = fInvPerm;
    
}

int TPZAnisotropicPermeability::ClassId() const {
    return Hash("TPZAnisotropicPermeability");
}
