//
// Created by victor on 12/09/2022.
//

#include "TPZCreateHybridizedMixedSpaces.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include "pzgeoelside.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZLagrangeMultiplierCS.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlowAniso.h"
#include "LCC_MultiphysicsInterfaceEl.h"
#include <sstream>
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZCompMeshTools.h"
#include "pzintel.h"


TPZCreateHybridizedMixedSpaces::TPZCreateHybridizedMixedSpaces(TPZGeoMesh *gmesh, std::set<int> &matids, std::set<int> &bcmatIds){
DefaultConstructor(gmesh,matids,bcmatIds);
ComputePeriferalMaterialIds();
}

TPZCreateHybridizedMixedSpaces::TPZCreateHybridizedMixedSpaces(TPZGeoMesh *gmesh, std::set<int> &matids, std::set<int> &fBCmatIds, int &wrapId,int &interfaceid,int &lagrangeId){
SetPeripheralsids(wrapId,interfaceid,lagrangeId);
DefaultConstructor(gmesh,matids,fBCmatIds);
}

void TPZCreateHybridizedMixedSpaces::SetPeripheralsids( int wrapId,int interfaceId,int lagrangeId){
    fWrapMatId = wrapId;
    fInterfaceMatId = interfaceId;
    fLagrangeMatId = lagrangeId;
}

void TPZCreateHybridizedMixedSpaces::DefaultConstructor(TPZGeoMesh *gmesh, std::set<int> &matids, std::set<int> &bcmatIds){
    fMaterialIds = matids;
    fBCMaterialIds = bcmatIds;
    fGeoMesh = gmesh;
    fDimension = gmesh->Dimension();
}

void TPZCreateHybridizedMixedSpaces::ConditioningGeomesh(){
    int64_t numelem = fGeoMesh->NElements();
    int meshDim = fGeoMesh->Dimension();
    for(int iel =0 ; iel < numelem; iel++){
        TPZGeoEl* gel = fGeoMesh->Element(iel);
        if(!gel || gel->HasSubElement() || gel->Dimension() != meshDim){
            continue;
        }
        int matId = gel->MaterialId();
        if(fMaterialIds.find(matId) == fMaterialIds.end()){
            continue;
        }
        int firstSide = gel->FirstSide(meshDim - 1);
        for(int iside = firstSide; iside < gel->NSides()-1; iside++){
            TPZGeoElSide gelside(gel,iside);
            if(gelside.HasNeighbour(fBCMaterialIds)){
                continue;
            }
            if(!gelside.HasNeighbour(fLagrangeMatId)){
                TPZGeoElBC gelbc(gelside,fLagrangeMatId);
            }
            TPZGeoElBC gelbcInt(gelside,fInterfaceMatId);
            TPZGeoElBC gelbcWrap(gelside,fWrapMatId);
#ifdef PZDEBUG
            {
                auto neighbour = gelside.Neighbour();
                if(neighbour.Element()->MaterialId() !=fWrapMatId){
                    DebugStop();
                }
                auto neighneigh = neighbour.Neighbour();
                if(neighneigh.Element()->MaterialId() != fInterfaceMatId){
                    DebugStop();
                }
                auto neighneighneigh = neighneigh.HasNeighbour(fLagrangeMatId);
                if(!neighneighneigh){
                    DebugStop();
                }
                if(neighneighneigh.Element()->MaterialId() != fLagrangeMatId){
                    DebugStop();
                }
            };
#endif
        }
    }
    {// Printing mesh
#ifdef PZDEBUG
        std::ofstream ofs1("gelCreation.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fGeoMesh,ofs1);
        return;
#endif
    }

}

TPZCompMesh* TPZCreateHybridizedMixedSpaces::CreateDiscHDivMesh(){
    TPZCompMesh *hdivMesh = new TPZCompMesh(fGeoMesh);
    hdivMesh->SetAllCreateFunctionsHDiv();
    hdivMesh->ApproxSpace().CreateDisconnectedElements(true);
    hdivMesh->SetDefaultOrder(fNormalFluxOrder);

    int meshDim = fGeoMesh->Dimension();
    for(auto it : fMaterialIds){
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(it,meshDim);
        hdivMesh->InsertMaterialObject(mat);
    }

    for(auto it : fBCMaterialIds){
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(it,meshDim - 1);
        hdivMesh->InsertMaterialObject(mat);
    }

    {
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(fWrapMatId, meshDim - 1);
        hdivMesh->InsertMaterialObject(mat);
    }

    hdivMesh->AutoBuild(fMaterialIds);
    int numEl = hdivMesh->NElements();
    for(int iel = 0; iel < numEl; iel++){
        auto *cel = hdivMesh->Element(iel);
        auto *intel = dynamic_cast<TPZInterpolatedElement*>(cel);
        auto *gel = cel->Reference();
        if(!cel){
            DebugStop();
        }
        cel->LoadElementReference();
        int firstSide = gel->FirstSide(meshDim - 1);
        for(int iside = firstSide; iside < gel->NSides()-1; iside++){
            TPZGeoElSide gelside(gel,iside);
            auto neigh = gelside.HasNeighbour(fBCMaterialIds);
            if(neigh){
                hdivMesh->ApproxSpace().CreateCompEl(neigh.Element(),*hdivMesh);
                auto *targetcel = hdivMesh->ApproxSpace().CreateCompEl(neigh.Element(),*hdivMesh);
                auto *targetintel = dynamic_cast<TPZInterpolatedElement*>(targetcel);
                targetintel->SetSideOrient(neigh.Side(),1);
                neigh.Element()->ResetReference();
            }
            else{
                neigh = gelside.Neighbour();
                auto *targetcel = hdivMesh->ApproxSpace().CreateCompEl(neigh.Element(),*hdivMesh);
                auto *targetintel = dynamic_cast<TPZInterpolatedElement*>(targetcel);
                targetintel->SetSideOrient(neigh.Side(),1);
                neigh.Element()->ResetReference();
            }
            intel->SetSideOrient(iside,1);
        }
        gel->ResetReference();
    }

    int64_t nc = hdivMesh->NConnects();
    for(int icon = 0; icon < nc; icon++){
        hdivMesh->ConnectVec()[icon].SetLagrangeMultiplier(fConfigHybridizedMixed.lagHDiv);
    }

    int orderplus = fPressureOrder - fNormalFluxOrder;
    TPZCompMeshTools::AdjustFluxPolynomialOrders(hdivMesh, orderplus); //Increases internal flux order by "hdivmais"
    hdivMesh->ExpandSolution();
    return hdivMesh;
}

TPZCompMesh* TPZCreateHybridizedMixedSpaces::CreateL2Mesh(){
    TPZCompMesh *L2Mesh = new TPZCompMesh(fGeoMesh);
    L2Mesh->SetAllCreateFunctionsContinuous();
    L2Mesh->ApproxSpace().CreateDisconnectedElements(true);
    L2Mesh->SetDefaultOrder(fPressureOrder);

    int meshDim = fGeoMesh->Dimension();
    for(auto it : fMaterialIds){
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(it,meshDim);
        L2Mesh->InsertMaterialObject(mat);
    }

    {
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(fLagrangeMatId, meshDim - 1);
        L2Mesh->InsertMaterialObject(mat);
    }

    L2Mesh->AutoBuild(fMaterialIds);

    L2Mesh->SetDimModel(meshDim-1);
    std::set<int> ids = {fLagrangeMatId};
    L2Mesh->SetDefaultOrder(fNormalFluxOrder);
    L2Mesh->AutoBuild(ids);
    L2Mesh->SetDefaultOrder(fPressureOrder);

    int64_t nc = L2Mesh->NConnects();
    for(int icon = 0; icon < nc; icon++){
        L2Mesh->ConnectVec()[icon].SetLagrangeMultiplier(fConfigHybridizedMixed.lagL2);
    }
    L2Mesh->ExpandSolution();
    return L2Mesh;
}

TPZCompMesh* TPZCreateHybridizedMixedSpaces::CreateConstantMesh(const int &lagNum){
    TPZCompMesh *constMesh = new TPZCompMesh(fGeoMesh);
    constMesh->SetAllCreateFunctionsDiscontinuous();
    constMesh->SetDefaultOrder(0);

    int meshDim = fGeoMesh->Dimension();
    for(auto it : fMaterialIds){
        TPZNullMaterial<STATE> *mat = new TPZNullMaterial<STATE>(it,meshDim);
        constMesh->InsertMaterialObject(mat);
    }

    constMesh->AutoBuild(fMaterialIds);

    int64_t nc = constMesh->NConnects();
    for(int icon = 0; icon < nc; icon++){
        constMesh->ConnectVec()[icon].SetLagrangeMultiplier(lagNum);
    }
    return constMesh;
}

void TPZCreateHybridizedMixedSpaces::AddMaterials(TPZMultiphysicsCompMesh *mcmesh, ProblemConfig& config ){

    int meshDim = mcmesh->Dimension();
#ifndef OPTMIZE_RUN_TIME
    TPZMixedDarcyFlowAniso* refmat = 0;
#else
    TPZMixedDarcyFlow* refmat = 0;
#endif
    
    for(auto it : fMaterialIds){
#ifndef OPTMIZE_RUN_TIME
        TPZMixedDarcyFlowAniso *mat = new TPZMixedDarcyFlowAniso(it,meshDim);
#else
        TPZMixedDarcyFlow *mat = new TPZMixedDarcyFlow(it,meshDim);
#endif
        mcmesh->InsertMaterialObject(mat);
        refmat = mat;
#ifndef OPTMIZE_RUN_TIME
        mat->SetPermeability(config.perm);
        mat->SetForcingFunction(fAnalyticSolution->ForceFunc(),5);
        mat->SetExactSol(fAnalyticSolution->ExactSolution(),5);
#else
        mat->SetConstantPermeability(1.);
#endif
    }

    for(auto it : fBCMaterialIds){
        TPZFMatrix<STATE> val1(1,1,1);
        TPZVec<STATE> val2(1,1);
        auto *mat = refmat->CreateBC(refmat,it,0,val1,val2);
        mcmesh->InsertMaterialObject(mat);
#ifndef OPTMIZE_RUN_TIME
        if (fAnalyticSolution.operator*().fExact != TLaplaceExample1::ENone) {
            mat->SetForcingFunctionBC(fAnalyticSolution->ExactSolution(),5);
        }
#endif
    }

    {
        TPZNullMaterialCS<STATE> *mat = new TPZNullMaterialCS<STATE>(fLagrangeMatId, meshDim - 1,1);
        mcmesh->InsertMaterialObject(mat);
    }

    {
        TPZNullMaterialCS<STATE> *mat = new TPZNullMaterialCS<STATE>(fWrapMatId, meshDim - 1,1);
        mcmesh->InsertMaterialObject(mat);
    }
    return;
}

void TPZCreateHybridizedMixedSpaces::AddInterfaceMaterial(TPZMultiphysicsCompMesh *mcmesh){

    int meshDim = mcmesh->Dimension();

    {
        TPZLagrangeMultiplierCS<STATE> *mat = new TPZLagrangeMultiplierCS<STATE>(fInterfaceMatId, meshDim - 1);
        mcmesh->InsertMaterialObject(mat);
    }
}


void TPZCreateHybridizedMixedSpaces::AddInterfaceComputationalElements(TPZMultiphysicsCompMesh *mcmesh){
    int numEl = mcmesh->NElements();

    auto fGeoMesh = mcmesh->Reference();
    fGeoMesh->ResetReference();
    mcmesh->LoadReferences();

    for(int iel = 0; iel < numEl; iel++){
        TPZCompEl *cel = mcmesh->Element(iel);
        if(!cel || !cel->Reference()){
            DebugStop();
        }
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        if(matid != fWrapMatId){
            continue;
        }
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neighLag = gelside.HasNeighbour(fLagrangeMatId);
        if(!neighLag){
            DebugStop();
        }
        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide cneigh = neighLag.Reference();

        TPZGeoElSide ginterface = gelside.Neighbour();
        if(ginterface.Element()->MaterialId() != fInterfaceMatId){
            DebugStop();
        }
        TPZMultiphysicsInterfaceElement *interface = new TPZMultiphysicsInterfaceElement(*mcmesh,ginterface.Element(),celside,cneigh);
    }
}

void TPZCreateHybridizedMixedSpaces::ComputePeriferalMaterialIds(int base)
{
    if(base < 2) base = 2;
    int max_matid = 0;
    int64_t nel = fGeoMesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel) continue;
        max_matid = std::max(max_matid,gel->MaterialId());
    }
    int remain = max_matid % base;
    int matid_base = max_matid-remain + base;

    fLagrangeMatId = matid_base+2*base;
    fWrapMatId = matid_base;
    fInterfaceMatId = matid_base + base;
}

void TPZCreateHybridizedMixedSpaces::Print(std::ostream &ofs){
    std::stringstream ss;
    ss << "\n";
    int numStars = 70;
    for(int istar = 0; istar < numStars ; istar++)
        ss <<"*";
    ss <<"\n";
    ss << __PRETTY_FUNCTION__ << "\nHybridized Mixed approximation space creation set up:\n\n";
    ss << "fSpaceType:\t" <<SpaceTypeName() << "\nNormal Flux Order: " << fNormalFluxOrder;
    ss << "\nPressure Order:\t" << fPressureOrder << "\nfHybridizeBC: " << fHybridizeBC;
    ss << "\n\nmatids:\n\t{";
    for(auto it : fMaterialIds){
        ss << it <<", ";
    }
    ss.seekp(-2,ss.cur); //return stringstream head by two positions, as consequence, this data will be erased.
    ss << "}\n";

    ss << "bcmatids:\n\t{";
    for(auto it : fBCMaterialIds){
        ss << it <<", ";
    }
    ss.seekp(-2,ss.cur); //return stringstream head by two positions, as consequence, this data will be erased.
    ss << "}\n";

    ss << "peripheralMatids:\n";
    ss << "\twrapid = " << fWrapMatId;
    ss << "\n\tinterfaceid = " << fInterfaceMatId;
    ss << "\n\tlagrangeid = " << fLagrangeMatId <<"\n";

    for(int istar = 0; istar < numStars ; istar++)
        ss <<"*";
    ss <<"\n";

    ofs << ss.str();
}

std::string TPZCreateHybridizedMixedSpaces::SpaceTypeName(){
    std::string spaceName;
    switch (fSpaceType) {
        case 0:
            spaceName = "ENone";
            break;
        case 1:
            spaceName = "EHybridizedMixed";
            break;
    }
    return spaceName;
}

TPZMultiphysicsCompMesh* TPZCreateHybridizedMixedSpaces::GenerateMesh(ProblemConfig& config){
    {
        ConditioningGeomesh();
        TPZCompMesh *hdivMesh = CreateDiscHDivMesh();
        TPZCompMesh *L2Mesh = CreateL2Mesh();
        TPZCompMeshTools::SetPressureOrders(hdivMesh, L2Mesh);
        TPZCompMesh *gspace   =   CreateConstantMesh(fConfigHybridizedMixed.lagg);
        TPZCompMesh *avgspace = CreateConstantMesh(fConfigHybridizedMixed.lagavg);

        TPZVec<TPZCompMesh*> meshvec = {hdivMesh,L2Mesh,gspace,avgspace};
        TPZMultiphysicsCompMesh *mcmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
        mcmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

        AddMaterials(mcmesh,config);
        mcmesh->BuildMultiphysicsSpace(meshvec);

        AddInterfaceMaterial(mcmesh);
        AddInterfaceComputationalElements(mcmesh);

        mcmesh->InitializeBlock();

        if(fShouldCondense) {
            GroupAndCondenseElements(mcmesh);
        }
        mcmesh->CleanUpUnconnectedNodes();

        return mcmesh;
    }
}

void TPZCreateHybridizedMixedSpaces::GroupAndCondenseElements(TPZMultiphysicsCompMesh *cmesh)
{
    /// same procedure as hybridize hdiv
    int64_t nel = cmesh->NElements();
    cmesh->LoadReferences();

    TPZCompEl *cel;
    TPZGeoEl *gel;
    TPZGeoElSide *gelside;
    TPZGeoElSide neigh;
    TPZCompEl *celTarget;

    std::set<int64_t> targetMatIds;
    targetMatIds.insert(fWrapMatId);
    targetMatIds.insert(fInterfaceMatId);
    for(auto it : fBCMaterialIds){
        targetMatIds.insert(it);
    }

    std::ofstream ofs1("checkGelmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fGeoMesh,ofs1);

    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (int64_t el = 0; el<nel; el++) {
        cel = cmesh->Element(el);
        if(!cel){
            continue;
        }
        if (cmesh->Element(el)->Dimension() != fDimension) {
            continue;
        }
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh);

        gel = cel->Reference();

        int firstSide;
        if (gel->Dimension() >= 0) {
            firstSide = gel->FirstSide(gel->Dimension() - 1);
        } else {
            DebugStop();
        }

        elgr->AddElement(cel);
        for (int iside = firstSide; iside < gel->NSides() - 1; iside++) {
            gelside = new TPZGeoElSide(gel, iside);
            neigh = *gelside;
            if(gelside->HasNeighbour(fBCMaterialIds)){
                neigh = gelside->Neighbour();
                if(fBCMaterialIds.find(neigh.Element()->MaterialId()) != fBCMaterialIds.end()){
                    celTarget = neigh.Element()->Reference();
                    elgr->AddElement(celTarget);
                }
            } else{
                neigh = gelside->Neighbour();
                std::vector<int> targetMatIds= {fWrapMatId,fInterfaceMatId};
                for (int ineigh = 0; ineigh < 2; ineigh++) {
                    if (neigh.Element()->MaterialId() != targetMatIds[ineigh]) {
                        DebugStop();
                    }
                    celTarget = neigh.Element()->Reference();
                    elgr->AddElement(celTarget);
                    neigh = neigh.Neighbour();
                }
            }
            gelside = NULL;
        }
    }
    cmesh->ComputeNodElCon();
    /*if(fSpaceType == EHybridizedMixed)
    {
        int64_t nconnects = cmesh->NConnects();
        for (int64_t ic = 0; ic<nconnects; ic++) {
            TPZConnect &c = cmesh->ConnectVec()[ic];
            if(c.LagrangeMultiplier() == fConfigHybridizedMixed.lagavg) c.IncrementElConnected();
        }
    }*/
    int neoNel = cmesh->NElements();
    for (int64_t el = nel; el < neoNel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (!elgr) {
            DebugStop();
        }
        TPZCondensedCompEl *cond = new TPZCondensedCompElT<STATE>(elgr);
        cond->SetKeepMatrix(false);
        }
}
