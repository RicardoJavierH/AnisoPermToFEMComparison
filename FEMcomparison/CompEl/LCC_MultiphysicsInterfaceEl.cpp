//
// Created by victor on 03/05/2021.
//

#include "LCC_MultiphysicsInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "pzelmat.h"
#include "pzinterpolationspace.h"
#include "TPZMaterial.h"
#include "pzmultiphysicselement.h"
#include "tpzintpoints.h"

#include "pzmultiphysicscompel.h"
#include "pzgeoel.h"

#include "pzgraphel.h"
#include "pzgraphmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "LCC_LagrangeMultiplier.h"
#include "TPZTimer.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("StiffnessInterface");
#endif

LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(): TPZMultiphysicsInterfaceElement()
{

}

LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, TPZCompElSide left, TPZCompElSide right) :
TPZMultiphysicsInterfaceElement(mesh, ref, left, right)
{

}

LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref):
TPZMultiphysicsInterfaceElement(mesh, ref)
{

}

/** @brief create a copy of the given element */
LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const LCC_TPZMultiphysicsInterfaceElement &copy):
TPZMultiphysicsInterfaceElement(mesh, copy)
{

}

/** @brief create a copy of the given element using index mapping */
LCC_TPZMultiphysicsInterfaceElement::LCC_TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const LCC_TPZMultiphysicsInterfaceElement &copy, std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t> & gl2lcElMap):
TPZMultiphysicsInterfaceElement(mesh, copy, gl2lcConMap,gl2lcElMap)
{

}
template<class TVar>
void LCC_TPZMultiphysicsInterfaceElement::CalcStiffT(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef) {
#ifdef FEMCOMPARISON_TIMER
    //extern double interfaceTime;
    //TPZTimer timer;
    //timer.start();
#endif
    
    if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
    TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (fLeftElSide.Element());
    TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
    TPZGeoEl *leftgel = leftel->Reference();
    TPZGeoEl *rightgel = rightel->Reference();

    InitializeElementMatrix(ek, ef);

    int NNodes = leftgel->NNodes();
    TPZManVector<int> leftNodeIndices;
    leftNodeIndices.Resize(NNodes,0);

    for (int iNode = 0; iNode < leftgel->NCornerNodes(); iNode++) {
        leftNodeIndices[iNode] = leftgel->NodeIndex(iNode);
    }

    std::pair<int, int> smallerNode;
    smallerNode.first = 0;
    smallerNode.second = leftNodeIndices[0];
    for(int iNode = 1; iNode < leftgel->NCornerNodes(); iNode++){
        if(leftNodeIndices[iNode] < smallerNode.second){
            smallerNode.first = iNode;
            smallerNode.second = leftNodeIndices[iNode];
        }
    }

    if(NNodes == 2){
        if(this->Dimension() != 1){
            DebugStop();
        }
        if(leftNodeIndices[0] < leftNodeIndices[1]){ // Matrix A/MinusA
            ChoosingOptimizedComputation( ek, ef,0);
        }
        else{                                        // Matrix B/MinusB
            ChoosingOptimizedComputation( ek, ef,1);
        }
    }
    else if(NNodes == 4){
        if(this->Dimension() != 2){
            DebugStop();
        }
        switch(smallerNode.first) {
            case 0:
                if(leftNodeIndices[1] < leftNodeIndices[3]){ // Case A;
                    ChoosingOptimizedComputation( ek, ef,0);
                }
                else{
                    ChoosingOptimizedComputation( ek, ef,1);
                }
                break;
            case 1:
                if(leftNodeIndices[2] < leftNodeIndices[0]){ // Case A;
                    ChoosingOptimizedComputation( ek, ef,2);
                }
                else{
                    ChoosingOptimizedComputation( ek, ef,3);
                }
                break;
            case 2:
                if(leftNodeIndices[3] < leftNodeIndices[1]){ // Case A;
                    ChoosingOptimizedComputation( ek, ef,4);
                }
                else{
                    ChoosingOptimizedComputation( ek, ef,5);
                }
                break;
            case 3:
                if(leftNodeIndices[0] < leftNodeIndices[2]){ // Case A;
                    ChoosingOptimizedComputation( ek, ef,6);
                }
                else{
                    ChoosingOptimizedComputation( ek, ef,7);
                }
                break;
            default:
                DebugStop();
        }
    }
    else if(NNodes == 3){
        if(this->Dimension() != 2){
            DebugStop();
        }
        switch(smallerNode.first) {
            case 0:
                if(leftNodeIndices[1] < leftNodeIndices[2]){ // Case A;
                    ChoosingOptimizedComputation( ek, ef,8);
                }
                else{
                    ChoosingOptimizedComputation( ek, ef,9);
                }
                break;
            case 1:
                if(leftNodeIndices[2] < leftNodeIndices[0]){ // Case A;
                    ChoosingOptimizedComputation( ek, ef,10);
                }
                else{
                    ChoosingOptimizedComputation( ek, ef,11);
                }
                break;
            case 2:
                if(leftNodeIndices[0] < leftNodeIndices[1]){ // Case A;
                    ChoosingOptimizedComputation( ek, ef,12);
                }
                else{
                    ChoosingOptimizedComputation( ek, ef,13);
                }
                break;
            default:
                DebugStop();
        }
    }
    else{
        ComputingCalcStiff(ek, ef);
    }
#ifdef FEMCOMPARISON_TIMER
    //timer.stop();
    //interfaceTime+=timer.seconds();
#endif
}
template<class TVar>
void LCC_TPZMultiphysicsInterfaceElement::ChoosingOptimizedComputation(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef, int matrixIndex){

    TPZMaterial *material = this->Material();
    if (!material) {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        ek.Reset();
        ef.Reset();
        return;
    }
    LCC_LagrangeMultiplier *lagMat = dynamic_cast<LCC_LagrangeMultiplier *>(material);
    if (!lagMat) {
        DebugStop();
    }

    TPZFMatrix<STATE> M;
    lagMat->GetStiffnessMatrix(M,matrixIndex);

    if(M.Cols() == 0){
        //FIXME: to review
        M.Resize(ek.fMat.Rows(), ek.fMat.Cols());
        ComputingCalcStiff(ek,ef);
        lagMat->FillStiffnessMatrix(ek.fMat,matrixIndex);
    }
    else{
#ifdef FEMCOMPARISON_DEBUG
        //DebugStop();
        ComputingCalcStiff(ek,ef);
        if(ek.fMat.Cols() != M.Cols() || ek.fMat.Rows() != M.Rows()){
            DebugStop();
        }
        STATE relDiff;
        for(int iRow = 0; iRow < M.Rows(); iRow++){
            for(int iCol = 0 ; iCol < M.Cols(); iCol++){
                relDiff = (M(iRow,iCol) - ek.fMat(iRow,iCol))/M(iRow,iCol);
                if(abs(M(iRow,iCol)) > 1.0e-16){
                    if(abs(relDiff) > 0.000001){
                        std::cout << "M: " << M(iRow,iCol) <<"; ek: " << ek.fMat(iRow,iCol) <<"\n";
                        DebugStop();
                    }
                }
            }
        }
#endif
        ek.fMat = M;
    }
}
template<class TVar>
void LCC_TPZMultiphysicsInterfaceElement::ComputingCalcStiff(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef) {

    TPZMaterial *material = this->Material();
    auto *matInterface =
        dynamic_cast<TPZMatInterfaceCombinedSpaces<TVar>*>(material);
     if(!material || !matInterface){
         PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
         ek.Reset();
         ef.Reset();
         return;
     }
    
    if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
    TPZMultiphysicsElement *leftel = dynamic_cast<TPZMultiphysicsElement *> (fLeftElSide.Element());
    TPZMultiphysicsElement *rightel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
    TPZGeoEl *leftgel = leftel->Reference();
    TPZGeoEl *rightgel = rightel->Reference();

#ifdef PZDEBUG
    if (!leftel || !rightel) {
        DebugStop();
    }
#endif
    
    std::map<int,TPZMaterialDataT<TVar>> datavecleft;
    std::map<int,TPZMaterialDataT<TVar>> datavecright;
    TPZMaterialDataT<TVar> data;
    
    InitMaterialData(data, datavecleft, datavecright);

    TPZManVector<TPZTransform<REAL>,6> leftcomptr, rightcomptr;
    leftel->AffineTransform(leftcomptr);
    rightel->AffineTransform(rightcomptr);

    for(auto &it : datavecleft){
        it.second.p = 0;
        int id = it.first;
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(leftel->Element(id));
        if (msp)
        {
            datavecleft[id].p =msp->MaxOrder();
        }
    }
//    for(auto &it : datavecright){
//        it.second.p = 0;
//        int id = it.first;
//        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(rightel->Element(id));
//        if (msp)
//        {
//            datavecright[id].p = msp->MaxOrder();
//        }
//    }

    TPZManVector<int> intleftorder;
    leftel->PolynomialOrder(intleftorder);
    TPZManVector<int> intrightorder;
    rightel->PolynomialOrder(intrightorder);
    int integrationorder = matInterface->GetIntegrationOrder(intleftorder, intrightorder);
    TPZGeoEl *gel = Reference();
    int dimension = gel->Dimension();
    int thisside = gel->NSides()-1;
    TPZFNMatrix<9,REAL> jac(dimension,dimension),axes(dimension,3), jacInv(dimension,dimension);

    TPZAutoPointer<TPZIntPoints> intrule = gel->CreateSideIntegrationRule(thisside, integrationorder);
    TPZManVector<REAL,3> Point(dimension), leftPoint(leftel->Dimension()), rightPoint(rightel->Dimension());
    TPZGeoElSide neighleft(fLeftElSide.Reference()), neighright(fRightElSide.Reference());
    TPZTransform<> trleft(dimension),trright(dimension);
    TPZGeoElSide gelside(this->Reference(),thisside);
    // compute the transformation between neighbours
    gelside.SideTransform3(neighleft, trleft);
    gelside.SideTransform3(neighright, trright);

    TPZTransform<> leftloctr = leftgel->SideToSideTransform(neighleft.Side(), leftgel->NSides()-1);
    TPZTransform<> rightloctr = rightgel->SideToSideTransform(neighright.Side(), rightgel->NSides()-1);
    // transform from the element to the interior of the neighbours
    trleft = leftloctr.Multiply(trleft);
    trright = rightloctr.Multiply(trright);


    int nintpoints = intrule->NPoints();
    for (int ip =0; ip<nintpoints; ip++) {
        REAL weight;
        data.intLocPtIndex = ip;
        intrule->Point(ip, Point, weight);
        ComputeRequiredData(data, Point);
        weight *= fabs(data.detjac);
        trleft.Apply(Point, leftPoint);
        leftel->ComputeRequiredData(leftPoint, leftcomptr, datavecleft);
        trright.Apply(Point, rightPoint);
        rightel->ComputeRequiredData(rightPoint, rightcomptr, datavecright);

        data.x = datavecleft.begin()->second.x;
        matInterface->ContributeInterface(data, datavecleft, datavecright, weight, ek.fMat, ef.fMat);
    }

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "\nH1-bound element:\n    Element id:" << leftgel->Id() << "\n    Geometric nodes: (";
        for (int i = 0; i < leftgel->NCornerNodes(); i++) {
            sout << leftgel->NodeIndex(i) << " ";
        }
        sout << ")\n    CONNECT INFO. BY GEOMETRIC ELEMENT:\n";
        TPZManVector<REAL, 3> xicenter(gel->Dimension(), 0.);
        TPZManVector<REAL> xcenter(3, 0.);
        for (int i = 0; i < leftgel->NSides(); i++) {
            leftgel->CenterPoint(i, xicenter);
            leftgel->X(xicenter, xcenter);
            sout << "        Cord = [" << xcenter;
            int iCon = i;
            if (leftel->NConnects()) {
                TPZConnect &con = leftel->Connect(iCon);
                sout << "] Side = " << i;
                if (1) {
                    sout << " Seqnumber = " << con.SequenceNumber();
                }
                sout << " Order = " << (int) con.Order() << " NState = " << (int) con.NState()
                     << " NShape " << con.NShape();

                if (1) {
                    sout << " IsLagrangeMult = " << (int) con.LagrangeMultiplier();
                }
                sout << " IsCondensed: " << (int) con.IsCondensed();

                sout << '\n';
            }
        }

        TPZFMatrix<STATE> Nodes;
        TPZFMatrix<STATE> X;
        int nNodes = leftgel->NNodes();
        int dim = leftgel->Dimension();
        int nPoints;
        switch(nNodes){
            case 2:
                nPoints = 2;
                break;
            case 3:
                nPoints = 3;
                break;
            case 4:
                nPoints = 10;
                break;
            default:
                DebugStop();
        }

        Nodes.Resize(nPoints,dim);
        X.Resize(nPoints,3);

        if(nNodes == 2){
            Nodes(0,0) = -1.;
            Nodes(0,0) = +1.;
        }

        if(nNodes == 3){
            Nodes(0,0) = 0.; Nodes(0,1) = 0.;
            Nodes(1,0) = 1.; Nodes(1,1) = 0.;
            Nodes(2,0) = 0.; Nodes(2,1) = 1.;
        }

        if(leftgel->NNodes() == 4) {
            Nodes(0,0) = -1.;  Nodes(0,1) = -1.;
            Nodes(1,0) = 1.;   Nodes(1,1) = -1.;
            Nodes(2,0) = 1.;   Nodes(2,1) = 1.;
            Nodes(3,0) = -1.;  Nodes(3,1) = 1.;
            Nodes(4,0) = 0.;   Nodes(4,1) = -1.;
            Nodes(5,0) = 1.;   Nodes(5,1) = 0.;
            Nodes(6,0) = 0.;   Nodes(6,1) = 1.;
            Nodes(7,0) = -1.;  Nodes(7,1) = 0.;
            Nodes(8,0) = 0.;   Nodes(8,1) = 0.;
            Nodes(9,0) = -0.5; Nodes(9,1) = -0.5;
        }


        TPZVec<STATE> xi(dim), x(3);
        for(int iNode = 0; iNode < nPoints; iNode++){
            for(int iCor = 0 ; iCor <dim ; iCor++){
                xi[0] = Nodes(iNode,iCor);
            }
            leftgel->X(xi,x);
            for(int iCor = 0 ; iCor <3 ; iCor++){
                X(iNode,iCor) = x[iCor];
            }
        }

        for(int i =0 ; i < nPoints; i++){
            sout << "        Node coord:  (" << X(i,0) << ", " << X(i,1) << ", " << X(i,2) << "): \n";
            leftel->AffineTransform(leftcomptr);
            for(int iCor = 0 ; iCor <dim ; iCor++){
                xi[0] = Nodes(i,iCor);
            }
            leftel->ComputeRequiredData(xi, leftcomptr, datavecleft);
            for(auto &it : datavecleft){
                sout << "            Space " << it.first<< ": " ;
                TPZFMatrix<REAL>  &phi = it.second.phi;
                for(int iPhi = 0; iPhi < phi.Rows(); iPhi++)
                    sout << phi(iPhi,0) << ", ";
                sout << "\n";
            }
        }

        datavecleft[0].axes.Print(sout);

        sout << "\nH(div)-bound element:\n    Element id:" << rightgel->Id() << "\n    Geometric nodes: (";
        for (int i = 0; i < rightgel->NCornerNodes(); i++) {
            sout << rightgel->NodeIndex(i) << " ";
        }
        sout << ")\n    CONNECT INFO. BY GEOMETRIC ELEMENT:\n";
        int nCon = rightel->NConnects();
        int nSides = gel->NSides();
        int firstSideToHaveConnect = 0;
        if (nSides != nCon) {
            firstSideToHaveConnect = gel->NSides() - 1;
        }
        for (int i = firstSideToHaveConnect; i < nSides; i++) {
            rightgel->CenterPoint(i, xicenter);
            rightgel->X(xicenter, xcenter);
            sout << "        Cord = [" << xcenter;

            int iCon = i - firstSideToHaveConnect;
            if (nCon) {
                TPZConnect &con = rightel->Connect(iCon);
                sout << "] Side = " << i;
                if (1) {
                    sout << " Seqnumber = " << con.SequenceNumber();
                }
                sout << " Order = " << (int) con.Order() << " NState = " << (int) con.NState()
                     << " NShape " << con.NShape();

                if (1) {
                    sout << " IsLagrangeMult = " << (int) con.LagrangeMultiplier();
                }
                sout << " IsCondensed: " << (int) con.IsCondensed();

                sout << '\n';
            }
        }

        sout << "    SHAPE FUNCTION VALUE PER NODE:\n";

        for(int i =0 ; i < nPoints; i++){
            sout << "        Node coord:  (" << X(i,0) << ", " << X(i,1) << ", " << X(i,2) << "): \n";
            for(int iCor = 0 ; iCor <dim ; iCor++){
                xi[0] = Nodes(i,iCor);
            }
            rightel->AffineTransform(rightcomptr);
            rightel->ComputeRequiredData(xi, rightcomptr, datavecright);
            for(auto &it : datavecright){
                sout << "            Space " << it.first<< ": " ;
                TPZFMatrix<REAL>  &phi = it.second.phi;
                for(int iPhi = 0; iPhi < phi.Rows(); iPhi++)
                    sout << phi(iPhi,0) << ", ";
                sout << "\n";
            }
        }
        datavecright[0].axes.Print(sout);

        LCC_LagrangeMultiplier *lagMat = dynamic_cast<LCC_LagrangeMultiplier *>(material);
        if (lagMat) {
            sout << "\nfMultiplier = " << lagMat->Multiplier() << "\n";
        }
        ek.fMat.Print(sout);
        sout
                << "----------------------------------------------------------------------------------------------------------------------------------------\n";
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif

}//CalcStiff
