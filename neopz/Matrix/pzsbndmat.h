/**
 * @file
 * @brief Contains TPZSBMatrix class which implements symmetric band matrices(hermitian, for the complex case. assumed to be
 upper triangular). Some functionalities depend on LAPACK.

 The functionalities that depend on LAPACK will result in runtime error if LAPACK is not linked to NeoPZ. Search for LAPACK in this header to 
 know which functions are affected by this dependency.
 LAPACK can be linked by setting USING_LAPACK=ON or USING_MKL=ON on CMake
when configuring the library.
 */

#ifndef TSBNDMATH
#define TSBNDMATH

#include "pzmatrix.h"
#include "pzfmatrix.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif
template<class TVar>
class TPZLapackEigenSolver;
/**
 * @brief Implements symmetric band matrices. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSBMatrix : public TPZMatrix<TVar>
{
    friend class TPZLapackEigenSolver<TVar>;
public:
    TPZSBMatrix() : TPZRegisterClassId(&TPZSBMatrix::ClassId),
    TPZMatrix<TVar>() , fStorage() { fBand = 0; this->fSymProp = SymProp::Herm;}
    TPZSBMatrix(const int64_t dim,const int64_t band );
    TPZSBMatrix(const TPZSBMatrix<TVar> &A ) = default;
    TPZSBMatrix(TPZSBMatrix<TVar> &&A ) = default;
    inline TPZSBMatrix<TVar>*NewMatrix() const override {return new TPZSBMatrix<TVar>{};}
    CLONEDEF(TPZSBMatrix)
    TPZSBMatrix &operator= (TPZSBMatrix<TVar> &&A ) = default;
    TPZSBMatrix &operator= (const TPZSBMatrix<TVar> &A ) = default;
    ~TPZSBMatrix() { Clear(); }
    
    int    PutVal(const int64_t row,const int64_t col,const TVar& element ) override;
    const TVar GetVal(const int64_t row,const int64_t col ) const override;
    
    TVar &operator()(int64_t row, int64_t col);
    
    /** @brief Sets symmetry property of current matrix (only hermitian/symmetric allowed)*/
    void SetSymmetry (SymProp sp) override;

    friend class TPZSBMatrix<float>;
    friend class TPZSBMatrix<double>;
    
    /// copy the values from a matrix with a different precision
    template<class TVar2>
    void CopyFromDiffPrecision(TPZSBMatrix<TVar2> &orig)
    {
        TPZMatrix<TVar>::CopyFromDiffPrecision(orig);
        fStorage.resize(orig.fStorage.size());
        int64_t nel = fStorage.size();
        for (int64_t el=0; el<nel; el++) {
            fStorage[el] = orig.fStorage[el];
        }
    }

    /** @brief Creates a copy from another TPZSBMatrix*/
    void CopyFrom(const TPZMatrix<TVar> *  mat) override
    {                                                           
        auto *from = dynamic_cast<const TPZSBMatrix<TVar> *>(mat);                
        if (from) {                                               
            *this = *from;                                          
        }                                                         
        else                                                      
            {                                                       
                PZError<<__PRETTY_FUNCTION__;                         
                PZError<<"\nERROR: Called with incompatible type\n."; 
                PZError<<"Aborting...\n";                             
                DebugStop();                                          
            }                                                       
    }
    
    /** @brief Computes z = beta * y + alpha * opt(this)*x */
    /** z and x cannot overlap in memory */
    void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;
    
    void Print(const char *name = NULL, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const override;
    //friend std::ostream & operator<< <>(std::ostream& out,const TPZSBMatrix<TVar>  &A); Leonardo removendo o '<>' antes do (std...
    template<class TT>friend std::ostream & operator<< (std::ostream& out,const TPZSBMatrix<TT>  &A);
    
    /** Fill the matrix with random values (non singular matrix) */
    void AutoFill(int64_t nrow, int64_t ncol, SymProp symmetric) override;
    
    
    /// Operadores com matrizes SKY LINE.
    // @{
    TPZSBMatrix operator+  (const TPZSBMatrix<TVar> &A ) const;
    TPZSBMatrix operator-  (const TPZSBMatrix<TVar> &A ) const;
    TPZSBMatrix &operator+=(const TPZSBMatrix<TVar> &A );
    TPZSBMatrix &operator-=(const TPZSBMatrix<TVar> &A );
    // @}
    TPZSBMatrix<TVar> operator*  (const TVar v ) const;
    TPZSBMatrix<TVar> &operator*=(const TVar v ) override;
    
    TPZSBMatrix<TVar> operator-() const { return operator*(-1.0); }
    
    /// Redimension the matrix keeping original elements.
    int Resize(const int64_t newDim ,const int64_t) override;
    
    /// Redimension the matrix and zeroes its elements
    int Redim(const int64_t newDim) {return Redim(newDim,newDim);}
    int Redim(const int64_t newRows ,const int64_t newCols) override;
    
    /// Zeroes the elements of the matrix
    int Zero() override;
    
    int64_t GetBand() const { return fBand; }
    int   SetBand(const int64_t newBand );
    
    /// To solve linear systems
    // @{
    
    /** @brief decompose the system of equations acording to the decomposition
      * scheme */
     virtual int Decompose(const DecomposeType dt) override {
       switch (dt) {
       case ELDLt:
         return Decompose_LDLt();
         break;
       case ECholesky:
         return Decompose_Cholesky();
         break;
       default:
         DebugStop();
         break;
       }
       return -1;
     }

    
    int SolveDirect( TPZFMatrix<TVar> &B , const DecomposeType dt) override {
        
        switch ( dt ) {
            case ECholesky:
                return( Solve_Cholesky( &B )  );
            case ELDLt:
                return( Solve_LDLt( &B )  );
            default:
                this->Error( "Solve  < Unknown decomposition type >" );
                break;
        }
        return ( 0 );
    }
    int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const override
    {
        if(this->fDecomposed != dt) DebugStop();
        switch ( dt ) {
            case ECholesky:
                return ( Subst_Forward(&F) && Subst_Backward(&F) );
            case ELDLt:
                return( Subst_LForward( &F ) && Subst_Diag( &F ) && Subst_LBackward( &F ) );
            default:
                this->Error( "Solve  < Unhandled decomposition type >" );
                break;
        }
        return ( 0 );
    }

    /**********************/
    /*** Solve Cholesky ***/
    //
    //  Se nao conseguir resolver por Cholesky retorna 0 e a matriz
    //   sera' modificada (seu valor perdera' o sentido).
    //
    int Solve_Cholesky( TPZFMatrix<TVar>* B )
    {
        return(
               ( !Decompose_Cholesky() )?  0 :( Subst_Forward( B ) && Subst_Backward( B ) )
               );
    }

    /******************/
    /*** Solve LDLt ***/

    int Solve_LDLt( TPZFMatrix<TVar>* B ) {
        
        return(
               ( !Decompose_LDLt() )? 0 :
               ( Subst_LForward( B ) && Subst_Diag( B ) && Subst_LBackward( B ) )
               );
    }

    


    //If LAPACK is available, it will use its implementation.
    int Decompose_Cholesky();  // Faz A = GGt.
    
    int Subst_Forward( TPZFMatrix<TVar>*B ) const;
    int Subst_Backward ( TPZFMatrix<TVar> *b ) const;

    int Decompose_LDLt();
    int Subst_LForward( TPZFMatrix<TVar> *B ) const;
    int Subst_LBackward( TPZFMatrix<TVar> *B ) const;
    int Subst_Diag( TPZFMatrix<TVar> *B ) const;
//    int Subst_Forward( TPZFMatrix<TVar>*B ) const;
//    int Subst_Backward( TPZFMatrix<TVar> *B ) const;
    
    // @}
    
    /*** @name Solve eigenvalues. Depends on LAPACK ***/
    /** @{ */
    
    /// Computes the eigenvalues and eigenvectors of the symmetric matrix
    // on exit the matrix contains the eigenvectors
    /** @brief Solves the Ax=w*x eigenvalue problem and calculates the eigenvectors. Depends on LAPACK.
     * @param w Stores the eigenvalues
     * @param Stores the correspondent eigenvectors
     * @return Returns info param from LAPACK(0 if executed correctly)
     */
    int SolveEigenProblem(TPZVec < CTVar > &w, TPZFMatrix < CTVar > &eigenVectors);
    /** @brief Solves the Ax=w*x eigenvalue problem and does NOT calculates the eigenvectors. Depends on LAPACK.
     * @param w Stores the eigenvalues
     * @return Returns info param from LAPACK(0 if executed correctly)
     */
    int SolveEigenProblem(TPZVec < CTVar > &w);
    /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and calculates the eigenvectors. Depends on LAPACK.
     * @param w Stores the eigenvalues
     * @param Stores the correspondent eigenvectors
     * @return Returns info param from LAPACK(0 if executed correctly)
     */
    int SolveGeneralisedEigenProblem(TPZSBMatrix< TVar > &B , TPZVec < CTVar > &w, TPZFMatrix < CTVar > &eigenVectors);
    /** @brief Solves the generalised Ax=w*B*x eigenvalue problem and does NOT calculates the eigenvectors. Depends on LAPACK.
     * @param w Stores the eigenvalues
     * @return Returns info param from LAPACK(0 if executed correctly)
     */
    int SolveGeneralisedEigenProblem(TPZSBMatrix< TVar > &B , TPZVec < CTVar > &w);
    
    /** @} */

int ClassId() const override;

protected:
    /** @brief Checks compatibility of matrices before Add/Subtract operations*/
    inline void CheckTypeCompatibility(const TPZMatrix<TVar>*A,
                                       const TPZMatrix<TVar>*B)const override
    {
        auto aPtr = dynamic_cast<const TPZSBMatrix<TVar>*>(A);
        auto bPtr = dynamic_cast<const TPZSBMatrix<TVar>*>(B);
        if(!aPtr || !bPtr || aPtr->fBand!=bPtr->fBand){
            PZError<<__PRETTY_FUNCTION__;
            PZError<<"\nERROR: incompatible matrices\n.Aborting...\n";
            DebugStop();
        }
    }
    inline TVar *&Elem() override
    {
        return fStorage.begin();
    }
    inline const TVar *Elem() const override
    {
        return fStorage.begin();
    }

    
    inline int64_t  Size() const override
    {
        return( this->Dim() * (fBand + 1) );
    }
private:
//    int  PutZero();
    //static int  Error(const char *msg1,const char* msg2="" ) ;
    int  Clear() override;
    void Copy (const TPZSBMatrix<TVar> & );
    
    int64_t Index(int64_t i, int64_t j) const
    {
#ifdef PZDEBUG
        if (i>j) {
            DebugStop();
        }
#endif
        return fBand+i-j+(fBand+1)*j;
    }
    
    TPZVec<TVar> &Storage() {
        return fStorage;
    }
private:
    
    TPZVec<TVar> fStorage;
    int64_t  fBand;
};

#endif
