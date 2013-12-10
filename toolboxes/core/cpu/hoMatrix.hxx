namespace Gadgetron
{

template <typename T> 
hoMatrix<T>::hoMatrix() : BaseClass()
{
}

template <typename T> 
hoMatrix<T>::hoMatrix(unsigned int rows, unsigned int cols) : BaseClass(cols, rows)
{
    this->fill(T(0));
}

template <typename T> 
hoMatrix<T>::hoMatrix(unsigned int rows, unsigned int cols, T* data, bool delete_data_on_destruct)
{
    std::vector<unsigned int> dim(2);
    dim[0] = sx;
    dim[1] = sy;
    this->create(dimensions,data,delete_data_on_destruct);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
hoMatrix<T>::~hoMatrix()
{

}

template <typename T> 
bool hoMatrix<T>::createMatrix(unsigned int rows, unsigned int cols)
{
    return this->createArray(cols, rows);
}

template <typename T> 
inline T& hoMatrix<T>::operator()(unsigned int r , unsigned int c)
{
    GADGET_DEBUG_CHECK_THROW(c<(*dimensions_)[0] && r<(*dimensions_)[1]);
    return accesser_[r][c];
}

template <typename T> 
inline const T& hoMatrix<T>::operator()(unsigned int r , unsigned int c) const
{
    GADGET_DEBUG_CHECK_THROW(c<(*dimensions_)[0] && r<(*dimensions_)[1]);
    return accesser_[r][c];
}

template <typename T> 
inline unsigned int hoMatrix<T>::rows() const
{
    return (*dimensions_)[1];
}

template <typename T> 
inline unsigned int hoMatrix<T>::cols() const
{
    return (*dimensions_)[0];
}

template <typename T> 
bool hoMatrix<T>::upperTri(const T& v)
{
    try
    {
        unsigned int r, c;
        for (r=0; r<(*dimensions_)[1]; r++)
        {
            for (c=0; c<(*dimensions_)[0]; c++)
            {
                if ( c > r )
                {
                    (*this)(r, c) = v;
                }
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in hoMatrix<T>::upperTri(const T& v) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool hoMatrix<T>::lowerTri(const T& v)
{
    try
    {
        unsigned int r, c;
        for (r=0; r<(*dimensions_)[1]; r++)
        {
            for (c=0; c<(*dimensions_)[0]; c++)
            {
                if ( r > c )
                {
                    (*this)(r, c) = v;
                }
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in hoMatrix<T>::lowerTri(const T& v) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool hoMatrix<T>::operator == (const Self& m) const
{
    GADGET_CHECK_RETURN_FALSE(this->dimensions_equal(&m));
    for ( size_t i=0; i<elements_; i++ )
    { 
        if (std::abs(data_[i]-m.data_[i])>DBL_EPSILON) 
        {
            return false;
        }
    }
    return true;
}

template <typename T> 
bool hoMatrix<T>::operator != (const Self& m) const
{
    return !(*this==m);
}

template <typename T> 
void hoMatrix<T>::print(std::ostream& os) const
{
    using namespace std;
    os.unsetf(std::ios::scientific);

    os << "hoMatrix : " << (*dimensions_)[1] << " " << (*dimensions_)[0] << " : " << std::string(typeid(T).name()) << endl;
    unsigned int r, c;
    for (r=0; r<(*dimensions_)[1]; r++) 
    {
        os << "r " << r << ":\t";
        for (c=0; c<(*dimensions_)[0]; c++)
        {
            os << setprecision(16) << (*this)(r,c) << "\t";
        }
        os << endl; 
    }
}

template <typename T> 
bool copyL2U(hoMatrix<T>& A)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

        unsigned int R = A.rows();
        unsigned int C = A.cols();

        unsigned int row, col;
        for(row=0; row<R; row++) 
        {
            for(col=0; col<row; col++ )
            {
                A(col, row) = A(row, col);
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in copyL2U(hoMatrix<T>& A) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool copyL2U(hoMatrix<T>& A, bool conj)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

        unsigned int R = A.rows();
        unsigned int row, col;

        if ( conj )
        {
            for(row=0; row<R; row++) 
            {
                for(col=0; col<row; col++ )
                {
                    A(col, row) = std::conj(A(row, col));
                }
            }
        }
        else
        {
            for(row=0; row<R; row++) 
            {
                for(col=0; col<row; col++ )
                {
                    A(col, row) = A(row, col);
                }
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in copyL2U(hoMatrix<T>& A, bool conj) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool copyU2L(hoMatrix<T>& A)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

        unsigned int R = A.rows();
        unsigned int C = A.cols();

        unsigned int row, col;
        for(row=0; row<R; row++) 
        {
            for(col=row+1; col<C; col++ )
            {
                A(col, row) = A(row, col);
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in copyU2L(hoMatrix<T>& A) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool copyU2L(hoMatrix<T>& A, bool conj)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

        unsigned int R = A.rows();
        unsigned int C = A.cols();

        unsigned int row, col;

        if ( conj )
        {
            for(row=0; row<R; row++) 
            {
                for(col=row+1; col<C; col++ )
                {
                    A(col, row) = std::conj(A(row, col));
                }
            }
        }
        else
        {
            for(row=0; row<R; row++) 
            {
                for(col=row+1; col<C; col++ )
                {
                    A(col, row) = A(row, col);
                }
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in copyU2L(hoMatrix<T>& A, bool conj) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool trans(const hoMatrix<T>& A, hoMatrix<T>& AT)
{
    try
    {
        if ( A.get_number_of_elements() == 0 ) return true;

        if ( !AT.dimensions_equal(&A) )
        {
            AT.createMatrix(A.rows(), A.cols());
        }

        int r, c;
        #pragma omp parallel for default(none) private(r, c) shared(A, AT)
        for ( c=0; c<(int)A.cols(); c++ )
        {
            for ( r=0; r<(int)A.rows(); r++ )
            {
                AT(c,r) = A(r,c);
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in trans(const hoMatrix<T>& A, hoMatrix<T>& AT) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool conjugatetrans(const hoMatrix<T>& A, hoMatrix<T>& AH)
{
    try
    {
        if ( A.get_number_of_elements() == 0 ) return true;

        if ( !AH.dimensions_equal(&A) )
        {
            AH.createMatrix(A.rows(), A.cols());
        }

        int r, c;
        #pragma omp parallel for default(none) private(r, c) shared(A, AH)
        for ( c=0; c<(int)A.cols(); c++ )
        {
            for ( r=0; r<(int)A.rows(); r++ )
            {
                AH(c,r) = std::conj(A(r,c));
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in conjugatetrans(const hoMatrix<T>& A, hoMatrix<T>& AH) ... ");
        return false;
    }
    return true;
}

// following matrix computation calls MKL functions
#ifdef USE_MKL

template<typename T> 
bool GeneralMatrixProduct_gemm(hoMatrix<T>& C, 
                            const hoMatrix<T>& A, bool transA, 
                            const hoMatrix<T>& B, bool transB)
{
    try
    {
        CBLAS_TRANSPOSE TA, TB;

        MKL_INT lda = A.cols();
        MKL_INT ldb = B.cols();
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        MKL_INT M = A.rows();
        MKL_INT K = A.cols();
        if ( transA )
        { 
            M = A.cols();
            K = A.rows();
        }

        MKL_INT N = B.cols();
        MKL_INT K2 = B.rows();
        if ( transB )
        { 
            N = B.rows();
            K2 = B.cols();
        }

        GADGET_CHECK_RETURN_FALSE(K==K2);
        if ( (C.rows()!=M) || (C.cols()!=N) )
        {
            GADGET_CHECK_RETURN_FALSE(C.createMatrix(M, N));
        }

        T* pC = C.begin();
        MKL_INT ldc = C.cols();

        T alpha(1), beta(0);

        if ( typeid(T)==typeid(float) )
        {
            if ( transA )
            {
                TA = CblasTrans;
            }
            else
            {
                TA = CblasNoTrans;
            }

            if ( transB )
            {
                TB = CblasTrans;
            }
            else
            {
                TB = CblasNoTrans;
            }

            if ( &A != &C )
            {
                cblas_sgemm(CblasRowMajor, TA, TB, M, N, K, 1, reinterpret_cast<const float*>(pA), lda, reinterpret_cast<const float*>(pB), ldb, 0, reinterpret_cast<float*>(pC), ldc);
            }
            else
            {
                hoMatrix<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_sgemm(CblasRowMajor, TransA, TransB, M, N, K, 1, reinterpret_cast<const float*>(pATmp), lda, reinterpret_cast<const float*>(pB), ldb, 0, reinterpret_cast<float*>(pC), ldc);
            }
        }
        else if ( typeid(T)==typeid(double) )
        {
            if ( transA )
            {
                TA = CblasTrans;
            }
            else
            {
                TA = CblasNoTrans;
            }

            if ( transB )
            {
                TB = CblasTrans;
            }
            else
            {
                TB = CblasNoTrans;
            }

            if ( &A != &C )
            {
                cblas_dgemm(CblasRowMajor, TA, TB, M, N, K, 1, reinterpret_cast<const double*>(pA), lda, reinterpret_cast<const double*>(pB), ldb, 0, reinterpret_cast<double*>(pC), ldc);
            }
            else
            {
                hoMatrix<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_dgemm(CblasRowMajor, TransA, TransB, M, N, K, 1, reinterpret_cast<const double*>(pATmp), lda, reinterpret_cast<const double*>(pB), ldb, 0, reinterpret_cast<double*>(pC), ldc);
            }
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            if ( transA )
            {
                TA = CblasConjTrans;
            }
            else
            {
                TA = CblasNoTrans;
            }

            if ( transB )
            {
                TB = CblasConjTrans;
            }
            else
            {
                TB = CblasNoTrans;
            }

            if ( &A != &C )
            {
                cblas_cgemm(CblasRowMajor, TA, TB, M, N, K, &alpha, pA, lda, pB, ldb, &beta, pC, ldc);
            }
            else
            {
                hoMatrix<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_cgemm(CblasRowMajor, TransA, TransB, M, N, K, &alpha, pATmp, lda, pB, ldb, &beta, pC, ldc);
            }
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            if ( transA )
            {
                TA = CblasConjTrans;
            }
            else
            {
                TA = CblasNoTrans;
            }

            if ( transB )
            {
                TB = CblasConjTrans;
            }
            else
            {
                TB = CblasNoTrans;
            }

            if ( &A != &C )
            {
                cblas_zgemm(CblasRowMajor, TA, TB, M, N, K, &alpha, pA, lda, pB, ldb, &beta, pC, ldc);
            }
            else
            {
                hoMatrix<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_zgemm(CblasRowMajor, TransA, TransB, M, N, K, &alpha, pATmp, lda, pB, ldb, &beta, pC, ldc);
            }
        }
        else
        {
            GADGET_ERROR_MSG("GeneralMatrixProduct_gemm : unsupported type " << typeid(T));
            return false;
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoMatrix<T>& C, const hoMatrix<T>& A, bool transA, const hoMatrix<T>& B, bool transB) ...");
        return false;
    }
    return true;
}

template<typename T> 
bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<T>& A, char uplo)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return true;
        GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

        int info;
        lapack_int n = (lapack_int)(A.rows());
        T* pA = A.begin();
        lapack_int lda = (lapack_int)(A.cols());

        if ( typeid(T)==typeid(float) )
        {
            info = LAPACKE_spotrf(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<float*>(pA), lda);
        }
        else if ( typeid(T)==typeid(double) )
        {
            info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<double*>(pA), lda);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            info = LAPACKE_cpotrf(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<MKL_Complex8*>(pA), lda);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            info = LAPACKE_zpotrf(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<MKL_Complex16*>(pA), lda);
        }
        else
        {
            GADGET_ERROR_MSG("CholeskyHermitianPositiveDefinite_potrf : unsupported type " << typeid(T));
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(info==0);

        if ( uplo == 'U' )
        {
            GADGET_CHECK_RETURN_FALSE(A.lowerTri(0));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(A.upperTri(0));
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in CholeskyHermitianPositiveDefinite_potrf(hoMatrix<T>& A, char uplo) ...");
        return false;
    }
    return true;
}

template<typename T> 
bool EigenAnalysis_syev_heev(hoMatrix<T>& A, hoMatrix<typename realType<T>::Type>& eigenValue)
{
    try
    {
        int M = (int)A.rows();
        GADGET_CHECK_RETURN_FALSE(A.cols() == M));

        if ( (eigenValue.rows()!=M) || (eigenValue.cols()!=1) )
        {
            GADGET_CHECK_RETURN_FALSE(D.createMatrix(M, 1));
        }

        int info;
        char jobz = 'V';
        char uplo = 'L';
        T* pA = A.begin();
        typename realType<T>::Type* pEV = eigenValue.begin();

        if ( typeid(T)==typeid(float) )
        {
            info = LAPACKE_ssyev(LAPACK_ROW_MAJOR, jobz, uplo, M, reinterpret_cast<float*>(pA), M, reinterpret_cast<float*>(pEV));
        }
        else if ( typeid(T)==typeid(double) )
        {
            info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplo, M, reinterpret_cast<double*>(pA), M, reinterpret_cast<double*>(pEV));
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            info = LAPACKE_cheev(LAPACK_ROW_MAJOR, jobz, uplo, M, reinterpret_cast<MKL_Complex8*>(pA), M, reinterpret_cast<float*>(pEV));
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            info = LAPACKE_zheev(LAPACK_ROW_MAJOR, jobz, uplo, M, reinterpret_cast<MKL_Complex16*>(pA), M, reinterpret_cast<double*>(pEV));
        }
        else
        {
            GADGET_ERROR_MSG("EigenAnalysis_syev_heev : unsupported type " << typeid(T));
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(info==0);
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in EigenAnalysis_syev_heev(hoMatrix<T>& A, hoMatrix<typename realType<T>::Type>& eigenValue) ... ");
        return false;
    }
    return true;
}

template<typename T> 
bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<T>& A)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return true;
        GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

        int info;
        char uplo = 'L';
        lapack_int n = (lapack_int)A.rows();
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.cols();

        if ( typeid(T)==typeid(float) )
        {
            info = LAPACKE_spotrf(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<float*>(pA), lda);
            GADGET_CHECK_RETURN_FALSE(info==0);

            info = LAPACKE_spotri(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<float*>(pA), lda);
            GADGET_CHECK_RETURN_FALSE(info==0);
        }
        else if ( typeid(T)==typeid(double) )
        {
            info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<double*>(pA), lda);
            GADGET_CHECK_RETURN_FALSE(info==0);

            info = LAPACKE_dpotri(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<double*>(pA), lda);
            GADGET_CHECK_RETURN_FALSE(info==0);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            info = LAPACKE_cpotrf(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<MKL_Complex8*>(pA), lda);
            GADGET_CHECK_RETURN_FALSE(info==0);

            info = LAPACKE_cpotri(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<MKL_Complex8*>(pA), lda);
            GADGET_CHECK_RETURN_FALSE(info==0);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            info = LAPACKE_zpotrf(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<MKL_Complex16*>(pA), lda);
            GADGET_CHECK_RETURN_FALSE(info==0);

            info = LAPACKE_zpotri(LAPACK_ROW_MAJOR, uplo, n, reinterpret_cast<MKL_Complex16*>(pA), lda);
            GADGET_CHECK_RETURN_FALSE(info==0);
        }
        else
        {
            GADGET_ERROR_MSG("SymmetricHermitianPositiveDefiniteInverse_potri : unsupported type " << typeid(T));
            return false;
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<T>& A) ...");
        return false;
    }
    return true;
}

template<typename T> 
bool TriangularInverse_trtri(hoMatrix<T>& A, char uplo)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return true;
        GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

        int info;
        char diag = 'N';
        lapack_int n = (lapack_int)A.rows();
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.cols();

        if ( typeid(T)==typeid(float) )
        {
            info = LAPACKE_strtri(LAPACK_ROW_MAJOR, uplo, diag, n, reinterpret_cast<float*>(pA), lda);
        }
        else if ( typeid(T)==typeid(double) )
        {
            info = LAPACKE_dtrtri(LAPACK_ROW_MAJOR, uplo, diag, n, reinterpret_cast<double*>(pA), lda);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            info = LAPACKE_ctrtri(LAPACK_ROW_MAJOR, uplo, diag, n, reinterpret_cast<MKL_Complex8*>(pA), lda);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            info = LAPACKE_ztrtri(LAPACK_ROW_MAJOR, uplo, diag, n, reinterpret_cast<MKL_Complex16*>(pA), lda);
        }
        else
        {
            GADGET_ERROR_MSG("TriangularInverse_trtri : unsupported type " << typeid(T));
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(info==0);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in TriangularInverse_trtri(hoMatrix<float>& A, char uplo) ...");
        return false;
    }
    return true;
}

template<typename T> 
bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<T>& A, hoMatrix<T>& b)
{
    try
    {
        if( A.get_number_of_elements()==0 ) return true;
        if( b.get_number_of_elements()==0 ) return true;
        GADGET_CHECK_RETURN_FALSE(A.rows()==b.rows());

        int info;
        char uplo = 'L';
        lapack_int n = (lapack_int)A.rows();
        lapack_int nrhs = (lapack_int)b.cols();
        T* pA = A.begin();
        lapack_int lda = (lapack_int)A.cols();
        T* pB = b.begin();
        lapack_int ldb = (lapack_int)b.cols();

        if ( typeid(T)==typeid(float) )
        {
            info = LAPACKE_sposv(LAPACK_ROW_MAJOR, uplo, n, nrhs, reinterpret_cast<float*>(pA), lda, reinterpret_cast<float*>(pB), ldb);
        }
        else if ( typeid(T)==typeid(double) )
        {
            info = LAPACKE_dposv(LAPACK_ROW_MAJOR, uplo, n, nrhs, reinterpret_cast<double*>(pA), lda, reinterpret_cast<double*>(pB), ldb);
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            info = LAPACKE_cposv(LAPACK_ROW_MAJOR, uplo, n, nrhs, reinterpret_cast<MKL_Complex8*>(pA), lda, reinterpret_cast<MKL_Complex8*>(pB), ldb);
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            info = LAPACKE_zposv(LAPACK_ROW_MAJOR, uplo, n, nrhs, reinterpret_cast<MKL_Complex16*>(pA), lda, reinterpret_cast<MKL_Complex16*>(pB), ldb);
        }
        else
        {
            GADGET_ERROR_MSG("SymmetricHermitianPositiveDefiniteLinearSystem_posv : unsupported type " << typeid(T));
            return false;
        }

        GADGET_CHECK_RETURN_FALSE(info==0);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<float>& A, hoMatrix<float>& b) ...");
        return false;
    }
    return true;
}

#endif // USE_MKL

}
