
namespace Gadgetron
{

template <typename T> 
hoMatrix<T>::hoMatrix() : BaseClass(1, 1)
{
}

template <typename T> 
hoMatrix<T>::hoMatrix(unsigned long long rows, unsigned long long cols) : BaseClass(cols, rows)
{
    this->fill(T(0));
}

template <typename T> 
hoMatrix<T>::hoMatrix(unsigned long long rows, unsigned long long cols, T* data, bool delete_data_on_destruct)
{
    std::vector<unsigned long long> dim(2);
    dim[0] = cols;
    dim[1] = rows;
    this->create(&dim,data,delete_data_on_destruct);
    GADGET_CHECK_THROW(this->init_accesser());
}

template <typename T> 
hoMatrix<T>::~hoMatrix()
{

}

template <typename T> 
hoMatrix<T>::hoMatrix(const hoMatrix<T>& a) : BaseClass(a)
{
}

template <typename T> 
hoMatrix<T>& hoMatrix<T>::operator=(const hoMatrix& rhs)
{
    if ( this == &rhs ) return *this;
    BaseClass::operator=(rhs);
    return *this;
}

template <typename T> 
bool hoMatrix<T>::createMatrix(unsigned long long rows, unsigned long long cols)
{
    return this->createArray(cols, rows);
}

template <typename T> 
bool hoMatrix<T>::createMatrix(unsigned long long rows, unsigned long long cols, T* data, bool delete_data_on_destruct)
{
    return this->createArray(cols, rows, data, delete_data_on_destruct);
}

template <typename T> 
inline T& hoMatrix<T>::operator()(long long r , long long c)
{
    GADGET_DEBUG_CHECK_THROW(c>=0 && r>=0 && c<(*dimensions_)[0] && r<(*dimensions_)[1]);
    return accesser_[r][c];
}

template <typename T> 
inline const T& hoMatrix<T>::operator()(long long r , long long c) const
{
    GADGET_DEBUG_CHECK_THROW(c>=0 && r>=0 && c<(*dimensions_)[0] && r<(*dimensions_)[1]);
    return accesser_[r][c];
}

template <typename T> 
inline unsigned long long hoMatrix<T>::rows() const
{
    if ( dimensions_->empty() ) return 0;
    return (*dimensions_)[1];
}

template <typename T> 
inline unsigned long long hoMatrix<T>::cols() const
{
    if ( dimensions_->empty() ) return 0;
    return (*dimensions_)[0];
}

template <typename T> 
bool hoMatrix<T>::upperTri(const T& v)
{
    try
    {
        unsigned long long r, c;
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
        unsigned long long r, c;
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
bool hoMatrix<T>::sumOverRow(hoNDArray<T>& res) const
{
    try
    {
        unsigned long long ROW = rows();
        unsigned long long COL = cols();

        if ( res.get_number_of_elements() != ROW )
        {
            res.create(ROW);
        }

        unsigned long long r, c;
        for ( r=0; r<ROW; r++ )
        {
            for ( c=0; c<COL; c++ )
            {
                res(r) += (*this)(r, c);
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in hoMatrix<T>::sumOverRow(hoNDArray<T>& r) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool hoMatrix<T>::sumOverCol(hoNDArray<T>& res) const
{
    try
    {
        unsigned long long ROW = rows();
        unsigned long long COL = cols();

        if ( res.get_number_of_elements() != COL )
        {
            res.create(COL);
        }

        long long c;

        //#pragma omp parallel for default(none) private(c) shared(ROW, COL, res)
        //for ( c=0; c<(long long)COL; c++ )
        //{
        //    for ( unsigned long long r=0; r<ROW; r++ )
        //    {
        //        res(c) += (*this)(r, c);
        //    }
        //}

        // #pragma omp parallel for default(none) private(c) shared(ROW, COL, res)
        for ( unsigned long long r=0; r<ROW; r++ )
        {
            for ( c=0; c<(long long)COL; c++ )
            {
                res(c) += (*this)(r, c);
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in hoMatrix<T>::sumOverCol(hoNDArray<T>& r) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool hoMatrix<T>::subMatrix(Self& res, unsigned long long startR, unsigned long long endR, unsigned long long startC, unsigned long long endC) const
{
    try
    {
        unsigned long long ROW = rows();
        unsigned long long COL = cols();

        GADGET_CHECK_RETURN_FALSE(startR>=0&&startR<ROW);
        GADGET_CHECK_RETURN_FALSE(startC>=0&&startC<COL);
        GADGET_CHECK_RETURN_FALSE(endR>=0&&endR<ROW);
        GADGET_CHECK_RETURN_FALSE(endC>=0&&endC<COL);
        GADGET_CHECK_RETURN_FALSE(endR>=startR);
        GADGET_CHECK_RETURN_FALSE(endC>=startC);

        GADGET_CHECK_RETURN_FALSE(res.createMatrix(endR-startR+1, endC-startC+1));

        unsigned long long r, c;
        for ( r=startR; r<=endR; r++ )
        {
            for ( c=startC; c<=endC; c++ )
            {
                res(r-startR, c-startC) = (*this)(r, c);
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in hoMatrix<T>::subMatrix(Self& res, unsigned long long startR, unsigned long long endR, unsigned long long startC, unsigned long long endC) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool hoMatrix<T>::operator == (const Self& m) const
{
    GADGET_CHECK_RETURN_FALSE(this->dimensions_equal(&m));
    for ( unsigned long long i=0; i<elements_; i++ )
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
    unsigned long long r, c;
    for (r=0; r<(*dimensions_)[1]; r++) 
    {
        os << "r " << r << ":\t";
        for (c=0; c<(*dimensions_)[0]; c++)
        {
            os << setprecision(10) << (*this)(r,c) << "\t";
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

        unsigned long long R = A.rows();
        unsigned long long C = A.cols();

        unsigned long long row, col;
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

        unsigned long long R = A.rows();
        unsigned long long row, col;

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

        unsigned long long R = A.rows();
        unsigned long long C = A.cols();

        unsigned long long row, col;
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

        unsigned long long R = A.rows();
        unsigned long long C = A.cols();

        unsigned long long row, col;

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

        if ( AT.rows()!=A.cols() || AT.cols()!=A.rows() )
        {
            AT.createMatrix(A.cols(), A.rows());
        }

        long long r, c;
        #ifdef GCC_OLD_FLAG
            #pragma omp parallel for default(none) private(r, c)
        #else
            #pragma omp parallel for default(none) private(r, c) shared(A, AT)
        #endif
        for ( c=0; c<(long long)A.cols(); c++ )
        {
            for ( r=0; r<(long long)A.rows(); r++ )
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

        if ( AH.rows()!=A.cols() || AH.cols()!=A.rows() )
        {
            AH.createMatrix(A.cols(), A.rows());
        }

        long long r, c;
        #ifdef GCC_OLD_FLAG
            #pragma omp parallel for default(none) private(r, c)
        #else
            #pragma omp parallel for default(none) private(r, c) shared(A, AH)
        #endif
        for ( c=0; c<(long long)A.cols(); c++ )
        {
            for ( r=0; r<(long long)A.rows(); r++ )
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

inline bool conjugatetrans(const hoMatrix<float>& A, hoMatrix<float>& AH)
{
    return trans(A, AH);
}

inline bool conjugatetrans(const hoMatrix<double>& A, hoMatrix<double>& AH)
{
    return trans(A, AH);
}

// following matrix computation calls MKL functions
#ifdef USE_MKL

template<typename T> 
bool GeneralMatrixProduct_gemm(hoNDArray<T>& C, 
                            const hoNDArray<T>& A, bool transA, 
                            const hoNDArray<T>& B, bool transB)
{
    try
    {
        CBLAS_TRANSPOSE TA, TB;

        MKL_INT lda = A.get_size(0);
        MKL_INT ldb = B.get_size(0);
        const T* pA = A.begin(); 
        const T* pB = B.begin(); 

        MKL_INT M = A.get_size(1);
        MKL_INT K = A.get_size(0);
        if ( transA )
        { 
            M = A.get_size(0);
            K = A.get_size(1);
        }

        MKL_INT N = B.get_size(0);
        MKL_INT K2 = B.get_size(1);
        if ( transB )
        { 
            N = B.get_size(1);
            K2 = B.get_size(0);
        }

        GADGET_CHECK_RETURN_FALSE(K==K2);
        if ( (C.get_size(1)!=M) || (C.get_size(0)!=N) )
        {
            C.create(N, M);
        }

        T* pC = C.begin();
        MKL_INT ldc = C.get_size(0);

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
                hoNDArray<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_sgemm(CblasRowMajor, TA, TB, M, N, K, 1, reinterpret_cast<const float*>(pATmp), lda, reinterpret_cast<const float*>(pB), ldb, 0, reinterpret_cast<float*>(pC), ldc);
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
                hoNDArray<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_dgemm(CblasRowMajor, TA, TB, M, N, K, 1, reinterpret_cast<const double*>(pATmp), lda, reinterpret_cast<const double*>(pB), ldb, 0, reinterpret_cast<double*>(pC), ldc);
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
                hoNDArray<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_cgemm(CblasRowMajor, TA, TB, M, N, K, &alpha, pATmp, lda, pB, ldb, &beta, pC, ldc);
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
                hoNDArray<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_zgemm(CblasRowMajor, TA, TB, M, N, K, &alpha, pATmp, lda, pB, ldb, &beta, pC, ldc);
            }
        }
        else
        {
            GADGET_ERROR_MSG("GeneralMatrixProduct_gemm : unsupported type " << typeid(T).name() );
            return false;
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GeneralMatrixProduct_gemm(hoNDArray<T>& C, const hoNDArray<T>& A, bool transA, const hoMatrix<T>& B, bool transB) ...");
        return false;
    }
    return true;
}

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
                hoNDArray<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_sgemm(CblasRowMajor, TA, TB, M, N, K, 1, reinterpret_cast<const float*>(pATmp), lda, reinterpret_cast<const float*>(pB), ldb, 0, reinterpret_cast<float*>(pC), ldc);
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
                hoNDArray<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_dgemm(CblasRowMajor, TA, TB, M, N, K, 1, reinterpret_cast<const double*>(pATmp), lda, reinterpret_cast<const double*>(pB), ldb, 0, reinterpret_cast<double*>(pC), ldc);
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
                hoNDArray<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_cgemm(CblasRowMajor, TA, TB, M, N, K, &alpha, pATmp, lda, pB, ldb, &beta, pC, ldc);
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
                hoNDArray<T> aTmp(A);
                T* pATmp = aTmp.begin();
                cblas_zgemm(CblasRowMajor, TA, TB, M, N, K, &alpha, pATmp, lda, pB, ldb, &beta, pC, ldc);
            }
        }
        else
        {
            GADGET_ERROR_MSG("GeneralMatrixProduct_gemm : unsupported type " << typeid(T).name() );
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

        MKL_INT info;
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
            GADGET_ERROR_MSG("CholeskyHermitianPositiveDefinite_potrf : unsupported type " << typeid(T).name());
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
        long long M = (long long)A.rows();
        GADGET_CHECK_RETURN_FALSE(A.cols() == M);

        if ( (eigenValue.rows()!=M) || (eigenValue.cols()!=1) )
        {
            GADGET_CHECK_RETURN_FALSE(eigenValue.createMatrix(M, 1));
        }

        MKL_INT info;
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
            GADGET_ERROR_MSG("EigenAnalysis_syev_heev : unsupported type " << typeid(T).name());
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
bool EigenAnalysis_syev_heev2(hoMatrix<T>& A, hoMatrix<T>& eigenValue)
{
    try
    {
        long long M = (long long)A.rows();
        GADGET_CHECK_RETURN_FALSE(A.cols() == M);

        if ( (eigenValue.rows()!=M) || (eigenValue.cols()!=1) )
        {
            GADGET_CHECK_RETURN_FALSE(eigenValue.createMatrix(M, 1));
        }

        hoMatrix<typename realType<T>::Type> D(M, 1);
        GADGET_CHECK_RETURN_FALSE(EigenAnalysis_syev_heev(A, D));
        GADGET_CHECK_RETURN_FALSE(eigenValue.copyFrom(D));
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in EigenAnalysis_syev_heev2(hoMatrix<T>& A, hoMatrix<T>& eigenValue) ... ");
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

        MKL_INT info;
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
            GADGET_ERROR_MSG("SymmetricHermitianPositiveDefiniteInverse_potri : unsupported type " << typeid(T).name());
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

        MKL_INT info;
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
            GADGET_ERROR_MSG("TriangularInverse_trtri : unsupported type " << typeid(T).name());
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

        MKL_INT info;
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
            GADGET_ERROR_MSG("SymmetricHermitianPositiveDefiniteLinearSystem_posv : unsupported type " << typeid(T).name());
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

template<typename T> 
bool SolveLinearSystem_Tikhonov(hoMatrix<T>& A, hoMatrix<T>& b, hoMatrix<T>& x, double lamda)
{
    GADGET_CHECK_RETURN_FALSE(b.rows()==A.rows());

    hoMatrix<T> AHA(A.cols(), A.cols());
    GADGET_CHECK_RETURN_FALSE(GeneralMatrixProduct_gemm(AHA, A, true, A, false));

    GADGET_CHECK_RETURN_FALSE(x.createMatrix(A.cols(), b.cols()));
    GADGET_CHECK_RETURN_FALSE(GeneralMatrixProduct_gemm(x, A, true, b, false));

    // apply the Tikhonov regularization
    // Ideally, we shall apply the regularization is lamda*maxEigenValue
    // However, computing the maximal eigenvalue is computational intensive
    // A natural alternative is to use the trace of AHA matrix, which is the sum of all eigen values
    // Since all eigen values are positive, the lamda*maxEigenValue is only ~10-20% different from lamda*sum(all eigenValues)
    // for more information, refer to:
    // Tikhonov A.N., Goncharsky A.V., Stepanov V.V., Yagola A.G., 1995, 
    // Numerical Methods for the Solution of Ill-Posed Problems, Kluwer Academic Publishers.

    unsigned long long col = AHA.cols();
    unsigned long long c;

    double trA = std::abs(AHA(0, 0));
    for ( c=1; c<col; c++ )
    {
        trA += std::abs(AHA(c, c));
    }

    double value = trA*lamda/col;
    for ( c=0; c<col; c++ )
    {
        AHA(c,c) = std::abs(AHA(c, c)) + value;
    }

    GADGET_CHECK_RETURN_FALSE(SymmetricHermitianPositiveDefiniteLinearSystem_posv(AHA, x));

    return true;
}

#endif // USE_MKL

//    //
//    // Instantiation
//    //
//
//    template class hoMatrix<float>;
//    template class hoMatrix<double>;
//    template class hoMatrix<GT_Complex8>;
//    template class hoMatrix<GT_Complex16>;
//
//    template EXPORTCPUCORE bool copyL2U(hoMatrix<float>& A);
//    template EXPORTCPUCORE bool copyL2U(hoMatrix<double>& A);
//    template EXPORTCPUCORE bool copyL2U(hoMatrix<GT_Complex8>& A, bool conj);
//    template EXPORTCPUCORE bool copyL2U(hoMatrix<GT_Complex16>& A, bool conj);
//
//    template EXPORTCPUCORE bool copyU2L(hoMatrix<float>& A);
//    template EXPORTCPUCORE bool copyU2L(hoMatrix<double>& A);
//    template EXPORTCPUCORE bool copyU2L(hoMatrix<GT_Complex8>& A, bool conj);
//    template EXPORTCPUCORE bool copyU2L(hoMatrix<GT_Complex16>& A, bool conj);
//
//    template EXPORTCPUCORE bool trans(const hoMatrix<float>& A, hoMatrix<float>& AT);
//    template EXPORTCPUCORE bool trans(const hoMatrix<double>& A, hoMatrix<double>& AT);
//    template EXPORTCPUCORE bool trans(const hoMatrix<GT_Complex8>& A, hoMatrix<GT_Complex8>& AT);
//    template EXPORTCPUCORE bool trans(const hoMatrix<GT_Complex16>& A, hoMatrix<GT_Complex16>& AT);
//
//    template EXPORTCPUCORE bool conjugatetrans(const hoMatrix<GT_Complex8>& A, hoMatrix<GT_Complex8>& AH);
//    template EXPORTCPUCORE bool conjugatetrans(const hoMatrix<GT_Complex16>& A, hoMatrix<GT_Complex16>& AH);
//
//#ifdef USE_MKL
//
//    template<> EXPORTCPUCORE bool GeneralMatrixProduct_gemm(hoMatrix<float>& C, const hoMatrix<float>& A, bool transA, const hoMatrix<float>& B, bool transB);
//    template<> EXPORTCPUCORE bool GeneralMatrixProduct_gemm(hoMatrix<double>& C, const hoMatrix<double>& A, bool transA, const hoMatrix<double>& B, bool transB);
//    template<> EXPORTCPUCORE bool GeneralMatrixProduct_gemm(hoMatrix<GT_Complex8>& C, const hoMatrix<GT_Complex8>& A, bool transA, const hoMatrix<GT_Complex8>& B, bool transB);
//    template<> EXPORTCPUCORE bool GeneralMatrixProduct_gemm(hoMatrix<GT_Complex16>& C, const hoMatrix<GT_Complex16>& A, bool transA, const hoMatrix<GT_Complex16>& B, bool transB);
//
//    template<> EXPORTCPUCORE bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<float>& A, char uplo);
//    template<> EXPORTCPUCORE bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<double>& A, char uplo);
//    template<> EXPORTCPUCORE bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<GT_Complex8>& A, char uplo);
//    template<> EXPORTCPUCORE bool CholeskyHermitianPositiveDefinite_potrf(hoMatrix<GT_Complex16>& A, char uplo);
//
//    template<> EXPORTCPUCORE bool EigenAnalysis_syev_heev(hoMatrix<float>& A, hoMatrix<float>& eigenValue);
//    template<> EXPORTCPUCORE bool EigenAnalysis_syev_heev(hoMatrix<double>& A, hoMatrix<double>& eigenValue);
//    template<> EXPORTCPUCORE bool EigenAnalysis_syev_heev(hoMatrix<GT_Complex8>& A, hoMatrix<float>& eigenValue);
//    template<> EXPORTCPUCORE bool EigenAnalysis_syev_heev(hoMatrix<GT_Complex16>& A, hoMatrix<double>& eigenValue);
//
//    template<> EXPORTCPUCORE bool EigenAnalysis_syev_heev2(hoMatrix<float>& A, hoMatrix<float>& eigenValue);
//    template<> EXPORTCPUCORE bool EigenAnalysis_syev_heev2(hoMatrix<double>& A, hoMatrix<double>& eigenValue);
//    template<> EXPORTCPUCORE bool EigenAnalysis_syev_heev2(hoMatrix<GT_Complex8>& A, hoMatrix<GT_Complex8>& eigenValue);
//    template<> EXPORTCPUCORE bool EigenAnalysis_syev_heev2(hoMatrix<GT_Complex16>& A, hoMatrix<GT_Complex16>& eigenValue);
//
//    template<> EXPORTCPUCORE bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<float>& A);
//    template<> EXPORTCPUCORE bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<double>& A);
//    template<> EXPORTCPUCORE bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<GT_Complex8>& A);
//    template<> EXPORTCPUCORE bool SymmetricHermitianPositiveDefiniteInverse_potri(hoMatrix<GT_Complex16>& A);
//
//    template<> EXPORTCPUCORE bool TriangularInverse_trtri(hoMatrix<float>& A, char uplo);
//    template<> EXPORTCPUCORE bool TriangularInverse_trtri(hoMatrix<double>& A, char uplo);
//    template<> EXPORTCPUCORE bool TriangularInverse_trtri(hoMatrix<GT_Complex8>& A, char uplo);
//    template<> EXPORTCPUCORE bool TriangularInverse_trtri(hoMatrix<GT_Complex16>& A, char uplo);
//
//    template<> EXPORTCPUCORE bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<float>& A, hoMatrix<float>& b);
//    template<> EXPORTCPUCORE bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<double>& A, hoMatrix<double>& b);
//    template<> EXPORTCPUCORE bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<GT_Complex8>& A, hoMatrix<GT_Complex8>& b);
//    template<> EXPORTCPUCORE bool SymmetricHermitianPositiveDefiniteLinearSystem_posv(hoMatrix<GT_Complex16>& A, hoMatrix<GT_Complex16>& b);
//
//    template<> EXPORTCPUCORE bool SolveLinearSystem_Tikhonov(hoMatrix<float>& A, hoMatrix<float>& b, hoMatrix<float>& x, double lamda);
//    template<> EXPORTCPUCORE bool SolveLinearSystem_Tikhonov(hoMatrix<double>& A, hoMatrix<double>& b, hoMatrix<double>& x, double lamda);
//    template<> EXPORTCPUCORE bool SolveLinearSystem_Tikhonov(hoMatrix<GT_Complex8>& A, hoMatrix<GT_Complex8>& b, hoMatrix<GT_Complex8>& x, double lamda);
//    template<> EXPORTCPUCORE bool SolveLinearSystem_Tikhonov(hoMatrix<GT_Complex16>& A, hoMatrix<GT_Complex16>& b, hoMatrix<GT_Complex16>& x, double lamda);
//
//#endif // USE_MKL

}
