
namespace Gadgetron
{

template <typename T>
hoMatrix<T>::hoMatrix() : BaseClass(1, 1)
{
}

template <typename T>
hoMatrix<T>::hoMatrix(size_t rows, size_t cols) : BaseClass(rows, cols)
{
    this->fill(T(0));
}

template <typename T>
hoMatrix<T>::hoMatrix(size_t rows, size_t cols, T* data, bool delete_data_on_destruct)
{
    std::vector<size_t> dim(2);
    dim[0] = rows;
    dim[1] = cols;
    this->create(dim,data,delete_data_on_destruct);
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
bool hoMatrix<T>::createMatrix(size_t rows, size_t cols)
{
    return this->createArray(rows, cols);
}

template <typename T>
bool hoMatrix<T>::createMatrix(size_t rows, size_t cols, T* data, bool delete_data_on_destruct)
{
    return this->createArray(rows, cols, data, delete_data_on_destruct);
}

template <typename T>
inline T& hoMatrix<T>::operator()(size_t r, size_t c)
{
    GADGET_DEBUG_CHECK_THROW(c>=0 && r>=0 && r<dimensions_[0] && c<dimensions_[1]);
    return accesser_[c][r];
}

template <typename T>
inline const T& hoMatrix<T>::operator()(size_t r, size_t c) const
{
    GADGET_DEBUG_CHECK_THROW(c>=0 && r>=0 && r<dimensions_[0] && c<dimensions_[1]);
    return accesser_[c][r];
}

template <typename T>
inline size_t hoMatrix<T>::rows() const
{
    if ( dimensions_.empty() ) return 0;
    return dimensions_[0];
}

template <typename T>
inline size_t hoMatrix<T>::cols() const
{
    if ( dimensions_.empty() ) return 0;
    return dimensions_[1];
}

template <typename T>
bool hoMatrix<T>::upperTri(const T& v)
{
    try
    {
        size_t r, c;
        for (r=0; r<dimensions_[0]; r++)
        {
            for (c=r+1; c<dimensions_[1]; c++)
            {
                (*this)(r, c) = v;
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoMatrix<T>::upperTri(const T& v) ... ");
        return false;
    }
    return true;
}

template <typename T>
bool hoMatrix<T>::lowerTri(const T& v)
{
    try
    {
        size_t r, c;
        for (c=0; c<dimensions_[1]; c++)
        {
            for (r=c+1; r<dimensions_[0]; r++)
            {
                (*this)(r, c) = v;
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoMatrix<T>::lowerTri(const T& v) ... ");
        return false;
    }
    return true;
}

template <typename T>
bool hoMatrix<T>::copyUpperTriToLower()
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(dimensions_[0]==dimensions_[1]);

        size_t r, c;
        for (r=0; r<dimensions_[0]; r++)
        {
            for (c=r+1; c<dimensions_[1]; c++)
            {
                (*this)(c, r)= (*this)(r, c);
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoMatrix<T>::copyUpperTriToLower() ... ");
        return false;
    }
    return true;
}

template <typename T>
bool hoMatrix<T>::copyLowerTriToUpper()
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(dimensions_[0]==dimensions_[1]);

        size_t r, c;
        for (c=0; c<dimensions_[1]; c++)
        {
            for (r=c+1; r<dimensions_[0]; r++)
            {
                (*this)(c, r)= (*this)(r, c);
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoMatrix<T>::copyUpperTriToLower() ... ");
        return false;
    }
    return true;
}

template <typename T>
bool hoMatrix<T>::sumOverRow(hoNDArray<T>& res) const
{
    try
    {
        size_t ROW = rows();
        size_t COL = cols();

        if ( res.get_number_of_elements() != ROW )
        {
            res.create(ROW);
        }

        T* pRes = res.begin();

        size_t r, c;

        for ( r=0; r<ROW; r++ )
        {
            pRes[r] = 0;
        }

        for ( c=0; c<COL; c++ )
        {
            for ( r=0; r<ROW; r++ )
            {
                // res(r) += (*this)(r, c);
                pRes[r] += this->data_[r+c*ROW];
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoMatrix<T>::sumOverRow(hoNDArray<T>& r) ... ");
        return false;
    }

    return true;
}

template <typename T>
bool hoMatrix<T>::sumOverCol(hoNDArray<T>& res) const
{
    try
    {
        size_t ROW = rows();
        size_t COL = cols();

        if ( res.get_number_of_elements() != COL )
        {
            res.create(COL);
        }

        T* pRes = res.begin();

        size_t r;
        long long c;

        for ( c=0; c<(long long)COL; c++ )
        {
            pRes[c] = 0;
        }

        //for ( r=0; r<ROW; r++ )
        //{
        //    for ( c=0; c<COL; c++ )
        //    {
        //        // res(c) += (*this)(r, c);
        //        pRes[c] += this->data_[r+c*ROW];
        //    }
        //}

        T* pCurr = NULL;
        T v(0);
        // #pragma omp parallel for default(none) private(c, r) shared(COL, ROW, pRes) if ( COL > 16 )
        for ( c=0; c<(long long)COL; c++ )
        {
            v = 0;
            pCurr = this->data_ + c*ROW;
            for ( r=0; r<ROW; r++ )
            {
                v += pCurr[r];
            }
            pRes[c] = v;
        }

        //size_t r, c;
        //for ( c=0; c<COL; c++ )
        //{
        //    T v = (*this)(0, c);
        //    for ( r=1; r<ROW; r++ )
        //    {
        //        v += (*this)(r, c);
        //        //v += this->data_[r+c*ROW];
        //    }
        //    res(c) = v;
        //}
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoMatrix<T>::sumOverCol(hoNDArray<T>& r) ... ");
        return false;
    }

    return true;
}

template <typename T>
bool hoMatrix<T>::subMatrix(Self& res, size_t startR, size_t endR, size_t startC, size_t endC) const
{
    try
    {
        size_t ROW = rows();
        size_t COL = cols();

        GADGET_CHECK_RETURN_FALSE(startR<ROW);
        GADGET_CHECK_RETURN_FALSE(startC<COL);
        GADGET_CHECK_RETURN_FALSE(endR<ROW);
        GADGET_CHECK_RETURN_FALSE(endC<COL);
        GADGET_CHECK_RETURN_FALSE(endR>=startR);
        GADGET_CHECK_RETURN_FALSE(endC>=startC);

        GADGET_CHECK_RETURN_FALSE(res.createMatrix(endR-startR+1, endC-startC+1));

        size_t r, c;
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
        GERROR_STREAM("Errors in hoMatrix<T>::subMatrix(Self& res, size_t startR, size_t endR, size_t startC, size_t endC) ... ");
        return false;
    }

    return true;
}

template <typename T>
bool hoMatrix<T>::setIdentity()
{
    try
    {
        size_t ROW = this->rows();
        size_t COL = this->cols();

        size_t N = ((ROW<COL) ? ROW : COL);

        this->fill(T(0));

        size_t r;
        for ( r=0; r<N; r++ )
        {
            (*this)(r, r) = T(1.0);
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoMatrix<T>::setIdentity() ... ");
        return false;
    }

    return true;
}

template <typename T>
bool hoMatrix<T>::normalize()
{
    try
    {
        T dist = std::abs(this->data_[0]);
        dist *= dist;

        unsigned int ii;
        for ( ii=1; ii<this->element_; ii++ )
        {
            T v = std::abs(this->data_[ii]);
            dist += v*v;
        }

        dist = std::sqrt(dist);

        if ( std::abs(dist) < DBL_EPSILON ) return false;

        for ( ii=0; ii<this->element_; ii++ )
        {
            this->data_[ii] /= dist;
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoMatrix<T>::normalize() ... ");
        return false;
    }

    return true;
}

template <typename T>
bool hoMatrix<T>::operator == (const Self& m) const
{
    GADGET_CHECK_RETURN_FALSE(this->dimensions_equal(m));
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

    os << "hoMatrix (row X col): " << this->rows() << " X " << this->cols() << " : " << std::string(typeid(T).name()) << endl;
    size_t r, c;
    for (r=0; r<dimensions_[0]; r++)
    {
        os << "r " << r << ":\t";
        for (c=0; c<dimensions_[1]; c++)
        {
            os << setprecision(10) << (*this)(r,c) << "\t";
        }
        os << endl;
    }
}

// --------------------------------------------------------------------------------------------------------

template <typename T>
hoMatrixReal<T>::hoMatrixReal() : BaseClass()
{
}

template <typename T>
hoMatrixReal<T>::hoMatrixReal(size_t rows, size_t cols) : BaseClass(rows, cols)
{
}

template <typename T>
hoMatrixReal<T>::hoMatrixReal(size_t rows, size_t cols, T* data, bool delete_data_on_destruct) : BaseClass(rows, cols, delete_data_on_destruct)
{
}

template <typename T>
hoMatrixReal<T>::~hoMatrixReal()
{
}

template <typename T>
hoMatrixReal<T>::hoMatrixReal(const hoMatrixReal<T>& a) : BaseClass(a)
{
}

template <typename T>
bool hoMatrixReal<T>::sort_ascending_along_row()
{
    try
    {
        size_t R = this->rows();
        size_t C = this->cols();

        size_t col;
        for(col=0; col<C; col++)
        {
            std::sort(data_+col*R, data_+(col+1)*R);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in hoMatrixReal<T>::sort_ascending_along_row() ... ");
        return false;
    }
    return true;
}

template <typename T>
bool hoMatrixReal<T>::sort_ascending_along_column()
{
    try
    {
        size_t R = this->rows();
        size_t C = this->cols();

        std::vector<T> buf(C);

        size_t col, row;
        for(row=0; row<R; row++)
        {
            for(col=0; col<C; col++)
            {
                buf[col] = data_[row + col*R];
            }

            std::sort(buf.begin(), buf.end());

            for(col=0; col<C; col++)
            {
                data_[row + col*R] = buf[col];
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in hoMatrixReal<T>::sort_ascending_along_column() ... ");
        return false;
    }
    return true;
}

// --------------------------------------------------------------------------------------------------------

template <typename T>
bool copyL2U(hoMatrix<T>& A)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(A.rows()==A.cols());

        size_t R = A.rows();

        size_t row, col;
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
        GERROR_STREAM("Errors in copyL2U(hoMatrix<T>& A) ... ");
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

        size_t R = A.rows();
        size_t row, col;

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
        GERROR_STREAM("Errors in copyL2U(hoMatrix<T>& A, bool conj) ... ");
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

        size_t R = A.rows();
        size_t C = A.cols();

        size_t row, col;
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
        GERROR_STREAM("Errors in copyU2L(hoMatrix<T>& A) ... ");
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

        size_t R = A.rows();
        size_t C = A.cols();

        size_t row, col;

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
        GERROR_STREAM("Errors in copyU2L(hoMatrix<T>& A, bool conj) ... ");
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
        #pragma omp parallel for default(none) private(r, c) shared(A, AT)
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
        GERROR_STREAM("Errors in trans(const hoMatrix<T>& A, hoMatrix<T>& AT) ... ");
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
        #pragma omp parallel for default(none) private(r, c) shared(A, AH)
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
        GERROR_STREAM("Errors in conjugatetrans(const hoMatrix<T>& A, hoMatrix<T>& AH) ... ");
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

// C = A*B
bool GeneralMatrixProduct(hoNDArray<float>& C, const hoNDArray<float>& A, bool transA, const hoNDArray<float>& B, bool transB);
bool GeneralMatrixProduct(hoNDArray<double>& C, const hoNDArray<double>& A, bool transA, const hoNDArray<double>& B, bool transB);
bool GeneralMatrixProduct(hoNDArray< std::complex<float> >& C, const hoNDArray< std::complex<float> >& A, bool transA, const hoNDArray< std::complex<float> >& B, bool transB);
bool GeneralMatrixProduct(hoNDArray< std::complex<double> >& C, const hoNDArray< std::complex<double> >& A, bool transA, const hoNDArray< std::complex<double> >& B, bool transB);

template<typename T>
bool GeneralMatrixProduct(hoMatrix<T>& C, const hoMatrix<T>& A, bool transA, const hoMatrix<T>& B, bool transB)
{
    try
    {
        hoNDArray<T> mC(C.get_dimensions(), C.begin(), false);
        hoNDArray<T> mA(A.get_dimensions(), const_cast<T*>(A.begin()), false);
        hoNDArray<T> mB(B.get_dimensions(), const_cast<T*>(B.begin()), false);

        Gadgetron::GeneralMatrixProduct(mC, mA, transA, mB, transB);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GeneralMatrixProduct(hoMatrix<T>& C, const hoMatrix<T>& A, bool transA, const hoMatrix<T>& B, bool transB) ...");
        return false;
    }
    return true;
}

}
