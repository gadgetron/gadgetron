#pragma once

#include "ho2DArray.h"
#include "complext.h"
#include <algorithm>
#include <iomanip>

namespace Gadgetron{

// the hoMatrix stores every column as the first dimension
// it has the column-wise storage
template <class T> class  hoMatrix : public ho2DArray<T>
{
public:

    typedef hoMatrix<T> Self;
    typedef ho2DArray<T> BaseClass;

    hoMatrix();
    hoMatrix(size_t rows, size_t cols);
    hoMatrix(size_t rows, size_t cols, T* data, bool delete_data_on_destruct = false);

    virtual ~hoMatrix();

    hoMatrix(const hoMatrix<T>& a);
    hoMatrix<T>& operator=(const hoMatrix& rhs);

    virtual bool createMatrix(size_t rows, size_t cols);
    virtual bool createMatrix(size_t rows, size_t cols, T* data, bool delete_data_on_destruct = false);

    T& operator()(size_t r , size_t c);
    const T& operator()(size_t r , size_t c) const;

    size_t rows() const;
    size_t cols() const;

    // assign the upper/lower triangle matrix as a fixed value
    bool upperTri(const T& v);
    bool lowerTri(const T& v);

    // copy upper triangle to the lower
    bool copyUpperTriToLower();
    bool copyLowerTriToUpper();

    // sum along row or col
    bool sumOverRow(hoNDArray<T>& res) const;
    bool sumOverCol(hoNDArray<T>& res) const;

    // get the sub matrix
    bool subMatrix(Self& res, size_t startR, size_t endR, size_t startC, size_t endC) const;

    // set the matrix to be identity
    bool setIdentity();

    // normalize the matrix, so the L2 norm of matrix is 1
    bool normalize();

    bool operator == (const Self& m) const;
    bool operator != (const Self& m) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;
    using BaseClass::accesser_;

};

/// for real matrix
template <class T> class hoMatrixReal : public hoMatrix<T>
{
public:

    typedef hoMatrixReal<T> Self;
    typedef hoMatrix<T> BaseClass;

    hoMatrixReal();
    hoMatrixReal(size_t rows, size_t cols);
    hoMatrixReal(size_t rows, size_t cols, T* data, bool delete_data_on_destruct = false);

    virtual ~hoMatrixReal();

    hoMatrixReal(const hoMatrixReal<T>& a);

    /// sort along the row direction (sort along the 1st dimension)
    bool sort_ascending_along_row();

    /// sort along the column direction (sort along the 2nd dimension)
    bool sort_ascending_along_column();

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;
    using BaseClass::accesser_;
};

}

#include "hoMatrix.hxx"
