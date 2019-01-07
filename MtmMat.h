#ifndef EX3_MTMMAT_H
#define EX3_MTMMAT_H


#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "MtmVec.h"

using std::size_t;

namespace MtmMath {

    template <typename T>
    class MtmMat : MtmVec<MtmVec<T>> {
    protected:
        const Dimensions dimensions; //Roi
    public:
        /*
         * Matrix constructor, dim_t is the dimension of the matrix and val is the initial value for the matrix elements
         */
        MtmMat<T>(Dimensions dim_t, const T& val=T()); //implemented

        //copy constructor
        MtmMat<T>(const MtmMat<T>& original) = default;

        //assignment
        MtmMat<T>& operator=(const MtmMat& original) = default;

        //destructor (Note: MtmVec matrix will destroy its own sub-lists)
        virtual ~MtmMat<T>() = default;

        //Roi :
        //operator[] inherited from MtmVec

        /* Matrix multiplication
         * creates new matrix
         * technique: when multiplying matrices A*B
         * row i of new matrix equals sum of scalar products of jth element of
         * row i of A * row j of B
         */
        friend MtmMat<T> operator*(const MtmMat<T> matrix1,
                const MtmMat<T> matrix2);

        //Matrix addition
        friend MtmMat<T> operator+(const MtmMat<T> matrix1,
                const MtmMat<T> matrix2);


    /*
     * Function that get function object f and uses it's () operator on each element in the matrix columns.
     * It outputs a vector in the size of the matrix columns where each element is the final output
     * by the function object's * operator
     */
        template <typename Func>
        MtmVec<T> matFunc(Func& f) const;

        /*
         * resizes a matrix to dimension dim, new elements gets the value val.
         */
        virtual void resize(Dimensions dim, const T& val=T());

        /*
         * reshapes matrix so linear elements value are the same without changing num of elements.
         */
        virtual void reshape(Dimensions newDim);

        /*
         * Performs transpose operation on matrix
         */
        virtual void transpose();

    };
}

//implementations begin here
using namespace MtmMath;

template<typename T>
MtmMat<T>::MtmMat(Dimensions dim_t, const T &val) : dimensions(dim_t),
    MtmVec<MtmVec<T>>(dim_t.getRow(), MtmVec<T>(dim_t.getCol(), val)) {
        //TODO: make sure vectors are horizontal!  Perhaps it would be prudent to transpose them here, for instance...
    }

template <typename T>
MtmMat<T> operator*(const MtmMat<T> matrix1, const MtmMat<T> matrix2){
    const size_t n = matrix1.dimensions.getCol();
    if (n != matrix2.dimensions.getRows()){
        throw MtmExceptions::DimensionMismatch(matrix1.dimensions,
                matrix2.dimensions);
    }

    const size_t numRows = matrix1.dimensions.getRow();
    const size_t numCols = matrix2.dimensions.getCol();
    MtmMat<T> answer(Dimensions(numRows, numCols), 0);
    for (int answerRow=0; answerRow < numRows; ++answerRow){
        for (int runner=0; runner<n; ++runner){
            answer[answerRow] = answer[answerRow]
                    + (matrix1[answerRow][runner] * matrix2[runner]);
        }
    }

    return answer;
}

template <typename T>
MtmMat<T> operator+(const MtmMat<T> matrix1, const MtmMat<T> matrix2){
    const size_t numRows = matrix1.dimensions.getRows();
    const size_t numCols = matrix1.dimensions.getCols();
    if (matrix2.dimensions.getRows() != numRows
        || matrix2.dimensions.getCols() != numCols){

        throw MtmExceptions::DimensionMismatch(matrix1.dimensions,
                matrix2.dimensions);
    }

    MtmMat<T> answer(Dimensions(numRows, numCols), 0);
    for (int row=0; row<numRows; ++row){
        answer[row] = matrix1[row] + matrix2[row]; //vector addition
    }
    return answer;
}


#endif //EX3_MTMMAT_H
